import os
import time
_t0_import = time.perf_counter()
import numpy as np
os.environ["JAX_ENABLE_X64"] = "1"     # MUST precede any jax import
os.environ["JAX_PLATFORMS"] = os.environ.get("JAX_PLATFORMS", "cpu")  # MUST precede any jax import; honors the driver script's export, else forces cpu

try:
    from kessler_perf_log import (log_call as _log_perf_call, log_event as _log_perf_event,
                                  set_enabled as _set_perf_enabled, flush_log as _flush_perf_log)
except ImportError:
    _log_perf_call = None
    _log_perf_event = None
    _set_perf_enabled = None
    _flush_perf_log = None

# This module is imported in the AtmosphereProcess constructor, outside every
# GPTL timer, so the import cost is logged here, under a_i:EAMxx::kessler::init.
_t0_jax_import = time.perf_counter()
import jax
import jax.numpy as jnp
if _log_perf_event is not None:
    _log_perf_event("a_i:EAMxx::kessler::init::py_import_jax", time.perf_counter() - _t0_jax_import)
jax.config.update("jax_enable_x64", True)

# JAX reports each compile stage (trace, lower to MLIR, XLA compile, persistent-cache
# read) through jax.monitoring. Log them per rank so a first call can be split into
# those stages. These fire only on a cache miss in this process, not on every call.
# Other debug output is set with environment variables (see run.py --jax-trace):
# JAX_LOG_COMPILES, JAX_EXPLAIN_CACHE_MISSES, JAX_DEBUG_LOG_MODULES.
if _log_perf_event is not None:
    def _jax_duration_listener(event, duration, **kwargs):
        fun = kwargs.get("fun_name")
        _log_perf_event("jax_monitor:" + event + (f"[{fun}]" if fun else ""), duration)

    def _jax_event_listener(event, **kwargs):   # counts only, e.g. cache hits/misses
        _log_perf_event("jax_monitor:" + event, 0.0)

    try:
        jax.monitoring.register_event_duration_secs_listener(_jax_duration_listener)
        jax.monitoring.register_event_listener(_jax_event_listener)
    except AttributeError:
        pass
# recommended: persistent compilation cache 
_cache_dir = os.environ.get(
    "JAX_COMPILATION_CACHE_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),".jax_cache"),
)
jax.config.update("jax_compilation_cache_dir", _cache_dir)
jax.config.update("jax_persistent_cache_min_compile_time_secs", 0)
jax.config.update("jax_persistent_cache_min_entry_size_bytes", -1)

# kessler_update is a low-arithmetic-intensity scheme: the host-facing
# bridge (contract 1) is dominated by the PCIe transfer -- measured 0.33x
# vs. serial Fortran, worse than not porting it at all. The device-resident
# bridge (contract 2) chains all three per-step phases on the GPU with no
# transfers between them and is the pattern this scheme strictly requires
# (measured 181x for the full 3-phase step). See USING_THE_BRIDGE.md and
# BRIDGE_CONTRACT2_notes.md in the laft-kessler-update project.
from bridge.kessler_update_init_bridge import kessler_update_init_bridge
from bridge.kessler_update_timestep_init_bridge import kessler_update_timestep_init_bridge_device
from bridge.kessler_update_run_bridge import kessler_update_run_bridge_device, to_device
from bridge.kessler_update_timestep_final_bridge import kessler_update_timestep_final_bridge_device

# kessler_jax.kessler_run_core expects arrays already in JAX (nz, ncol)
# layout -- kessler_run_bridge is what does the (ncol, nz) <-> (nz, ncol)
# layout conversion (in-jit, PATH C: pure H2D/D2H with the axis reversal
# fused into the jitted kernel by XLA), so EAMxx's native (ncol, nz)
# arrays go through the bridge, not kessler_jax.kessler_run directly.
from bridge.kessler_init_bridge import kessler_init_bridge
from bridge.kessler_run_bridge import kessler_run_bridge, kessler_run_bridge_device

# --- py_backend=device ---------------------------------------------------------
# With kessler::py_backend=device, EAMxx hands run()/update() NumPy arrays
# whose data pointers are DEVICE addresses (create_py_field<Device> in
# share/field/field_pyutils.hpp). They only carry the pointer, shape and
# strides: reading them as NumPy (jnp.asarray, theta[...] = ...) dereferences
# GPU memory on the host and segfaults. Upstream EAMxx's cld_fraction_cupy.py
# follows the same convention (it re-wraps them with cupy). Here, on a GPU,
# each one is wrapped in a __cuda_array_interface__ object, which jnp.asarray
# imports on the device without a copy, and the results are written back with
# device-to-device cudaMemcpy (libcudart ships with JAX's CUDA plugin, so no
# cupy is needed). The arrays must be C-contiguous: EAMxx passes contiguous
# copies of the strided tracers (KesslerMicrophysics::m_py_contig). On a CPU
# build device pointers are host pointers, so the host path below is used.
# Selected once by set_backend(), from KesslerMicrophysics::initialize_impl,
# before warmup(), which then compiles for the path that was selected.
_device_path = False
_cudart = None


def set_backend(py_backend):
    global _device_path, _cudart
    _device_path = (py_backend == "device" and jax.default_backend() == "gpu")
    if _device_path and _cudart is None:
        import ctypes
        _cudart = ctypes.CDLL("libcudart.so.12")
        _cudart.cudaMemcpy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int]
        _cudart.cudaMemcpy.restype = ctypes.c_int
        _cudart.cudaDeviceSynchronize.restype = ctypes.c_int
        _cudart.cudaMalloc.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_size_t]
        _cudart.cudaMalloc.restype = ctypes.c_int
        _cudart.cudaFree.argtypes = [ctypes.c_void_p]
        _cudart.cudaFree.restype = ctypes.c_int
    print(f"[kessler.py] py_backend={py_backend}, jax backend={jax.default_backend()}, "
          f"device path={'on' if _device_path else 'off'}", flush=True)


# --- kessler::py_run_variant ---------------------------------------------------
# How kessler_run's sub-cycle loop runs. "current": kessler_jax as translated
# (vmap of a per-column while_loop: every column does the slowest column's
# passes). Anything else: that variant of kessler_jax_experimental (see its
# README.md): auto (pallas on a GPU, map otherwise), pallas, split, map,
# split_map, vmap. All give the same results as "current" bit for bit in
# driver_kessler_py.cpp. Selected once by set_run_variant(), from
# KesslerMicrophysics::initialize_impl, before warmup().
_run_variant = "current"


def set_run_variant(variant):
    global _run_variant
    if variant != "current":
        from kessler_jax_experimental.kessler_run_variants import VARIANTS
        if variant not in VARIANTS:
            raise ValueError(f"kessler.py: unknown py_run_variant {variant!r}; "
                             f"one of current, {', '.join(VARIANTS)}")
    _run_variant = variant
    print(f"[kessler.py] py_run_variant={variant}, jax backend={jax.default_backend()}", flush=True)


def _kessler_run_device(**kw):
    """kessler_run_bridge_device, or the selected variant's (same arguments
    and outputs: device arrays in EAMxx's (ncol, nz) order)."""
    if _run_variant == "current":
        return kessler_run_bridge_device(**kw)
    from kessler_jax_experimental.kessler_run_variants import kessler_run_bridge_device_variant
    return kessler_run_bridge_device_variant(_run_variant, **kw)


def _kessler_run_host(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk, theta, qv, qc, qr,
                      precl, relhum, scheme_name, errmsg, errflg, lv, pref, rhoqr):
    """kessler_run_bridge (host NumPy arrays in and out), on the selected
    variant: the same steps as bridge/kessler_run_bridge.py (H2D, the
    device bridge, one batched D2H)."""
    if _run_variant == "current":
        return kessler_run_bridge(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, z, pk, theta, qv, qc, qr,
                                  precl, relhum, scheme_name, errmsg, errflg, lv, pref, rhoqr)
    out = _kessler_run_device(
        ncol=ncol, nz=nz, dt=dt, lyr_surf=lyr_surf, lyr_toa=lyr_toa, cpair=cpair, rair=rair,
        rho=to_device(rho), z=to_device(z), pk=to_device(pk), theta=to_device(theta),
        qv=to_device(qv), qc=to_device(qc), qr=to_device(qr), precl=to_device(precl),
        relhum=to_device(relhum), errflg=errflg, lv=lv, pref=pref, rhoqr=rhoqr)
    theta_o, qv_o, qc_o, qr_o, precl_o, relhum_o, errflg_o, lv, pref, rhoqr = out
    theta_h, qv_h, qc_h, qr_h, precl_h, relhum_h, errflg_h = jax.device_get(
        (theta_o, qv_o, qc_o, qr_o, precl_o, relhum_o, errflg_o))
    return (np.asarray(theta_h, dtype=np.float64), np.asarray(qv_h, dtype=np.float64),
            np.asarray(qc_h, dtype=np.float64), np.asarray(qr_h, dtype=np.float64),
            np.asarray(precl_h, dtype=np.float64), np.asarray(relhum_h, dtype=np.float64),
            scheme_name, errmsg, np.asarray(errflg_h).item(), lv, pref, rhoqr)


class _DeviceArray:
    """__cuda_array_interface__ view of an EAMxx device field. Reads only
    the NumPy array's metadata, never its data."""
    def __init__(self, a):
        if not a.flags.c_contiguous:
            raise ValueError("kessler.py: py_backend=device needs C-contiguous fields")
        self.__cuda_array_interface__ = {
            "shape": a.shape, "typestr": a.dtype.str, "strides": None, "version": 3,
            "data": (a.__array_interface__["data"][0], False)}


def _import(a):
    """An EAMxx device field as a JAX array on the GPU (no host copy)."""
    return jnp.asarray(_DeviceArray(a), dtype=jnp.float64)


def _copy_back(pairs):
    """Copy each (EAMxx device field, finished JAX array) pair on the GPU, then wait."""
    import ctypes
    for dst, src in pairs:
        cai = src.__cuda_array_interface__
        if dst.size != src.size or cai.get("strides") is not None:
            raise ValueError("kessler.py: JAX output does not match its EAMxx field")
        err = _cudart.cudaMemcpy(ctypes.c_void_p(dst.__array_interface__["data"][0]),
                                 ctypes.c_void_p(cai["data"][0]), dst.nbytes, 3)  # 3: device to device
        if err != 0:
            raise RuntimeError(f"kessler.py: cudaMemcpy failed (cudaError {err})")
    err = _cudart.cudaDeviceSynchronize()
    if err != 0:
        raise RuntimeError(f"kessler.py: cudaDeviceSynchronize failed (cudaError {err})")


def _device_copies(arrays, ptrs):
    """For warmup(): a GPU copy of each host array, handed back as the same
    kind of object EAMxx passes with py_backend=device (a NumPy array whose
    data pointer is a device address, never read on the host). The device
    pointers are appended to ptrs, for cudaFree."""
    import ctypes
    out = []
    for a in arrays:
        p = ctypes.c_void_p()
        err = _cudart.cudaMalloc(ctypes.byref(p), a.nbytes)
        if err != 0:
            raise RuntimeError(f"kessler.py: cudaMalloc failed (cudaError {err})")
        ptrs.append(p.value)
        err = _cudart.cudaMemcpy(p, ctypes.c_void_p(a.ctypes.data), a.nbytes, 1)  # 1: host to device
        if err != 0:
            raise RuntimeError(f"kessler.py: cudaMemcpy failed (cudaError {err})")
        # Wrapping the address reads nothing (like create_py_field<Device>)
        buf = (ctypes.c_char * a.nbytes).from_address(p.value)
        out.append(np.frombuffer(buf, dtype=a.dtype).reshape(a.shape))
    return out

# Setup a few global variables that we set values for with init() and then use in run()
latvap = None
pref = None
rhoqr = None
gravit = None
cpair = None
rair = None
errmsg = ""
errflg = 0
scheme_name = "kessler"

# latvap, P0, rhoqr, gravit, cpair, rair
def init(lv_in, pref_in, rhoqr_in, gravit_in, cpair_in, rair_in):
    global latvap, pref, rhoqr, gravit, cpair, rair
    global errmsg, errflg

    errmsg, errflg, latvap, pref, rhoqr = kessler_init_bridge(lv_in, pref_in, rhoqr_in, errmsg, errflg, latvap, pref, rhoqr)

    # cpair/rair are spatially-uniform physical constants (not per-level
    # MODULE state), so they are just stored here rather than threaded
    # through kessler_init_bridge's INOUT pattern.
    cpair = float(cpair_in)
    rair = float(rair_in)

    # kessler_update carries its own MODULE variable `gravit` (INOUT,
    # threaded from here through every kessler_update_timestep_final call
    # in update() below) -- separate from kessler's own MODULE vars
    # (latvap/pref/rhoqr) above.
    errmsg, errflg, gravit = kessler_update_init_bridge(gravit_in, "", 0, gravit)

# Compiles everything run() and update() call, so the first real step doesn't have
# to. Called once from KesslerMicrophysics::initialize_impl() after init() and
# set_backend() when kessler::py_warmup is true, inside GPTL's a_i:EAMxx::kessler::init.
#
# It calls run() and update() themselves on dummy arrays, so each jit sees exactly
# the arguments it will see in the model (shapes, dtypes, weak types, static args,
# and where the arrays live), and the real calls hit its in-memory cache. Calling
# jit(...).lower().compile() instead would not fill that cache. It calls them twice
# because the module variables (latvap, pref, rhoqr, gravit) come back from the
# first call as JAX arrays rather than Python floats, a different signature.
#
# On the device path the dummies are GPU copies passed as EAMxx passes its fields
# (_device_copies), so they go through _import and _copy_back like the real ones.
# That matters: arrays imported from device pointers are committed to the GPU,
# NumPy arrays given to jnp.asarray are not, and jit keys its compiled code on
# that too. Warmed with host arrays, the device path recompiled all four bridges
# and update()'s add/multiply on the first step, and kessler_run and
# timestep_final again on the second (driver_kessler_py.cpp --py-warmup).
#
# Dummy state: no rain or cloud (qc = qr = 0), so every column needs one
# sub-step. The compiled code depends only on shapes, not values.
def warmup(ncol, nz, lyr_surf, lyr_toa):
    global latvap, pref, rhoqr, gravit, errmsg, errflg, scheme_name

    t0 = time.perf_counter()
    # Backend start-up (XLA client, device discovery), otherwise paid by the first jnp call
    devices = jax.devices()
    if _log_perf_event is not None:
        _log_perf_event("a_i:EAMxx::kessler::init::py_backend_init", time.perf_counter() - t0)
    try:
        ncpu = len(os.sched_getaffinity(0))
    except AttributeError:
        ncpu = os.cpu_count()
    print(f"[kessler.py] warmup: ncol={ncol} nz={nz} devices={devices} cpus_in_affinity={ncpu} "
          f"XLA_FLAGS={os.environ.get('XLA_FLAGS', '')!r} cache={_cache_dir}", flush=True)

    saved = (latvap, pref, rhoqr, gravit, errmsg, errflg, scheme_name)
    device_ptrs = []
    if _set_perf_enabled is not None:
        _set_perf_enabled(False)
    try:
        dt = 1800.0   # traced (not static), so any Python float gives the same compile
        shape = (ncol, nz)
        z_mid = np.broadcast_to(500.0 * np.arange(nz, 0, -1, dtype=np.float64), shape).copy()  # level 0 = TOA
        for _ in range(2):
            run_args = [np.ones(shape),                    # rho
                        z_mid,                             # zm
                        np.ones(shape),                    # pk
                        np.full(shape, 300.0),             # theta
                        np.full(shape, 1.0e-3),            # qv
                        np.zeros(shape),                   # qc
                        np.zeros(shape),                   # qr
                        np.zeros(ncol),                    # precl
                        np.zeros(shape)]                   # relhum
            update_args = [z_mid,
                           np.ones(shape),                 # pk
                           np.full(shape, 300.0),          # theta
                           np.zeros(ncol),                 # phis
                           np.full(shape, 300.0),          # temp
                           np.zeros(shape),                # temp_prev
                           np.zeros(shape),                # temp_tend
                           np.zeros(shape)]                # st_energy
            if _device_path:
                run_args = _device_copies(run_args, device_ptrs)
                update_args = _device_copies(update_args, device_ptrs)
            run(ncol, nz, dt, lyr_surf, lyr_toa, *run_args)
            update(ncol, nz, dt, *update_args)
    finally:
        latvap, pref, rhoqr, gravit, errmsg, errflg, scheme_name = saved
        if device_ptrs:
            # run()/update() have waited for their copies, and the JAX arrays
            # imported from these buffers are gone with run_args/update_args
            run_args = update_args = None
            _cudart.cudaDeviceSynchronize()
            for p in device_ptrs:
                _cudart.cudaFree(p)
        if _set_perf_enabled is not None:
            _set_perf_enabled(True)
    if _log_perf_event is not None:
        _log_perf_event("a_i:EAMxx::kessler::init::py_warmup", time.perf_counter() - t0)

# calls the kessler_run function from https://github.com/NCAR/llm-fortran-modernization/tree/main/fortran2jax-kessler/_officialJAX
# these arrays should be automatically updated back in EAMxx if everything is setup correctly.

def run(ncol, nz, dt, lyr_surf, lyr_toa, rho, zm, pk, theta, qv, qc, qr, precl, relhum):

    global latvap, pref, rhoqr, cpair, rair
    global errmsg, errflg, scheme_name

    if _device_path:
        return _run_device(ncol, nz, dt, lyr_surf, lyr_toa, rho, zm, pk, theta, qv, qc, qr, precl, relhum)

    # Compute — kessler_run_bridge (contract 1, host-facing; or the
    # py_run_variant's equivalent, _kessler_run_host) handles the
    # (ncol, nz) <-> (nz, ncol) layout conversion in-jit and does its own
    # H2D/D2H; called directly on EAMxx's native (ncol, nz) arrays.
    #
    # Timed the same way as fortran_bridge/kessler_eamxx_bridge_main.F90's
    # system_clock wrap around its `kessler_run` call, and logged under the
    # same label, so py_run and F90_run rows line up in the perf CSV.
    _t0_run = time.perf_counter()
    theta_out, qv_out, qc_out, qr_out, precl_out, relhum_out, scheme_name, errmsg, errflg, latvap, pref, rhoqr = _kessler_run_host(ncol, nz, dt, lyr_surf, lyr_toa, cpair, rair, rho, zm, pk, theta, qv, qc, qr, precl, relhum, scheme_name, errmsg, errflg, latvap, pref, rhoqr)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_run", ncol, nz, dt, time.perf_counter() - _t0_run)
        except Exception:
            pass

    # theta/qv/qc/qr/precl/relhum are zero-copy views into EAMxx's field
    # buffers (created in create_py_field()); the C++ caller (py_module_call)
    # discards this function's return value, so results must be written
    # back in place for them to reach EAMxx at all.
    #
    # KNOWN PERF ISSUE (tracked, not yet fixed): this write-back block
    # ("writeback" in kessler_perf_log.csv) costs about as much as the
    # kessler_run JAX call itself on the GPU production run (~0.0087s/call
    # vs ~0.0099s/call, over 576 calls/rank) -- a real, measured chunk of
    # a:EAMxx::kessler::run's total that isn't the JAX computation at all.
    #
    # Ruled out: GPU-resident destination memory. Re-ran the identical GPU
    # case with py_backend=host (theta/qv/qc/qr/precl/relhum become plain
    # host arrays instead of pybind11 views over Kokkos device memory) --
    # writeback cost was unchanged (5.024s vs 5.026s total over the run),
    # and so was a:EAMxx::kessler::run's overall total (87.768s vs 87.001s).
    # So this isn't about host vs. device memory on either side.
    #
    # Leading hypothesis instead: fixed per-call overhead from doing SIX
    # separate JAX-array -> NumPy conversions/assignments, rather than a
    # data-volume/bandwidth cost. ~28MB across the six arrays should take
    # ~2ms at ordinary memory bandwidth; we're seeing ~8.7ms/call (~1.45ms
    # per array), which is in the range of per-call dispatch overhead
    # (JAX array materialization, buffer-protocol negotiation), not copy
    # time -- and that would explain why it doesn't move with py_backend.
    #
    # Proposed fix (not yet implemented -- holding off for now): theta, qv,
    # qc, qr, and relhum are all the same (ncol, nz) shape. Stack them into
    # one array before returning from kessler_run_core (e.g.
    # jnp.stack([...], axis=0)), materialize that ONE stacked array here
    # instead of five separate ones, then slice it apart for the five
    # assignments below. That cuts six separate JAX->NumPy crossings down
    # to two (the stacked block + precl), which should cut this cost
    # roughly proportionally if the fixed-per-call-overhead hypothesis is
    # right. Would need a matching change in kessler_jax/kessler_run.py's
    # return values (and kessler_run's unpacking) to un-stack before this
    # runs -- not just a change here.

    _t0_writeback = time.perf_counter()
    theta[...]  = theta_out
    qv[...]     = qv_out
    qc[...]     = qc_out
    qr[...]     = qr_out
    precl[...]  = precl_out
    relhum[...] = relhum_out
    if _log_perf_call is not None:
        try:
            _log_perf_call("writeback", ncol, nz, dt, time.perf_counter() - _t0_writeback)
        except Exception:
            pass

def _run_device(ncol, nz, dt, lyr_surf, lyr_toa, rho, zm, pk, theta, qv, qc, qr, precl, relhum):
    """run() for py_backend=device on a GPU: the same jitted bridge
    (kessler_run_bridge is to_device -> kessler_run_bridge_device ->
    device_get; or the py_run_variant's), on EAMxx's device fields, with no
    host transfers."""
    global latvap, pref, rhoqr, errflg

    _t0_run = time.perf_counter()
    theta_o, qv_o, qc_o, qr_o, precl_o, relhum_o, errflg_o, latvap, pref, rhoqr = _kessler_run_device(
        ncol=ncol, nz=nz, dt=dt, lyr_surf=lyr_surf, lyr_toa=lyr_toa, cpair=cpair, rair=rair,
        rho=_import(rho), z=_import(zm), pk=_import(pk), theta=_import(theta),
        qv=_import(qv), qc=_import(qc), qr=_import(qr), precl=_import(precl), relhum=_import(relhum),
        errflg=errflg, lv=latvap, pref=pref, rhoqr=rhoqr)
    jax.block_until_ready((theta_o, qv_o, qc_o, qr_o, precl_o, relhum_o))
    errflg = int(errflg_o)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_run", ncol, nz, dt, time.perf_counter() - _t0_run)
        except Exception:
            pass

    _t0_writeback = time.perf_counter()
    _copy_back([(theta, theta_o), (qv, qv_o), (qc, qc_o), (qr, qr_o),
                (precl, precl_o), (relhum, relhum_o)])
    if _log_perf_call is not None:
        try:
            _log_perf_call("writeback", ncol, nz, dt, time.perf_counter() - _t0_writeback)
        except Exception:
            pass

def update(ncol, nz, dt, zm, pk, theta, phis, temp, temp_prev, temp_tend, st_energy):
    global gravit, cpair
    global errmsg, errflg

    # Contract 2 (device-resident): upload each array ONCE (to_device is a
    # pure H2D copy -- no host-side permute; the bridge_device wrappers
    # reverse axes inside their jitted region), chain timestep_init -> run ->
    # timestep_final entirely on the GPU with no transfers between phases,
    # then fetch ONCE at the end. Per-phase errflg is left on the device
    # (contract 2's own convention) and fetched together with the arrays,
    # not checked per call. cpair is a spatially-uniform scalar, so it
    # is passed through as-is (no device transfer needed).
    # py_backend=device: import EAMxx's device fields instead (no H2D)
    _in = _import if _device_path else to_device
    temp_d  = _in(temp)
    theta_d = _in(theta)
    exner_d = _in(pk)
    zm_d    = _in(zm)
    phis_d  = _in(phis)
    zeros_d = jnp.zeros((ncol, nz), dtype=jnp.float64)

    # Per-phase timing below mirrors fortran_bridge/kessler_eamxx_bridge_update.F90's
    # system_clock wrap around kessler_update_timestep_init/_run/_timestep_final
    # (same three labels, so py_run and F90_run rows line up in the perf CSV).
    # Unlike the Fortran calls, these three are intentionally left un-synced
    # on the device -- contract 2's whole point is to chain them with no D2H
    # between phases and fetch once at the end (see the batched device_get
    # below) -- so perf_counter() here times host-side dispatch, not device
    # execution. A true per-phase device time would need jax.block_until_ready()
    # after each call, which would force a sync point between phases and
    # defeat that chaining.
    _t0 = time.perf_counter()
    temp_prev_d, ttend_t_d, e1 = kessler_update_timestep_init_bridge_device(
        temp=temp_d, temp_prev=zeros_d, ttend_t=zeros_d, errflg=0)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_update_timestep_init", ncol, nz, dt, time.perf_counter() - _t0)
        except Exception:
            pass

    _t0 = time.perf_counter()
    ttend_t_d, e2 = kessler_update_run_bridge_device(
        nz=nz, ncol=ncol, dt=dt, theta=theta_d, exner=exner_d,
        temp_prev=temp_prev_d, ttend_t=ttend_t_d, errflg=0)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_update_run", ncol, nz, dt, time.perf_counter() - _t0)
        except Exception:
            pass

    # Apply the temperature tendency (CAM-SIMA's apply_tendency_of_air_temperature):
    # kessler_update_run only accumulates ttend_t, and none of the kessler_update
    # phases write temp. Done before timestep_final, which computes st_energy
    # from temp. temp_prev_d/ttend_t_d come back in temp's (ncol, nz) layout.
    # HOMME derives its FT forcing from the change in T_mid across physics, so
    # this is what delivers Kessler's heating to the dycore. Same label as the
    # Fortran bridge's kessler_apply_T_tend step.
    _t0 = time.perf_counter()
    temp_new_d = temp_prev_d + ttend_t_d * dt
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_apply_T_tend", ncol, nz, dt, time.perf_counter() - _t0)
        except Exception:
            pass

    _t0 = time.perf_counter()
    st_energy_d, e3, gravit = kessler_update_timestep_final_bridge_device(
        nz=nz, cpair=cpair, temp=temp_new_d, zm=zm_d, phis=phis_d,
        st_energy=zeros_d, errflg=0, gravit=gravit)
    if _log_perf_call is not None:
        try:
            _log_perf_call("kessler_update_timestep_final", ncol, nz, dt, time.perf_counter() - _t0)
        except Exception:
            pass

    if _device_path:
        # Outputs stay on the GPU: copy them into EAMxx's device fields.
        jax.block_until_ready((temp_new_d, temp_prev_d, ttend_t_d, st_energy_d))
        errflg = max(int(e1), int(e2), int(e3))
        errmsg = "" if errflg == 0 else "kessler_update: bad time splitting"
        _copy_back([(temp, temp_new_d), (temp_prev, temp_prev_d),
                    (temp_tend, ttend_t_d), (st_energy, st_energy_d)])
        return

    # ONE batched D2H fetch for all three phases' outputs and errflg scalars.
    temp_out, temp_prev_out, temp_tend_out, st_energy_out, e1, e2, e3 = jax.device_get(
        (temp_new_d, temp_prev_d, ttend_t_d, st_energy_d, e1, e2, e3))

    errflg = max(int(e1), int(e2), int(e3))
    errmsg = "" if errflg == 0 else "kessler_update: bad time splitting"

    # temp/temp_prev/temp_tend/st_energy are zero-copy views into EAMxx's field
    # buffers (same as run()'s writeback) -- results must be written back
    # in place for them to reach EAMxx.
    temp[...] = temp_out
    temp_prev[...] = temp_prev_out
    temp_tend[...] = temp_tend_out
    st_energy[...] = st_energy_out

# Called once per rank from KesslerMicrophysics::finalize_impl() -- flushes
# this rank's in-memory perf totals (accumulated by log_call() in run()/
# kessler_run() above) out to the perf-log CSV.
def finalize():
    if _flush_perf_log is not None:
        try:
            _flush_perf_log()
        except Exception:
            pass

if _log_perf_event is not None:
    _log_perf_event("a_i:EAMxx::kessler::init::py_import", time.perf_counter() - _t0_import)
