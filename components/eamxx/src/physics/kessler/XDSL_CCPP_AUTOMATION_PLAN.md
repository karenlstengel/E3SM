# Automating EAMxx Bridge Generation with xdsl-ccpp

This document is a forward-looking design proposal, not a record of work already done (see
`README.md` for that). It addresses: given both the CPU and GPU (OpenACC) Fortran versions of a
scheme like Kessler, how would `xdsl-ccpp` need to be used and extended to automatically generate
all of the code required to call it from EAMxx, including the EAMxx `AtmosphereProcess` C++
interface itself, not just the Fortran bridge/cap layer? Investigated 2026-07-29 by reading the
`xdsl-ccpp` generator source (`xdsl_ccpp/transforms/`) in addition to its docs. No code changes were
made; nothing here has been implemented.

---

## What already exists in xdsl-ccpp that's reusable

The generator has more GPU-directive infrastructure than its own docs (`multilanguage_limitations.md`)
let on. `xdsl_ccpp/transforms/gpu_data_pass.py` and `gpu_ccpp_cap_pass.py` read `memory_space`
metadata annotations and automatically insert `!$acc`/`!$omp target` data-movement directives
(`enter/exit data`, `update device/self`) at the correct lifecycle-phase boundaries, including
handling for cross-phase hoisting and per-scheme "diverged" variables. This is exactly the machinery
that produced the directives seen in `generated_bridge/Kessler_ccpp_cap.F90`.

However, this GPU-directive machinery is wired only into the plain `ccpp_cap.py` path -- the one
that expects a Fortran host module (`type = module`, e.g. the never-generated
`eamxx_kessler_host_mod`) -- and is **not** connected to `xdsl_ccpp/transforms/cpp_interop.py`'s
chost-cap path, which is what actually generates `Kessler_ccpp_chost_cap.F90`/`.h` and is what the
EAMxx C++ interface calls. That's the structural reason the chost cap has zero GPU-directive support
today.

`cpp_interop.py` itself (the chost-cap generator, pass name `generate-cpp-cap`) is a solid,
well-factored per-lifecycle emitter: it already does kind mapping (`kind_phys` -> `real(c_double)`/
`double`), DDT flattening, and automatic `ncol`/`nz` injection, all driven off the completed IR. A
new EAMxx-targeted printer should build on this same completed-IR state rather than re-deriving it
from the raw `.meta` files.

---

## Phase A -- Teach the metadata format about CPU/GPU scheme variants

This directly fixes bugs #2/#3 from `README.md` at the tool level instead of via a hand-patched meta
fork.

1. **Fix the upstream drift first.** Update
   `GPU_ports/atmospheric_physics/schemes/kessler/kessler_update.meta` so it actually matches its own
   `.F90` (add the missing `ncol`/`nz` entries on `kessler_update_timestep_init` and
   `kessler_update_timestep_final`, in the real subroutine's argument order, including the
   `errflg`-before-`errmsg` ordering on `timestep_final`). Do the same check for `kessler.meta`, even
   though `kessler_run`/`kessler_init` were confirmed identical between the CPU and GPU_ports trees.
2. **Add a variant tag to the metadata format.** Extend `[ccpp-table-properties]` with a new
   property, e.g. `variant = openacc`, parallel to the existing `array_layout` and `language`
   properties, so the CPU and GPU_ports `.meta` files for the same scheme can both be fed to the
   generator in one invocation without one overwriting the other.
3. **Extend `suite_cap.py`'s call emission.** The code that currently emits a single hardcoded
   `call kessler_update_timestep_init(...)` needs to detect when two variant tables for the same
   scheme+lifecycle have different argument lists, and emit a `#ifdef <directive-flag>` /
   `#else` branch automatically -- mechanizing the fix from the README's to-do list, generically,
   for any future scheme divergence, not just Kessler.

Effort: moderate. No new printer required. This alone would make `ccpp_xdsl` capable of generating a
correct, dual-variant chost cap for Kessler in one command, with no hand-editing of generated output.

---

## Phase B -- Make the chost cap's device-pointer contract explicit

Today, "the host always hands the cap an already-resident device pointer" is true only because EAMxx
happens to behave that way and the GPU_ports scheme happens to use `!$acc ... deviceptr(...)`
clauses internally -- nothing in the metadata says this is guaranteed. Proposed fix: add a host-meta
property (e.g. `gpu_pointer_mode = deviceptr`, sibling to the existing `array_layout`) that tells
`cpp_interop.py` to emit zero data-staging directives at the chost-cap boundary and simply pass
pointers through as-is.

This turns `multilanguage_limitations.md` section 2's current state -- "the generated code provides
no help with this, hope your scheme happens to use `deviceptr`" -- into a documented,
generator-checked contract. It also creates a place for the generator to flag a mismatch at
generation time (e.g. `memory_space = device` declared on an argument, but the scheme's own
directives don't use `deviceptr`), instead of only failing silently at runtime.

---

## Phase C -- A new "EAMxx AtmosphereProcess" printer

This is the actual "generate the whole C++ interface" ask. Nothing like it exists in xdsl-ccpp today
-- `multilanguage_plan.md` only ever generates a flat C header plus a thin C++ ergonomics wrapper
(`Kessler_chost.hpp`), never a framework-integrated host class. It would be a new backend module,
e.g. `xdsl_ccpp/backend/print_eamxx_process.py`, consuming the same completed IR that
`cpp_interop.py` already builds.

### Mechanical part (low risk -- the IR already has what's needed)

- `initialize_impl` / `run_impl` / `finalize_impl` bodies that call the generated
  `Kessler_chost_physics_*` entry points in the correct lifecycle order. This is a direct readout of
  the suite's lifecycle function list, which `cpp_interop.py` already computes internally.

### Requires genuinely new, EAMxx-specific metadata vocabulary

- `create_requests()`'s `add_field<Required/Updated/Computed>` / `add_tracer` calls need each host
  variable classified beyond what CCPP's `intent` already captures: is this a Field-Manager-registered
  field, a tracer, or purely local/derived data never seen by the host's field manager? What
  `FieldLayout` tags does it need (`COL`, `LEV`)? (`units` is already present in the metadata and
  needs no extension.)
- Buffer management (`requested_buffer_size_in_bytes` / `init_buffers`, `ATMBufferManager`) and the
  Kokkos transpose glue (`params_helpers` / `params_computed` structs, `h_*` vs `f_*` view selection)
  require the generator to understand EAMxx's buffer-manager and Kokkos-View APIs specifically. There
  is no existing analog anywhere in the tool for this. This is the highest-effort, most bespoke piece
  of the whole plan, and the first version would likely still need per-process hand-tuning (pack
  sizes, buffer counts) even once generated.

### Should stay hand-written -- permanently, not just for now

- **`add_invariant_check` / `add_postcondition_check` physical bounds** (e.g. `qv` clamped to
  `[1e-13, 0.2]`). This is an EAMxx QA convention with no CCPP equivalent. Encoding "what bounds are
  physically sane for this variable" into generic scheme metadata would be scope creep well beyond
  what a bridge-code generator should take on.
- **Energy-fixer boundary-flux zeroing** -- EAMxx-integration bookkeeping unrelated to Kessler
  itself.

The claim that originally sat here -- that all of the host-side physics derivation
(`PF::exner_function`, `calculate_theta_from_T`, `calculate_dz`, `calculate_z_int`/`calculate_z_mid`)
should stay permanently hand-written because it's "not part of Kessler's CCPP-described interface at
all" turned out to be wrong. See the next section.

---

## Revision (2026-07-29): the real suite XML changes the "hand-written" assessment

The assessment above was based only on the two-scheme suite (`kessler`, `kessler_update`) described
by `xdsl-cpp/examples/kessler/scheme/kessler_suite.xml`, the ad hoc suite definition used to generate
this bridge. The actual upstream suite definition,
`atmospheric_physics/suites/suite_kessler.xml`, lists 20 schemes across two groups:

```
physics_before_coupler:
  calc_exner, temp_to_potential_temp, calc_dry_air_ideal_gas_density,
  wet_to_dry_water_vapor, wet_to_dry_cloud_liquid_water, wet_to_dry_rain,
  kessler,
  potential_temp_to_temp, dry_to_wet_water_vapor, dry_to_wet_cloud_liquid_water, dry_to_wet_rain,
  kessler_update,
  qneg, geopotential_temp,
  check_energy_zero_fluxes, check_energy_scaling, check_energy_chng,
  sima_state_diagnostics, kessler_diagnostics
physics_after_coupler:
  thermo_water_update, check_energy_scaling, dycore_energy_consistency_adjust,
  apply_tendency_of_air_temperature, sima_tend_diagnostics
```

`eamxx_kessler_process_interface.cpp` already tracks this exact list via inline
`// <scheme>name</scheme>` comments -- whoever wrote the bridge was working through this real suite
XML scheme-by-scheme. Reading `schemes/utilities/state_converters.F90` (which implements
`calc_exner`, `temp_to_potential_temp`, `calc_dry_air_ideal_gas_density`, and all the wet/dry
conversion schemes) splits the old, single "physics derivation" bucket into three categories that
need different treatment, plus a fourth for `geopotential_temp` specifically.

### Category 1 -- Real CCPP schemes, formulas match, kept hand-fused purely for performance

`calc_exner_run` (`exner = (pmid/ref_pres)**(rair/cpair)`) and `temp_to_potential_temp_run`
(`theta = temp/exner`) are essentially the same math as `PF::exner_function`/
`PF::calculate_theta_from_T`, just expressed as standalone CCPP schemes instead of Kokkos device
functions. One real subtlety: `calc_exner_run` takes per-column, composition-dependent `rair`/`cpair`
as arguments, while EAMxx's `PF::exner_function(p_mid)` uses fixed dry-air constants internally -- so
the current hand-written preprocessing kernel and the real CCPP scheme are not even bit-identical
today; the CCPP version is arguably *more* thermodynamically self-consistent, since it reuses the
same composition-dependent `cpair`/`rair` fed into Kessler downstream. These genuinely could be
generated and called as separate bridge calls -- the reason not to is pure performance: doing so
would add two more Fortran round-trips (with column-major transposes each way) for cheap per-column
arithmetic that's currently fused into one Kokkos `parallel_for` alongside the rest of the
preprocessing. That is a legitimate engineering tradeoff, not a generator limitation, so the original
"not part of Kessler's CCPP interface at all" framing was wrong for this subset.

### Category 2 -- Nominally generatable; needs manual/metadata verification of matching conventions, not an assumption

`calc_dry_air_ideal_gas_density_run` computes `rho = pmiddry/(rair*temp)` using **dry** mid-level
pressure, not EAMxx's mass-weighted `rho = pseudo_density/(g*dz)` derivation -- a genuinely different
vertical-coordinate approach requiring a `pmiddry` quantity EAMxx does not currently compute.
`create_requests()` even has a commented-out line, `add_field<Required>("pseudo_density_dry", ...)`,
suggesting this was started and abandoned. Whether the CCPP scheme's ideal-gas density and EAMxx's
mass-weighted density agree closely enough to be interchangeable is a real open physics question, not
a code-generation question, and should not be assumed either way without checking.

The wet/dry mixing-ratio converters (`wet_to_dry_water_vapor`, `dry_to_wet_water_vapor`, etc.) looked
like a similar unaddressed gap at first glance, since the bridge has no comment justifying skipping
them (every other skipped scheme in Category 3 below has one). **This has since been confirmed not to
be a bug for Kessler specifically**: both `kessler_run` and `kessler_update`'s Fortran already declare
`qv`/`qc`/`qr` with standard name `water_vapor_mixing_ratio_wrt_dry_air` (etc.), and this matches what
EAMxx already provides for those fields, so no wet/dry conversion is needed in this particular case --
the standard names agree on both sides of the host/scheme boundary.

That agreement will not hold for every future scheme, though. CCPP's standard-name convention is
exactly the mechanism that should catch a mismatch (a scheme expecting
`water_vapor_mixing_ratio_wrt_dry_air` will simply fail to match a host variable declared
`_wrt_moist_air`, forcing either an explicit conversion scheme in the suite or a host-side fix) -- but
that protection only works if the host `.meta` file's standard name is *honestly* the basis the host
variable is actually on. A host author who mislabels a wet-basis field as dry-basis produces a
silent, high-confidence-looking match that is simply wrong, and no automated check catches a
mislabeling problem, since matching is defined as standard-name equality, not a semantic audit. So
the right practice going forward is: confirm the moisture basis (and any other convention a standard
name implies) by hand for each new field the first time it's wired up, and/or add this to a
host-to-scheme meta consistency check (see the suite-coverage idea below) rather than assuming
agreement by default.

### Category 3 -- Real schemes, structurally superseded by EAMxx's own centralized infrastructure

`qneg`, `check_energy_zero_fluxes`/`check_energy_scaling`/`check_energy_chng`,
`sima_state_diagnostics`/`kessler_diagnostics`/`sima_tend_diagnostics`, `thermo_water_update`,
`dycore_energy_consistency_adjust`, `apply_tendency_of_air_temperature` are all real, generatable CCPP
schemes, but each has a comment in `eamxx_kessler_process_interface.cpp` mapping it to an
EAMxx-generic mechanism instead: postcondition/invariant field checks substitute for `qneg`,
`output_fields.yml`-driven diagnostics substitute for the `sima_*`/`kessler_diagnostics` schemes, and
a single centralized energy-fixer `AtmosphereProcess` (wrapping the whole process list, not per-scheme)
substitutes for the `check_energy_*` family. Calling the generated versions of these *in addition to*
EAMxx's centralized equivalents would risk double-applying corrections -- e.g. two independent energy
adjustments on the same budget. This is a durable, architectural reason to exclude them, not a
metadata gap, and a smarter generator should never try to paper over it. Two of these
(`thermo_water_update`, and the `dycore_energy_consistency_adjust`/`apply_tendency_of_air_temperature`
pair) are marked with "I think" / "TODO ?" in the existing comments -- genuinely unresolved
uncertainty in the current bridge, independent of code generation.

### `geopotential_temp` -- a fourth, distinct case: linked but dead

`geopotential_temp` is in the suite, and its `.F90` is already linked into
`kessler/CMakeLists.txt`'s source list, but it is never actually called anywhere -- `z_mid` is instead
derived by a hand-written Kokkos team-parallel-scan (`calculate_z_int`/`calculate_z_mid`). Unlike
Category 3's entries, there is no comment establishing this as a deliberate, verified substitution --
it reads as linked dead code rather than a documented architectural decision.

### A cheaper, higher-value piece of tooling than the full `AtmosphereProcess` printer

Reading the real suite XML suggests a "suite-coverage checker": a small script that diffs the real
suite XML's `<scheme>` list against the `<scheme>` tracking comments already present in
`eamxx_kessler_process_interface.cpp`, and flags any suite entry with no corresponding comment at
all. That check alone would have flagged the wet/dry conversion question mechanically (it turned out
to be a non-issue for Kessler, but the bridge gave no evidence of having checked), and would catch the
same class of silent gap for every future scheme brought in this way. It requires no new xdsl-ccpp
printer and no new metadata vocabulary -- just a script comparing two lists of scheme names -- and
would be worth building well before Phase C.

---

## Recommended order

Do Phase A and Phase B first. They are a direct, generalizable fix for a bug already found and
documented in `README.md`, require no new printer, and immediately make `ccpp_xdsl` produce a
correct dual-variant chost cap for Kessler with no manual patching of generated output. The
suite-coverage checker above is a similarly cheap, high-value addition and should happen around the
same time -- it needs no generator changes at all. Phase C is the larger "generate the whole
interface" ask, but it is an order of magnitude more implementation effort, and even fully built would
still leave some hand-written code inside `run_impl`: the QA-check bounds and energy-fixer bookkeeping
permanently (see above), and the Category 1 preprocessing kernels (`calc_exner`/
`temp_to_potential_temp`) by deliberate performance choice rather than necessity. "Automatically
generate all of the code" is realistically achievable for the bridge/cap layer and the
lifecycle-orchestration skeleton, but not for the entire file.
