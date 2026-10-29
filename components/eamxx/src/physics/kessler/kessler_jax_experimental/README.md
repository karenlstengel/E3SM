# kessler_jax_experimental

Prototype restructurings of the sub-cycle loop in `kessler_jax/kessler_run.py`
(the LAFT translation). `kessler_jax/` itself is unchanged; everything here
reproduces `kessler_run_core` bit for bit and changes only how the rain
sub-cycle loop is run. The model uses it only when `kessler::py_run_variant`
is set to something other than `current` (see `../README.md`); the standalone
driver exercises it directly (`../driver_kessler_py.cpp --variants`) and
through `kessler.py` (`--run-variant`).

## The problem

`kessler_run_core` runs each column's sub-cycle as a `lax.while_loop` and
`jax.vmap`s it over the columns. Under `vmap`, JAX turns that into one loop
that runs until every column is done, masking the finished ones, so every
column does as many passes as the slowest one. Dry columns need 1 pass and
rainy ones 17 to 144; at 8% rainy columns that is 17.7x the necessary work
(`../driver_kessler_py.cpp --passes`). On the GPU each loop iteration also
costs a host round trip for the loop condition and several kernel launches.

## Files

- `kessler_run_variants.py`: `kessler_run_core_variant(variant, ...)` and
  `kessler_run_bridge_device_variant(variant, ...)`, with the same arguments
  and outputs as `kessler_run_core` and `bridge/kessler_run_bridge.py`'s
  `kessler_run_bridge_device`, plus a static `variant`. The setup and finish
  code is copied from `kessler_run_core`, and the per-column loop body from
  its `_column`, unchanged.
- `kessler_run_pallas.py`: the sub-cycle as one Pallas kernel (GPU).

## Variants

| variant        | what it does                                                         | use        |
|----------------|----------------------------------------------------------------------|------------|
| `auto`         | `pallas` on a GPU, `map` otherwise                                   | default    |
| `vmap`         | the current structure, rebuilt from the copied pieces (a check)      |            |
| `split`        | one pass on all columns, then the vmapped loop on only the columns still running, gathered into a fixed-size buffer (`lax.switch` over a few sizes, one compile) | |
| `map`          | `lax.map` over columns: each column runs its own loop, as in the Fortran | CPU    |
| `split_map`    | `split` with `lax.map` on the gathered columns                       |            |
| `pallas`       | one Pallas program per column, levels spread over its threads, state in registers (the C++ D1 structure) | GPU |
| `pallas_lanes` | one Pallas thread per column, loop over levels (too few threads at ne16) |        |

Pallas notes (Triton backend on A100, JAX 0.6.2): values cannot be sliced,
so the level k+1 neighbour goes through memory (store, `debug_barrier`,
gather-load, `debug_barrier`); blocks are powers of two, so the 58 levels are
padded to 64; float64 `exp`/`pow`/`sqrt` match XLA's bit for bit; float64
`dot` crashes Triton and is not used. Tuning knobs (read at trace time):
`KESSLER_PALLAS_NUM_WARPS` (default 2), `KESSLER_PALLAS_SHIFT` (`memory`,
default, or `select`: an exact masked sum, slower), `KESSLER_SPLIT_DIVISORS`.

## Results (driver, synthetic columns, 2026-10-29, deg0053)

`kessler_run_bridge_device` only (no EAMxx data movement), s/call, min of
10 (5 on CPU), 8% rainy columns. All variants bit for bit with the current
code at 3%, 8% and 12% rain on GPU ne16/ne60, and at 8% on CPU ne16/ne60.
GPU run-to-run spread is about 20%.

| grid     | current  | split    | map      | pallas   | speedup (best) | C++ D1 |
|----------|---------:|---------:|---------:|---------:|---------------:|-------:|
| GPU ne16 | 10.9e-3  | 5.8e-3   | -        | 1.5e-3   | 7.3x (pallas)  | 1.15e-3 |
| GPU ne60 | 97.9e-3  | 27.6e-3  | -        | 7.4e-3   | 13.2x (pallas) | 7.56e-3 |
| CPU ne16 | 33.9e-3  | 9.4e-3   | 3.4e-3   | -        | 9.8x (map)     | 3.70e-3 |
| CPU ne60 | 0.511    | 0.138    | 48.6e-3  | -        | 10.5x (map)    | 52.3e-3 |

(C++ D1: the C++ branch's team-per-column `kessler_run`, before I2.)

Through `kessler.py` (`../check_kessler_py_bfb.sh`, 2026-10-29): every
variant on every `py_backend` path matches `current` bit for bit over 3
chained model-like steps, repeats itself bit for bit (no race in the Pallas
barrier shift showed up), and passes the finite/non-negative checks; GPU and
CPU, ne16 and ne60, rain fractions 0 to 0.12.

## Not done yet

- Model runs. `kessler.py` selects a variant with `kessler::py_run_variant`
  (default `current`, i.e. `kessler_jax`), see `../README.md`; checked bit
  for bit through `kessler.py` in the driver, not yet in a model run.
- Rain checked only on the driver's synthetic columns, not real model state.
- The kernel's fixed cost (padding, transposes, ~0.45 ms at ne16) is not tuned.
