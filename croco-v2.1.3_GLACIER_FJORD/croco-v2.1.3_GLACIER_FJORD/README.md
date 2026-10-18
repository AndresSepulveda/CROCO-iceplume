# GLACIER, FJORD and PLUME test cases for croco-v2.1.3

Files to add to, or replace in, a standard croco-v2.1.3 distribution. Copy the `OCEAN/` and
`TEST_CASES/` directories of this archive over the root of the CROCO source tree:

```bash
cp -r OCEAN TEST_CASES /path/to/croco-v2.1.3/
```

## New files

| File | Content |
|---|---|
| `OCEAN/iceplume.F` | ICEPLUME subglacial plume model (port of ROMS-ICEPLUME) |
| `OCEAN/iceplume.h` | ICEPLUME shared variables |
| `OCEAN/iceplume_const.h` | ICEPLUME physical constants |
| `OCEAN/iceplume_opkd.F90` | ODE solver (ODEPACK DLSODE) used by ICEPLUME |
| `OCEAN/fjord.h` | FJORD test case parameters |
| `TEST_CASES/croco.in.GLACIER`, `croco.in.GLACIER_nbq` | GLACIER input (hydrostatic / NBQ) |
| `TEST_CASES/croco.in.Fjord`, `croco.in.Fjord_nbq` | FJORD input, river or no head source (hydrostatic / NBQ) |
| `TEST_CASES/croco.in.Fjord_glacier`, `croco.in.Fjord_glacier_nbq` | FJORD input with `FJORD_GLACIER` (hydrostatic / NBQ) |
| `TEST_CASES/croco.in.Plume` | PLUME input (explicitly resolved plume, NBQ) |

## Replaced files

| File | Changes |
|---|---|
| `OCEAN/Makefile` | Adds `iceplume.F` and `iceplume_opkd.F90` to the sources |
| `OCEAN/cppdefs.h` | GLACIER, FJORD and PLUME test cases and their options (including `GLACIER_NBQ`, `FJORD_NBQ`) |
| `OCEAN/cppdefs_dev.h` | ICEPLUME option consistency; ICEPLUME allowed with NBQ/M3FAST (HMass coupling only) |
| `OCEAN/param.h` | Grid size, point sources and passive tracers for GLACIER, FJORD and PLUME |
| `OCEAN/read_inp.F` | `iceplume:` keyword; input file names of GLACIER, FJORD (hydrostatic / NBQ) and PLUME |
| `OCEAN/step.F` | Call to the plume model |
| `OCEAN/step2d.F`, `OCEAN/step3d_uv2.F`, `OCEAN/omega.F`, `OCEAN/pre_step3d.F`, `OCEAN/step3d_t.F` | Coupling of plume volume and tracer fluxes (hydrostatic) |
| `OCEAN/step3d_fast.F` | Coupling of plume fluxes in the 3D fast mode (NBQ and M3FAST): barotropic transport, surface velocity and per-level NBQ momentum at the glacier face |
| `OCEAN/ana_grid.F` | Analytical grids, bathymetry and masks of GLACIER, FJORD and PLUME |
| `OCEAN/ana_initial.F` | Initial T, S and passive tracers of GLACIER, FJORD and PLUME |
| `OCEAN/analytical.F` | FJORD: wind stress, open-boundary tide, exchange flow and T/S, analytical T/S profiles; PLUME: bottom discharge shape (`ana_psource`) |

## Use

In the configuration `cppdefs.h`: `#undef REGIONAL` and `#define GLACIER`, `#define FJORD` or
`#define PLUME`.
FJORD options (`FJORD_SILL`, `FJORD_TIDES`, `FJORD_EXCHANGE`, `FJORD_RIVER`, `FJORD_GLACIER`,
`FJORD_WIND_UP`, `FJORD_WIND_DOWN`, `FJORD_GREENLAND`, `FJORD_NBQ`) are set in the FJORD block;
parameters are in `OCEAN/fjord.h`. The input file is chosen automatically by `read_inp.F`.

## Non-hydrostatic option (NBQ)

`GLACIER_NBQ` (GLACIER block) and `FJORD_NBQ` (FJORD block) activate CROCO's non-hydrostatic
pseudo-compressible solver (`NBQ`, with `NBQ_PRECISE`). The `_nbq` input files add the keyword
`time_stepping_nbq` (reduced sound speed 150 m/s, larger than sqrt(g H) for both fjords) and use
NDTFAST = 40 (dtfast = 0.75 s) to satisfy the acoustic CFL; outputs carry the suffix `_nbq`.

What NBQ changes: the plume itself (rise, entrainment, melt, neutral-buoyancy depth) is still
computed by buoyant plume theory, a 1D sub-grid model whose vertical velocities (order 1 m/s) are
internal to the parameterization. NBQ makes the *resolved* response non-hydrostatic: w becomes a
prognostic variable instead of being diagnosed from continuity, which matters near the glacier face
where the plume fluxes converge and diverge strongly in one cell. Only the HMass coupling
(`PSOURCE`) is supported with NBQ; `PSOURCE_MASS` (VMass) stops the compilation.
Explicitly resolving the plume would require dropping ICEPLUME and grids of a few metres.

Tests (gfortran, serial):

| Configuration | Length | Result |
|---|---|---|
| GLACIER, hydrostatic vs NBQ | 10 d | Both stable; plume mass balance closed (< 1e-4 m3/s). NBQ: weaker and smoother resolved w at the plume cell (2-day mean max 0.019 vs 0.059 m/s), neutral buoyancy deeper (-16.6 vs -11.0 m), melt 9.26 vs 9.39 m3/s, similar surface outflow (0.43 vs 0.49 m/s). |
| FJORD + glacier + Greenland, sill, tide, exchange: hydrostatic vs NBQ | 6 h | Same tide at the mouth and same plume (neutral buoyancy -163 m, melt 38.5-38.9 m3/s). Max w near the glacier 0.12 m/s with NBQ vs 0.22 m/s hydrostatic; the hydrostatic w shows grid-scale alternating cells below the plume detrainment that NBQ removes. |
| FJORD + river with NBQ | 2 h | Runs. |

Cost: NBQ is several times slower than the hydrostatic run (3D fast mode and smaller dtfast).
Caveat: the NBQ runs are short (6 h for FJORD); small local overshoots in the plume cell
(S = 34.63 for an ambient maximum of 34.60 after 6 h) appear as in the hydrostatic case.

## PLUME: explicitly resolved plume (no ICEPLUME)

Benchmark for the ICEPLUME parameterization. Near field of the GLACIER channel as a vertical slice
(300 x 3 x 60 points): dx = dy = 5 m, 60 uniform levels (dz = 3.3 m), h = 200 m, glacier (masked
column) at the west end, closed basin 1.5 km long. Same stratification (T = 4 °C, halocline at 50 m)
and same discharge per unit width as GLACIER (100 m3/s over a 220-m sheet = 0.45 m2/s), released
through the two bottom levels of the ice face as point sources (fresh, 0 °C, dye = 1). NBQ dynamics,
GLS k-epsilon sub-grid mixing, dt = 0.25 s, NDTFAST = 6, reduced sound speed 50 m/s; 1 hour by default
(about 20 min on one core).

Setup notes: WENO5 advection with NBQ at this resolution, and the NBQ open boundary, were unstable
(spurious free-surface oscillations of several metres); the default advection schemes and a closed
basin are therefore used. The discharge raises the level by about 1.1 m per hour (0.5 % of the depth).

Result (mean of minutes 30-60) against buoyant plume theory with the ICEPLUME equations (line plume,
no melt): the wall plume, 10-15 m wide, reaches the surface in about 5 minutes. Below the halocline
resolved and theoretical plumes agree (at -150 m: w = 0.83 vs 0.89 m/s, volume flux 4.8 vs 4.8 m2/s);
above it the resolved plume accelerates to about 1.5 m/s and entrains 30-40 % less, so the surface
dilution is lower (discharge fraction 0.04 vs 0.027). Caveats: 2D slice (no 3D turbulence), plume only
2-3 cells wide, k-epsilon closure and a short closed-basin run; a 3D, finer benchmark is the next step.

Checked: copying these files over a clean croco-v2.1.3 tree gives a source tree identical to the
tested one.
