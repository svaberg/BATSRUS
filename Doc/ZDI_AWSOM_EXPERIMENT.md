# ZDI Ingestion and AWSoM Boundary Experiments (codex/zdi)

This document describes the ZDI-specific additions on branch `codex/zdi`, how to
validate the ingested field using shell `#SAVEPLOT` output, and how to run the
current AWSoM boundary-condition experiments (`baseline`, `clamp`, `nudge`).

## Scope

This branch adds:

- direct reading of a Donati-style ZDI coefficient file (3 coefficient sets:
  radial/poloidal/toroidal) in the main BATSRUS repo (no `util/` changes)
- AWSoM `#USERINPUT` controls for ZDI ingestion and boundary forcing
- new user plot variables for `#SAVEPLOT` (`zdibr`, `zdibphi`, `zdibtheta`, and
  poloidal/toroidal splits)
- experimental AWSoM inner-boundary forcing modes (`clamp`, `nudge`) for
  tangential components, plus optional radial-component forcing
- Python tools to read `idl_ascii` shell plots and make magnetic-field movies

## ZDI Coefficient File Format (supported)

The reader currently supports the Donati-style text files used in TOUPIES/ZDI
workflows, e.g. `coeff-Mel25-5.dat`.

Expected structure:

1. Header text line
2. Integer header line (3 integers)
3. Block 1 coefficients (triangular `(l,m)` list)
4. Blank line
5. Block 2 coefficients (triangular `(l,m)` list)
6. Blank line
7. Block 3 coefficients (triangular `(l,m)` list)

Each coefficient row is read as:

- `l  m  real  imag`

The three blocks are interpreted as:

- block 1 = radial (`alpha`)
- block 2 = poloidal (`beta`)
- block 3 = toroidal (`gamma`)

The reader validates triangular ordering (`l = 1..lmax`, `m = 0..l`) and infers
`lmax` from the row count.

## ZDI Conventions Implemented (matched to ZDIpy)

The evaluator in `/Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/src/ModZdiMagnetogram.f90`
is matched to `ZDIpy/core/magneticGeom.py` with the following conventions:

- normalized spherical-harmonic basis factor
  - `sqrt((2l+1)/(4*pi) * (l-m)!/(l+m)!)`
- tangential basis uses `1/(l+1)` factor
- Donati `nPotential == -3` coefficient conjugation handling is applied
- evaluator returns components as:
  - `Br`
  - `Bphi` (azimuth/longitude component)
  - `Btheta` (co-latitude component; positive toward increasing theta)
- therefore, latitude component is:
  - `Blat = -Btheta`

This matches the ZDIpy convention where plotting often uses `Blat = -Bclat`.

## AWSoM User Input: ZDI Controls

The new controls live in the AWSoM user-module path (`srcUser/ModUserAwsom.f90`).

### `#ZDIMAGNETOGRAM`

```txt
#ZDIMAGNETOGRAM
T                       UseZdiMagnetogram
Param/CORONA/coeff.dat  NameZdiCoeffFile
1.0                     ZdiFieldScaleIo
0.0                     ZdiLonShiftDeg
```

Parameters:

- `UseZdiMagnetogram`: enable direct ZDI coefficient ingestion
- `NameZdiCoeffFile`: Donati-style ZDI coeff file
- `ZdiFieldScaleIo`: multiplicative scale applied to evaluated ZDI field
- `ZdiLonShiftDeg`: longitude shift applied before ZDI evaluation

### `#ZDIBOUNDARY`

```txt
#ZDIBOUNDARY
T         UseZdiBoundary
F         UseZdiBoundaryRadial
clamp     TypeZdiBoundary
cosine    TypeZdiRamp
1.0       ZdiBcStrength
1.0       ZdiBcScale
1         ZdiRampIterStart
30        ZdiRampIterStop
-1.0      ZdiRampStart
-1.0      ZdiRampStop
```

Parameters:

- `UseZdiBoundary`: enable ZDI boundary forcing in AWSoM inner BC
- `UseZdiBoundaryRadial`:
  - `F` (default): only tangential (`Bphi`, `Btheta`) is forced/mixed
  - `T`: radial component is also mixed/clamped toward ZDI target (full 3-vector)
- `TypeZdiBoundary`:
  - `off`: disabled
  - `clamp`: direct mixing toward target each step
  - `nudge`: damped/relaxation-style mixing (`ZdiBcStrength` controls rate)
- `TypeZdiRamp`: `none`, `linear`, or `cosine`
- `ZdiBcStrength`: forcing strength (used especially by `nudge`)
- `ZdiBcScale`: scaling factor applied to the boundary target (separate from
  `ZdiFieldScaleIo` used for the evaluated ZDI field)
- `ZdiRampIterStart`, `ZdiRampIterStop`: iteration-based ramp window
- `ZdiRampStart`, `ZdiRampStop`: time-based ramp window (negative disables)

Note: `UseZdiBoundaryRadial` is a new line in the `#ZDIBOUNDARY` block. Older
local `PARAM.in` files that use the earlier `codex/zdi` syntax must be updated.

## Shell Diagnostics and `#SAVEPLOT` Variables

AWSoM `user_set_plot_var` now exposes the following ZDI diagnostic variables for
`#SAVEPLOT`:

- `zdibr`, `zdibphi`, `zdibtheta`
- `zdibphip`, `zdibthetap` (poloidal-only tangential contribution)
- `zdibphit`, `zdibthetat` (toroidal-only tangential contribution)
- aliases `blon` and `blat` are available (with `blat = -btheta` convention)

Example shell plot (IDL ASCII) at fixed radius:

```txt
#SAVEPLOT
1                       nPlotFile
shl VAR idl_ascii       StringPlot
1                       DnSavePlot
-1.0                    DtSavePlot
HGR                     TypeCoord
1.0                     rMin
1.0                     rMax
0.0                     LonMin
360.0                   LonMax
10.0                    dLon
-90.0                   LatMin
90.0                    LatMax
10.0                    dLat
br bphi btheta zdibr zdibphi zdibtheta  NameVars
{default}               NamePars
```

## Plot Format Backend Status (shell `#SAVEPLOT`)

Tested on this branch with shell plots (`shl VAR ...`):

- `idl_ascii`: works (used for diagnostics and movie tooling)
- `tec`: works (writes Tecplot `.dat` shell files)
- `tcp`: **does not work** for shell plots in this path (`unknown TypeFile=tcp`)
- `hdf`: parser/writer path exists and run reports saving, but no shell HDF file
  was observed on disk in local tests (needs more investigation)
- `netcdf` / `nc`: no shell `#SAVEPLOT` backend found in current parser/writer path

If you want a Tecplot shell output, use e.g.:

```txt
shl VAR tec             StringPlot
```

## Python Tooling (IDL ASCII shell reader + movie generator)

Added scripts:

- `/Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/Scripts/ZDI/bats_idl_ascii.py`
- `/Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/Scripts/ZDI/make_shell_component_movie.py`

### Reader summary example

```bash
python3 /Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/Scripts/ZDI/bats_idl_ascii.py \
  /Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/run_test_awsom_zdi_ingest/SC/IO2/shl_var_6_n00000001.out
```

### 3-panel movie (`br`, `bphi`, `btheta`) with symlog color bars

```bash
python3 /Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/Scripts/ZDI/make_shell_component_movie.py \
  /Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/run_test_awsom_zdi_ingest/SC/IO2 \
  --output /Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/zdi_movie_assets/baseline_br_bphi_btheta.mp4 \
  --frame-dir /Users/dagfev/Documents/SWMFsoftware/BATSRUS-zdi/zdi_movie_assets/baseline_br_bphi_btheta_frames \
  --layout 1x3 \
  --fps 10 \
  --title "AWSoM ZDI baseline (free tangential field)"
```

The movie script:

- reads shell `idl_ascii` files (`shl_var_*.out`)
- extracts `br`, `bphi`, `btheta` (or custom variables)
- uses symlog color normalization (`SymLogNorm`) with per-component symmetric limits
- adds annotations (`step`, `time`, source file, and optional mean mismatch to ZDI
  target when `zdib*` columns are present)
- resamples frames to the first shell grid if a grid mismatch is detected
- uses `ffmpeg` to write `.mp4`

## AMR Notes

The shell plot output (`shl VAR ...`) is sampled onto the explicit lon/lat grid
requested in `#SAVEPLOT`, so in normal usage it is already a fixed grid even when
AMR is active in the volume solution.

However, the Python reader/movie tools still support resampling because:

- different runs may use different `dLon/dLat`
- you may compare outputs produced with different shell plot settings
- future diagnostics may emit varying shell grids

Current implementation uses bilinear interpolation in lon-lat space (with periodic
longitude handling).

## Recommended Science Comparison Workflow (current branch)

1. **Baseline** (free tangential field): `UseZdiBoundary=F` or `TypeZdiBoundary=off`
2. **Clamp** (tangential target): `UseZdiBoundary=T`, `TypeZdiBoundary=clamp`,
   `UseZdiBoundaryRadial=F`
3. **Nudge** (tangential target): `UseZdiBoundary=T`, `TypeZdiBoundary=nudge`,
   `UseZdiBoundaryRadial=F`
4. **Full clamp / damped full clamp (end-state target)**:
   - set `UseZdiBoundaryRadial=T`
   - compare coronal structure against the historical/free-perpendicular setup

Diagnostics to compare:

- shell maps of `br`, `bphi`, `btheta`
- mismatch to ZDI targets (`zdibr`, `zdibphi`, `zdibtheta`)
- coronal structure changes (open/closed topology proxies, current sheets,
  wind speed and density patterns, etc., depending on study focus)

## Known Limitations / Next Steps

- The `hdf` shell `#SAVEPLOT` path needs investigation (reported saved, no file seen).
- `tcp` shell output appears unsupported in the tested path.
- Cross-run movie color scaling is currently per-run; if strict side-by-side
  visual comparison is needed, add fixed color scales shared across runs.
- Full 3-component ZDI boundary forcing (`UseZdiBoundaryRadial=T`) is implemented
  but should be tested more extensively for stability and long-run behavior.
