# ECN Spray G3 – OpenFOAM v14 baseline (rebuilt from `code/`)

Target: Lien et al., Fuel 357 (2024) 129708 – G3 results, Figs. 11, 14, 15.
Your original case in `../code/` is untouched.

## Run

```sh
./Allrun smoke      # serial, ~20 steps: catches dictionary/keyword errors in a minute
./Allrun 4          # full run on 4 cores, then reconstruct + VTK export
python3 ../postprocess/g3_validate.py . --label "M1 Cd0.60"
```

Output goes to `validation/`: `fig11_penetration.png`, `fig14_liquid_mass_ratio.png`,
`sim_metrics.csv`, `summary.txt` (RMSD and liquid mass ratio at end of injection).

## What changed vs `code/`

| Item | old | new | why |
|---|---|---|---|
| Geometry | injector on x=0 wall, spray +x, gravity -z | injector at top centre, spray **-z** (ECN / Lien convention) | axial penetration = -z, matches ECN data |
| Mesh | 40³ uniform, 2.5 mm | multi-grading, **1.0 mm** core (762 k cells); M0/M2 levels in `blockMeshDict` | 2.5 mm under-resolves gas entrainment -> penetration too short |
| Droplet size at nozzle | Rosin-Rammler d=50 µm, 10–165 µm | **Lien Eq. 36** (number CDF, 1.6–23 µm, D32 ≈ 10 µm) as truncated Rosin-Rammler, Q 0 | the paper's Case 1; old sizes ~10x too large |
| Plume direction | 37° | 35° | Lien Table 3 |
| Cone half-angle | 3.5° | 15° | Lien: 30° full cone (calibration knob) |
| Cd | 0.52 (~250 m/s) | 0.60 (~215 m/s) | start of the Cd sweep 0.52–0.73 |
| Drag | sphereDrag | sphereDrag (distortedSphereDrag available) | distortion needs the oscillation eq. from a breakup model |
| Turbulent dispersion | none | stochasticDispersionRAS | RANS needs it, otherwise plumes too narrow |
| Breakup | ReitzKHRT | **PilchErdman** | v14 ReitzKHRT, ReitzDiwakar, TAB, ETAB crash on step 1 (FPE/NaN) with 1.6–23 µm droplets; PilchErdman passes and is the bag/multimode correlation Lien's breakup approach builds on (Nagaoka & Kawamura). `--breakup none` = ablation |
| k-ε C1 | 1.44 | 1.55 | round-jet correction for sprays |
| parcels/s per hole | 2.5e6 | 1e7 | smoother statistics for PLV penetration |
| rho0 / Cp0 | 1000 / 4187 (water) | 633 / 2300 (iso-octane @ 363 K) | cosmetic: sprayCloud overwrites these from the IC8H18 liquid model |
| Output | ascii, nothing extra | ascii.gz, `cloudInfo` function object, writeCellCentres + VTK | needed for the ECN metrics |

## Calibration knobs (edit with the generator, not by hand)

```sh
python3 tools/make_cloudProperties.py --Cd 0.73          # injection velocity
python3 tools/make_cloudProperties.py --thetaOuter 10    # plume width
python3 tools/make_cloudProperties.py --B1 10 --CRT 0.6  # breakup
python3 tools/make_cloudProperties.py --sizes blob       # Lien Case 2 (165 µm blobs)
```
Note: `Allrun` regenerates `cloudProperties` from the script **defaults** – for a sweep,
change `PARAMS` in the script (or delete that line from Allrun).

## Not verified against a live v14 install

These keywords were not in your old case. If the smoke test complains, the log names
the bad keyword – send it to me:
`distortedSphereDrag`, `stochasticDispersionRAS`, `kEpsilonCoeffs`, the `cloudInfo`
function object (lib `liblagrangianFunctionObjects.so`), `Q 0` in RosinRammler,
the `#calc` lines in blockMeshDict, `foamToVTK -useTimeName -noInternal`.

## Open issues / to investigate later

### ReitzKHRT floating point exception with small droplets (possible upstream bug → contribution)

**Symptom (OpenFOAM 14, build `14-7b05503f98a8`):** `foamRun` aborts with `Floating point exception`
in `Foam::ReitzKHRT<...>::update(...)` on the **first time step**, right after the first 8 parcels
are injected.

**Reproduce:**
```sh
CLOUD_ARGS="--breakup ReitzKHRT" sh Allrun smoke
```

**What we know (smoke-test bisection, `tools/bisect_smoke.sh`, `tools/bisect_breakup.sh`):**
- Fails independent of KHRT constants (B0/B1/CRT = 0.55/7/0.42 *and* the defaults 0.61/55/0.1),
  drag model (sphereDrag / distortedSphereDrag), oscillation eq. on/off, dispersion on/off.
- Also fails with `FOAM_SIGFPE=false` → genuine NaN/inf, not just a trapped harmless FPE.
- Same model ran fine in the old case with Rosin-Rammler 10–165 µm droplets → triggered by the
  small injected droplets (Lien Eq. 36, 1.6–23 µm, We ≈ 15).
- ReitzDiwakar, TAB, ETAB also crash on step 1; PilchErdman and `none` pass.
- First injected parcels carry very little mass (flow-rate profile starts at 0).

**Suspects to check in `src/lagrangian/spray/submodels/BreakupModel/ReitzKHRT/ReitzKHRT.C`:**
- RT branch for tiny drops: drag acceleration `Urel/tMom` is huge (~1e7–1e8 m/s²), so
  `lambdaRT < d` already at step 1 → RT path active immediately.
- KH mass-shedding / child-parcel branch: `ms/averageParcelMass`, `d/d0`, and the cubic solve for
  the parent diameter (acos/cbrt/sqrt of out-of-range arguments).
- `trajectory = U/mag(U)` and other divisions without `vSmall` guards.

**Plan:** debug build (`WM_COMPILE_OPTION=Debug`) or `Info` prints in `update()` on the failing
parcel → find the offending expression → minimal test case (single injector, fixed 2 µm droplets,
aachenBomb-style) → report on bugs.openfoam.org with a patch.
