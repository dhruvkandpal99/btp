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
| Drag | sphereDrag | distortedSphereDrag | stock analogue of Lien's distortion-corrected drag |
| Turbulent dispersion | none | stochasticDispersionRAS | RANS needs it, otherwise plumes too narrow |
| KH-RT B0 / B1 / CRT | 0.61 / 55 / 0.1 | 0.55 / 7 / 0.42 | Lien Eqs. 33–35 evaluated at We ≈ 15 |
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
