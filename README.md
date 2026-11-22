# BTP – ECN Spray G3 in OpenFOAM v14 (baseline for PINN spray surrogate)

| Folder | What |
|---|---|
| `case_G3/` | current case (see `case_G3/README.md` for setup, changes and run instructions) |
| `postprocess/` | `g3_validate.py` (ECN/Lien-style penetration, liquid mass ratio, RMSD) + ECN G3 data as CSV |
| `code/` | original first attempt (kept for reference) |

Only code is tracked; results, mesh, logs, papers and raw ECN `.mat` files are ignored (`.gitignore`).

```sh
cd case_G3
./Allrun smoke          # quick dictionary check
./Allrun 4              # full run
python3 ../postprocess/g3_validate.py . --label "M1 Cd0.60"
```
