#!/usr/bin/env python3
"""
ECN Spray G3 validation for the OpenFOAM v14 case  (targets: Lien et al. 2024, Figs. 11, 14, 15)

Needs only numpy + matplotlib.

    python3 g3_validate.py <case> [--label "M1 Cd0.60"] [--ref lien_case1_liquid.csv ...]

Expects in <case> (Allrun produces all of these):
    VTK/lagrangian/cloud/*_<time>.vtk   parcels   (foamToVTK -ascii -useTimeName -noInternal)
    <time>/IC8H18[.gz]                   fuel vapour mass fraction
    0/Cz[.gz]  (or 0/C)                  cell-centre z   (foamPostProcess -func writeCellCentres -time 0)
    log.foamRun                          optional, for liquid mass every step

Definitions (ECN "follow the fuel", as used by Lien et al.):
  * Liquid penetration: parcel liquid volume is smeared over one CFD cell (box kernel, width = cell size),
    projected along the line of sight -> PLV map [mm^3 liquid / mm^2]; penetration = largest AXIAL
    distance from the injector tip where PLV >= 2e-4 mm3/mm2 (the ECN/Lien "low" threshold).
  * Vapour penetration: largest axial distance where Y_IC8H18 >= 1e-3 (ECN CFD standard), also
    reported as max(vapour, liquid) because schlieren sees the liquid region too.
  * Liquid mass ratio: sum(nParticle * rho * pi/6 d^3) / 10 mg.
  * Also the old metrics (max parcel tip, 95 % / 99 % liquid mass) for continuity with BTP1.
"""
import argparse, glob, gzip, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M_INJ = 10.0e-6          # kg, total injected mass
EOI = 0.78               # ms


# ----------------------------------------------------------------------------- readers
def _open(path):
    if os.path.exists(path):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    if os.path.exists(path + ".gz"):
        return gzip.open(path + ".gz", "rt")
    return None


def read_foam_internal(path):
    """internalField of an OpenFOAM ascii vol field -> np.array (N,) or (N,3) or float if uniform."""
    fh = _open(path)
    if fh is None:
        return None
    txt = fh.read()
    fh.close()
    m = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
    if m:
        v = m.group(1).strip("() ").split()
        return float(v[0]) if len(v) == 1 else np.array(v, float)
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", txt)
    if not m:
        raise ValueError(f"cannot parse internalField in {path}")
    kind, n = m.group(1), int(m.group(2))
    body = txt[m.end():]
    if kind == "scalar":
        arr = np.array(body[:body.index(")")].split(), float)
    else:  # vector
        end = body.index(")\n;") if ")\n;" in body else body.rindex(")")
        arr = np.array(body[:end].replace("(", " ").replace(")", " ").split(), float).reshape(-1, 3)
    if len(arr) != n:
        raise ValueError(f"{path}: expected {n} values, got {len(arr)}")
    return arr


def read_vtk_lagrangian(path):
    """Legacy ascii VTK polydata written by foamToVTK -> dict(points=(N,3), <field>=array)."""
    with open(path) as fh:
        tok = fh.read().split()
    out, i, n = {}, 0, 0
    while i < len(tok):
        t = tok[i]
        if t == "POINTS":
            n = int(tok[i + 1]); i += 3
            out["points"] = np.array(tok[i:i + 3 * n], float).reshape(n, 3); i += 3 * n
        elif t == "FIELD":
            nf = int(tok[i + 2]); i += 3
            for _ in range(nf):
                name, nc, nt = tok[i], int(tok[i + 1]), int(tok[i + 2]); i += 4
                a = np.array(tok[i:i + nc * nt], float); i += nc * nt
                out[name] = a.reshape(nt, nc) if nc > 1 else a
        elif t == "SCALARS":
            name = tok[i + 1]; nc = int(tok[i + 3]) if tok[i + 3].isdigit() else 1
            i += 4 if tok[i + 3].isdigit() else 3
            if tok[i] == "LOOKUP_TABLE": i += 2
            out[name] = np.array(tok[i:i + n * nc], float); i += n * nc
        elif t == "VECTORS":
            name = tok[i + 1]; i += 3
            out[name] = np.array(tok[i:i + 3 * n], float).reshape(n, 3); i += 3 * n
        else:
            i += 1
    if "points" not in out:
        out["points"] = np.zeros((0, 3))
    return out


def time_of(path):
    m = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", os.path.basename(path).replace(".vtk", ""))
    return float(m[-1]) if m else None


def parse_log_mass(logfile):
    if not os.path.exists(logfile):
        return None
    t, m, cur = [], [], None
    with open(logfile) as fh:
        for line in fh:
            if line.startswith("Time ="):
                cur = float(line.split("=")[1].split()[0].rstrip("s"))
            elif "Current mass in system" in line and cur is not None:
                t.append(cur * 1e3); m.append(float(line.split("=")[1]))
    return (np.array(t), np.array(m)) if t else None


# ----------------------------------------------------------------------------- metrics
def box_blur(img, w_px):
    """uniform (box) filter of width w_px pixels in both directions = parcel volume spread over one CFD cell.
    A box (not Gaussian) kernel avoids long tails that would push the PLV threshold beyond the last parcel."""
    w = int(round(w_px))
    if w <= 1:
        return img
    for ax in (0, 1):
        c = np.cumsum(np.pad(img, [(w, w) if i == ax else (0, 0) for i in (0, 1)]), axis=ax)
        hi = c.take(range(w + w // 2, w + w // 2 + img.shape[ax]), axis=ax)
        lo = c.take(range(w // 2, w // 2 + img.shape[ax]), axis=ax)
        img = (hi - lo) / w
    return img


def plv_map(xyz, vol_mm3, view_deg, pix, sigma, half_width=60.0, depth=110.0):
    """Projected liquid volume [mm3/mm2] on a (lateral s, axial a) grid; view along direction view_deg in x-y."""
    th = np.radians(view_deg)
    s = (-xyz[:, 0] * np.sin(th) + xyz[:, 1] * np.cos(th)) * 1e3   # lateral coordinate perpendicular to view
    a = -xyz[:, 2] * 1e3                                            # axial distance from tip [mm]
    se = np.arange(-half_width, half_width + pix, pix)
    ae = np.arange(0.0, depth + pix, pix)
    H, _, _ = np.histogram2d(a, s, bins=[ae, se], weights=vol_mm3)
    return box_blur(H, sigma / pix) / pix ** 2, ae, se


def liquid_metrics(vtk, view_deg, pix, sigma, thr):
    P = read_vtk_lagrangian(vtk)
    xyz = P["points"]
    if len(xyz) == 0:
        return dict(n=0, mass=0.0)
    d = P.get("d")
    nP = P.get("nParticle", np.ones(len(xyz)))
    rho = P.get("rho", np.full(len(xyz), 633.0))
    vol = nP * np.pi / 6.0 * d ** 3                      # m3 liquid per parcel
    mass = vol * rho
    a = -xyz[:, 2] * 1e3
    o = np.argsort(a); cm = np.cumsum(mass[o]) / mass.sum()
    plv, ae, se = plv_map(xyz, vol * 1e9, view_deg, pix, sigma)
    rows = np.where((plv >= thr).any(axis=1))[0]
    pen = ae[rows.max() + 1] if len(rows) else 0.0
    widths = {}
    for z in (15, 25, 35):
        r = plv[min(int(z / pix), plv.shape[0] - 1)] >= thr
        widths[z] = (se[np.where(r)[0].max() + 1] - se[np.where(r)[0].min()]) if r.any() else 0.0
    return dict(n=len(xyz), mass=mass.sum(), pen_plv=pen, tip=a.max(),
                pen95=a[o][np.searchsorted(cm, 0.95)], pen99=a[o][min(np.searchsorted(cm, 0.99), len(a) - 1)],
                w15=widths[15], w25=widths[25], w35=widths[35],
                d32=(nP * d ** 3).sum() / (nP * d ** 2).sum() * 1e6)


def vapour_pen(case, tdir, cz, ythr):
    Y = read_foam_internal(os.path.join(case, tdir, "IC8H18"))
    if Y is None or cz is None:
        return np.nan
    if np.isscalar(Y):
        return 0.0 if Y < ythr else -cz.min() * 1e3
    m = Y >= ythr
    return -cz[m].min() * 1e3 if m.any() else 0.0


def rmsd(t_sim, y_sim, t_exp, y_exp, t0=0.1, t1=1.5):
    sel = (t_exp >= max(t0, t_sim.min())) & (t_exp <= min(t1, t_sim.max())) & (y_exp > 0)
    if sel.sum() < 3:
        return np.nan, np.nan
    ys = np.interp(t_exp[sel], t_sim, y_sim)
    r = np.sqrt(np.mean((ys - y_exp[sel]) ** 2))
    return r, 100 * r / np.mean(y_exp[sel])


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("--ecn", default=os.path.join(HERE, "ecn_G3_penetration.csv"))
    ap.add_argument("--label", default=None)
    ap.add_argument("--view", type=float, default=0.0, help="line-of-sight azimuth [deg] (0 = looking along +y)")
    ap.add_argument("--pix", type=float, default=0.1, help="PLV pixel size [mm]")
    ap.add_argument("--sigma", type=float, default=1.0, help="box-kernel width for parcel smearing [mm] (= CFD cell size)")
    ap.add_argument("--plv-thr", type=float, default=2e-4, help="liquid PLV threshold [mm3/mm2]")
    ap.add_argument("--y-thr", type=float, default=1e-3, help="vapour mass-fraction threshold")
    ap.add_argument("--ref", nargs="*", default=[], help="digitised reference curves: csv time_ms,value  (name = legend)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    case = os.path.abspath(a.case)
    label = a.label or os.path.basename(case)
    out = a.out or os.path.join(case, "validation")
    os.makedirs(out, exist_ok=True)

    vtks = sorted(glob.glob(os.path.join(case, "VTK", "lagrangian", "*", "*.vtk")), key=time_of)
    if not vtks:
        sys.exit("no VTK/lagrangian/*/*.vtk found - run: foamToVTK -ascii -useTimeName -noInternal")
    cz = read_foam_internal(os.path.join(case, "0", "Cz"))
    if cz is None:
        C = read_foam_internal(os.path.join(case, "0", "C"))
        cz = C[:, 2] if C is not None else None
    if cz is None:
        print("warning: no 0/Cz - vapour penetration skipped (run foamPostProcess -func writeCellCentres -time 0)")
    tdirs = {}
    for dname in os.listdir(case):
        try:
            tdirs[round(float(dname), 9)] = dname
        except ValueError:
            pass

    rows = []
    for v in vtks:
        t = time_of(v)
        if t is None:
            continue
        L = liquid_metrics(v, a.view, a.pix, a.sigma, a.plv_thr)
        td = tdirs.get(round(t, 9))
        vp = vapour_pen(case, td, cz, a.y_thr) if td else np.nan
        rows.append([t * 1e3, L.get("pen_plv", 0), vp, max(vp, L.get("pen_plv", 0)) if np.isfinite(vp) else np.nan,
                     L["mass"] / M_INJ, L.get("tip", 0), L.get("pen95", 0), L.get("pen99", 0),
                     L.get("w15", 0), L.get("w25", 0), L.get("w35", 0), L.get("d32", np.nan), L["n"]])
        print(f"t={t*1e3:6.3f} ms  liq(PLV)={rows[-1][1]:6.2f}  vap={vp:6.2f}  mass ratio={rows[-1][4]:5.3f}  n={L['n']}")
    S = np.array(rows)
    hdr = "time_ms,liq_pen_plv_mm,vap_pen_Y_mm,vap_pen_combined_mm,liquid_mass_ratio,tip_mm,pen95_mm,pen99_mm,w15_mm,w25_mm,w35_mm,D32_um,nParcels"
    np.savetxt(os.path.join(out, "sim_metrics.csv"), S, delimiter=",", header=hdr, comments="", fmt="%.5g")

    E = np.genfromtxt(a.ecn, delimiter=",", names=True)
    te = E["time_ms"]
    r_l = rmsd(S[:, 0], S[:, 1], te, E["liq_len_lowthr_avgmap_mm"])
    r_v = rmsd(S[:, 0], S[:, 3], te, E["vap_len_indiv_mm"]) if np.isfinite(S[:, 3]).any() else (np.nan, np.nan)
    i_eoi = np.argmin(abs(S[:, 0] - EOI))
    summary = (f"{label}\n  liquid penetration RMSD (0.1-1.5 ms) = {r_l[0]:.2f} mm ({r_l[1]:.1f} % of mean)\n"
               f"  vapour penetration RMSD (0.1-1.5 ms) = {r_v[0]:.2f} mm ({r_v[1]:.1f} % of mean)\n"
               f"  liquid mass ratio at EOI ({S[i_eoi,0]:.2f} ms) = {S[i_eoi,4]:.3f}   (Lien G3: ~0.45)\n")
    print(summary)
    open(os.path.join(out, "summary.txt"), "w").write(summary)

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    refs = [(os.path.splitext(os.path.basename(r))[0], np.loadtxt(r, delimiter=",", ndmin=2)) for r in a.ref]

    # --- Fig. 11-style: penetration
    fig, ax = plt.subplots(figsize=(7.5, 5), dpi=150)
    m = te <= 1.6
    ax.plot(te[m], E["liq_len_lowthr_avgmap_mm"][m], "k-", lw=2, label="ECN liquid (PLV 2e-4)")
    ax.fill_between(te[m], (E["liq_len_lowthr_indiv_mm"] - E["liq_len_lowthr_std_mm"])[m],
                    (E["liq_len_lowthr_indiv_mm"] + E["liq_len_lowthr_std_mm"])[m], color="k", alpha=0.12, lw=0)
    ax.plot(te[m], E["vap_len_indiv_mm"][m], "k--", lw=2, label="ECN vapour (schlieren)")
    ax.plot(S[:, 0], S[:, 1], "o-", color="tab:blue", ms=4, label=f"sim liquid (PLV) RMSD {r_l[1]:.1f}%")
    ax.plot(S[:, 0], S[:, 3], "s--", color="tab:red", ms=4, label=f"sim vapour (Y>1e-3) RMSD {r_v[1]:.1f}%")
    ax.plot(S[:, 0], S[:, 5], ":", color="tab:green", lw=1.2, label="sim max parcel tip")
    for name, R in refs:
        ax.plot(R[:, 0], R[:, 1], "-", lw=1.2, label=name)
    ax.axvline(EOI, color="gray", ls=":", lw=1)
    ax.set(xlim=(0, 1.6), ylim=(0, 85), xlabel="time aSOI [ms]", ylabel="axial penetration [mm]",
           title=f"ECN Spray G3 - {label}")
    ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig11_penetration.png"))

    # --- Fig. 14-style: liquid mass ratio
    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)
    lm = parse_log_mass(os.path.join(case, "log.foamRun"))
    if lm is not None:
        ax.plot(lm[0], lm[1] / M_INJ, "-", color="tab:blue", lw=1.5, label="sim (log, every step)")
    ax.plot(S[:, 0], S[:, 4], "o", color="tab:blue", ms=4, label="sim (parcels at write times)")
    for name, R in refs:
        if "mass" in name.lower():
            ax.plot(R[:, 0], R[:, 1], "k-", lw=1.5, label=name)
    ax.axvline(EOI, color="gray", ls=":", lw=1); ax.axhline(0.45, color="gray", ls="--", lw=0.8)
    ax.text(EOI + 0.02, 0.47, "Lien G3 @ EOI ~ 0.45", fontsize=8, color="gray")
    ax.set(xlim=(0, 1.6), ylim=(0, 1), xlabel="time aSOI [ms]", ylabel="liquid mass / injected mass [-]",
           title=f"Liquid mass ratio - {label}")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig14_liquid_mass_ratio.png"))
    print(f"wrote {out}/sim_metrics.csv, summary.txt, fig11_penetration.png, fig14_liquid_mass_ratio.png")


if __name__ == "__main__":
    main()
