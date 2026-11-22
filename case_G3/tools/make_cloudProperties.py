#!/usr/bin/env python3
"""
Generate constant/cloudProperties for the ECN Spray G3 case (OpenFOAM v14, sprayCloud).

All the knobs you will touch during calibration are in PARAMS below (or on the
command line), so a sweep is one line, e.g.

    python3 tools/make_cloudProperties.py --Cd 0.73
    python3 tools/make_cloudProperties.py --sizes blob          # Lien Case 2
    python3 tools/make_cloudProperties.py --thetaOuter 10 --B1 10

Geometry convention (same as Lien et al. 2024 / ECN):
    injector tip at the origin, injector axis pointing down -z,
    8 plumes, each tilted `plumeAngle` degrees from the axis,
    plume azimuths phi_k = phi0 + 45 k.
"""
import argparse, math, os, sys

PARAMS = dict(
    # ---- operating condition (ECN Spray G3) ----
    T0=363.0,              # fuel temperature [K]
    massTotalAll=10.0e-6,  # total injected mass, all 8 holes [kg]  (10 mg)
    duration=7.8e-4,       # hydraulic injection duration [s]
    # trapezoidal rate shape (relative values; OpenFOAM rescales to massTotal)
    tOpen=5.0e-5,          # end of needle-opening ramp [s]
    tClose=7.0e-4,         # start of needle-closing ramp [s]

    # ---- nozzle / plume geometry ----
    dOuter=1.65e-4,        # hole exit diameter [m] (165 um)
    Cd=0.60,               # discharge coeff. -> U = mdot/(Cd*rho*A).  0.52 ~ 250 m/s, 0.60 ~ 215 m/s, 0.73 ~ 177 m/s
    plumeAngle=35.0,       # plume direction from injector axis [deg]  (Lien Table 3)
    thetaOuter=15.0,       # cone HALF-angle per plume [deg] (Lien: 30 deg full cone)
    phi0=22.5,             # azimuth of first plume [deg]
    holeRadius=0.5e-3,     # radial position of hole exits [m]
    tipZ=-0.5e-3,          # axial position of hole exits [m] (just inside first cell row)

    # ---- parcels ----
    parcelsPerSecond=1.0e7,  # per hole (~7800 parcels/hole over 0.78 ms)

    # ---- breakup: Reitz KH-RT, constants = Lien Eqs. 33-35 evaluated at We~15 ----
    B0=0.55, B1=7.0, CRT=0.42, Ctau=1.0, msLimit=0.2, WeberLimit=6.0,

    # ---- sub-model switches (debugging / ablation) ----
    drag="sphereDrag",                     # or distortedSphereDrag (only matters with a breakup model that solves the oscillation eq.)
    dispersion="stochasticDispersionRAS",  # or none
    breakup="PilchErdman",                 # none | PilchErdman | ReitzKHRT | ReitzDiwakar | TAB | ETAB
                                           # v14 smoke tests with Lien's 1.6-23 um droplets: only none and PilchErdman
                                           # pass; ReitzKHRT/ReitzDiwakar/TAB/ETAB crash on step 1 (FPE / NaN)
    oscillation="no",                      # solveOscillationEq (distortedSphereDrag needs yes)

    # ---- initial droplet size ----
    sizes="lien",          # "lien" (Case 1), "lienTable" (exact tabulated CDF), "blob" (Case 2)
)

# Lien et al. 2024, Eq. (36):  d(x) = 1.569 exp(1.823 x) + 5.185e-13 exp(30.88 x)  [um], x ~ U(0,1)
# (the paper prints "5.185"; the 1e-13 factor is needed to reproduce their Fig. 3 CDF - checked point by point)
# -> number-based CDF.  Truncated Rosin-Rammler fit (Q=0): d=3.27 um, n=0.99, [1.569, 23.07] um, max CDF error 0.017
LIEN_RR = dict(min=1.569e-6, max=2.307e-5, d=3.27e-6, n=0.99)
LIEN_TABLE = [  # (diameter [m], number CDF) sampled from Eq. 36
    (1.5690e-06, 0.000), (1.7187e-06, 0.050), (1.8828e-06, 0.100), (2.2593e-06, 0.200),
    (2.7111e-06, 0.300), (3.2532e-06, 0.400), (3.9038e-06, 0.500), (4.6845e-06, 0.600),
    (5.6224e-06, 0.700), (6.7730e-06, 0.800), (7.5190e-06, 0.850), (8.7032e-06, 0.900),
    (1.0087e-05, 0.930), (1.1719e-05, 0.950), (1.4486e-05, 0.970), (1.6568e-05, 0.980),
    (1.9347e-05, 0.990), (2.1072e-05, 0.995), (2.3071e-05, 1.000),
]


def size_block(kind, ind="            "):
    if kind == "lien":
        r = LIEN_RR
        return (f"{ind}sizeDistribution\n{ind}{{\n"
                f"{ind}    // Lien et al. 2024 Eq.36 (number-based) as truncated Rosin-Rammler, Q 0 = number distribution\n"
                f"{ind}    type        RosinRammler;\n{ind}    Q           0;\n"
                f"{ind}    min         {r['min']:.4e};\n{ind}    max         {r['max']:.4e};\n"
                f"{ind}    d           {r['d']:.4e};\n{ind}    n           {r['n']};\n{ind}}}\n")
    if kind == "lienTable":
        rows = "\n".join(f"{ind}        ({d:.4e} {c:.3f})" for d, c in LIEN_TABLE)
        return (f"{ind}sizeDistribution\n{ind}{{\n"
                f"{ind}    // exact Eq.36 CDF; use only if your v14 build has the tabulatedCumulative distribution\n"
                f"{ind}    type        tabulatedCumulative;\n{ind}    Q           0;\n"
                f"{ind}    distribution\n{ind}    (\n{rows}\n{ind}    );\n{ind}}}\n")
    if kind == "blob":
        return (f"{ind}sizeDistribution\n{ind}{{\n"
                f"{ind}    // blob method: droplets = nozzle diameter (Lien Case 2)\n"
                f"{ind}    type        fixedValue;\n{ind}    value       1.65e-4;\n{ind}}}\n")
    sys.exit(f"unknown sizes option {kind}")


def hole_block(k, p):
    phi = math.radians(p["phi0"] + 45.0 * k)
    a = math.radians(p["plumeAngle"])
    pos = (p["holeRadius"] * math.cos(phi), p["holeRadius"] * math.sin(phi), p["tipZ"])
    dirn = (math.sin(a) * math.cos(phi), math.sin(a) * math.sin(phi), -math.cos(a))
    f = lambda v: "(" + " ".join(f"{c: .6f}" for c in v) + ")"
    T, to, tc = p["duration"], p["tOpen"], p["tClose"]
    return f"""        hole{k + 1}
        {{
            type            coneInjection;
            SOI             0;
            massTotal       {p['massTotalAll'] / 8:.6e};
            uniformParcelSize volume;
            injectionMethod disc;
            flowType        flowRateAndDischarge;
            dInner          0;
            dOuter          {p['dOuter']:.4e};
            duration        {T:.4e};
            position        {f(pos)};
            direction       {f(dirn)};
            parcelsPerSecond {p['parcelsPerSecond']:.4e};
            flowRateProfile table
            (
                (0          0)
                ({to:.4e} 1)
                ({tc:.4e} 1)
                ({T:.4e} 0)
            );
            Cd              constant {p['Cd']};
            thetaInner      constant 0;
            thetaOuter      constant {p['thetaOuter']};
{size_block(p['sizes'])}        }}
"""


def build(p):
    holes = "".join(hole_block(k, p) for k in range(8))
    rho_l = 633.0  # iso-octane at 363 K (only a placeholder: sprayCloud overwrites rho/Cp from the IC8H18 liquid model)
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  14
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      cloudProperties;
}}
// GENERATED by tools/make_cloudProperties.py - edit the script, not this file.
// sizes={p['sizes']}  Cd={p['Cd']}  thetaOuter={p['thetaOuter']}  plumeAngle={p['plumeAngle']}
// KHRT: B0={p['B0']} B1={p['B1']} CRT={p['CRT']}
// drag={p['drag']} dispersion={p['dispersion']} breakup={p['breakup']} oscillation={p['oscillation']}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

type        sprayCloud;

cpuLoad     true;

solution
{{
    coupled         true;
    transient       yes;
    cellValueSourceCorrection on;
    maxCo           0.3;

    sourceTerms
    {{
        schemes
        {{
            rho             explicit 1;
            U               explicit 1;
            Yi              explicit 1;
            h               explicit 1;
            radiation       explicit 1;
        }}
    }}

    interpolationSchemes
    {{
        rho             cell;
        U               cellPoint;
        mu              cell;
        T               cell;
        Cp              cell;
        kappa           cell;
        p               cell;
    }}

    integrationSchemes
    {{
        U               Euler;
        T               analytical;
    }}
}}

constantProperties
{{
    T0              {p['T0']};
    rho0            {rho_l};     // iso-octane @ 363 K (overwritten by IC8H18 liquid properties)
    Cp0             2300;    // iso-octane @ 363 K [J/kg/K] (overwritten by IC8H18 liquid properties)
    constantVolume  false;
}}

subModels
{{
    particleForces
    {{
        {p['drag']};
    }}

    injectionModels
    {{
{holes}    }}

    dispersionModel         {p['dispersion']};

    patchInteractionModel   standardWallInteraction;

    heatTransferModel       RanzMarshall;

    compositionModel        singlePhaseMixture;

    phaseChangeModel        liquidEvaporationBoil;

    surfaceFilmModel        none;

    atomisationModel        none;

    breakupModel            {p['breakup']};

    stochasticCollisionModel none;

    radiation               off;

    standardWallInteraction
    {{
        type                rebound;
    }}

    RanzMarshall
    {{
        BirdCorrection      true;
    }}

    singlePhaseMixture
    {{
        phases
        (
            liquid
            {{
                IC8H18              1;
            }}
        );
    }}

    liquidEvaporationBoil
    {{
        enthalpyTransfer    enthalpyDifference;
        activeLiquids       ( IC8H18 );
    }}

    // breakup model coefficients (only the one selected by breakupModel is used)
    ReitzKHRT
    {{
        solveOscillationEq  {p['oscillation']};
        B0                  {p['B0']};
        B1                  {p['B1']};
        Ctau                {p['Ctau']};
        CRT                 {p['CRT']};
        msLimit             {p['msLimit']};
        WeberLimit          {p['WeberLimit']};
    }}

    PilchErdman
    {{
        solveOscillationEq  {p['oscillation']};
        B1                  0.375;
        B2                  0.2274;
    }}

    ReitzDiwakar
    {{
        solveOscillationEq  {p['oscillation']};
        Cbag                6;
        cb                  0.785;
        Cstrip              0.5;
        Cs                  10;
    }}

    TAB
    {{
        solveOscillationEq  yes;
        y0                  0;
        yDot0               0;
        SMDCalculationMethod method2;
    }}

    ETAB
    {{
        solveOscillationEq  yes;
        Cmu                 10;
        Comega              8;
        WeCrit              12;
        k1                  0.2;
        k2                  0.2;
        WeTransition        100;
    }}
}}

cloudFunctions
{{}}

// ************************************************************************* //
"""


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k, v in PARAMS.items():
        ap.add_argument("--" + k, type=type(v), default=v)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "constant", "cloudProperties"))
    a = vars(ap.parse_args())
    out = a.pop("out")
    if a["drag"] == "distortedSphereDrag" and a["oscillation"] != "yes":
        sys.exit("distortedSphereDrag needs --oscillation yes")
    with open(out, "w") as fh:
        fh.write(build(a))
    print(f"wrote {os.path.normpath(out)}  (sizes={a['sizes']}, Cd={a['Cd']}, thetaOuter={a['thetaOuter']})")
