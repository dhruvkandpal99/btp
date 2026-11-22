#!/bin/sh
# Find which sub-model causes the floating point exception: runs the smoke test
# once per configuration and prints PASS/FAIL. Run from case_G3:  sh tools/bisect_smoke.sh
cd "$(dirname "$0")/.." || exit 1
OLDKH="--B0 0.61 --B1 55 --CRT 0.1"
try() {   # try <name> <generator flags>
    printf "%-34s" "$1"
    if CLOUD_ARGS="$2" sh Allrun smoke > log.bisect_$1 2>&1; then echo PASS; else echo FAIL; fi
}
try A_sphereDrag_noOscillation  "--breakup ReitzKHRT --drag sphereDrag --oscillation no"
try B_oldKHRTconstants          "--breakup ReitzKHRT $OLDKH"
try C_noBreakup                 "--drag sphereDrag --oscillation no --breakup none"
try D_noDispersion              "--breakup ReitzKHRT --dispersion none"
try E_oldKH_sphereDrag          "--breakup ReitzKHRT --drag sphereDrag --oscillation no $OLDKH"
python3 tools/make_cloudProperties.py > /dev/null   # restore defaults
echo "details: log.bisect_* (Allrun output) and log.smoke (last test)"
