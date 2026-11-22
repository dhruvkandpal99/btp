#!/bin/sh
# Which breakup models survive the smoke test with Lien's small droplets?
# Run from case_G3:  sh tools/bisect_breakup.sh
cd "$(dirname "$0")/.." || exit 1
try() {   # try <name> <generator flags> [env]
    printf "%-26s" "$1"
    if env $3 CLOUD_ARGS="$2" sh Allrun smoke > log.bisect_$1 2>&1; then
        if grep -qi "nan" log.smoke; then echo "PASS but NaN in log"; else echo PASS; fi
    else echo FAIL; fi
}
try PilchErdman        "--breakup PilchErdman"
try ReitzDiwakar       "--breakup ReitzDiwakar"
try TAB                "--breakup TAB --oscillation yes"
try ETAB               "--breakup ETAB --oscillation yes"
# KHRT with FPE trapping off: tells us whether KHRT produces real NaNs or just trips the trap
try KHRT_noFPEtrap     "--breakup ReitzKHRT" "FOAM_SIGFPE=false"
python3 tools/make_cloudProperties.py > /dev/null   # restore defaults
echo "details: log.bisect_*"
