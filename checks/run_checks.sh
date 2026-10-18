#!/usr/bin/env bash
# RUN_CHECKS  Runs check scripts headless and files their logs/CSVs under
# analysis_outputs/<stage>/checks_<scorer>/.
#   usage: checks/run_checks.sh <stage> "<scorers>" "<scripts>"
#   e.g.   checks/run_checks.sh after_bugfixes "orig patched" "t1_lambda_recovery t3_parameter_sweep t4_granule_bootstrap"
# t2 does not score, so run it with a single scorer ("orig").
set -u
stage=$1; scorers=${2:-"orig patched"}; scripts=${3:-"t1_lambda_recovery t2_false_positives_and_bias t3_parameter_sweep t4_granule_bootstrap"}
M="/c/Program Files/MATLAB/R2025b/bin/matlab.exe"
cd "$(dirname "$0")"
for sc in $scorers; do
  out="../analysis_outputs/$stage/checks_$sc"; mkdir -p "$out"
  for s in $scripts; do
    echo "=== $stage $sc $s start $(date +%T)"
    GM_SCORER=$sc "$M" -batch "$s" > "$out/$s.log" 2>&1; echo "=== $s exit $? $(date +%T)"
    mv -f ${s}_*.csv "$out/" 2>/dev/null || true
  done
done
exit 0
