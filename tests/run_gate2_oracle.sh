#!/usr/bin/env bash
# Disciplined Gate-2 validation driver.
#
# Enforces the validation contract: clean rebuild from source, freshly synced
# liboqp.so, primary PySCF oracle executed (NOT skipped) on the low-symmetry C1
# all-distinct-atom case, with build_rc / oracle_rc / skipped / failed / errored
# all reported. A skipped primary oracle is treated as RED.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "== Gate-2 oracle driver =="

# 0. remove stale artifacts so nothing can be mistaken for a fresh run
rm -f /tmp/hess_nuc_blocks.txt /tmp/hess1_selftest.out /tmp/hess1_sep.out

# 1. clean rebuild from source
ninja -C build oqp > /tmp/gate2_build.log 2>&1
build_rc=$?
echo "build_rc=$build_rc"
if [ $build_rc -ne 0 ]; then
  echo "BUILD FAILED:"; grep -iE "error|FAILED" /tmp/gate2_build.log | head
  echo "GATE2=RED (build)"; exit 1
fi

# 2. freshly sync the runtime library actually loaded by pyoqp (cffi dlopen)
cp build/source/liboqp.so lib/liboqp.so
echo "synced lib/liboqp.so (sha=$(sha256sum lib/liboqp.so | cut -c1-12))"

# 3. run the primary oracle; skip-as-failure via OQP_ORACLE_REQUIRED=1
OQP_ORACLE_REQUIRED=1 OPENQP_ROOT="$ROOT" OMP_NUM_THREADS=1 \
  python3 -m unittest -v tests.test_hess_nuc_oracle > /tmp/gate2_oracle.log 2>&1
oracle_rc=$?

# 4. parse unittest summary counts
line="$(grep -aE '^(OK|FAILED)' /tmp/gate2_oracle.log | tail -1)"
ran="$(grep -aE '^Ran [0-9]+ test' /tmp/gate2_oracle.log | tail -1)"
skipped=$(printf '%s' "$line" | grep -oE 'skipped=[0-9]+' | grep -oE '[0-9]+'); skipped=${skipped:-0}
failed=$(printf '%s'  "$line" | grep -oE 'failures=[0-9]+' | grep -oE '[0-9]+'); failed=${failed:-0}
errored=$(printf '%s' "$line" | grep -oE 'errors=[0-9]+'   | grep -oE '[0-9]+'); errored=${errored:-0}

echo "oracle_rc=$oracle_rc"
echo "$ran"
echo "skipped=$skipped failed=$failed errored=$errored"

# 5. verdict: green only if rc==0 AND nothing skipped/failed/errored
if [ "$oracle_rc" -eq 0 ] && [ "$skipped" -eq 0 ] && [ "$failed" -eq 0 ] && [ "$errored" -eq 0 ]; then
  echo "GATE2=GREEN"
  exit 0
fi
echo "---- oracle log tail ----"
grep -aE 'mismatch|not proven|parity|disagree|charge factor|base V_C|AssertionError|Error:|FAIL' /tmp/gate2_oracle.log | head -20
echo "GATE2=RED"
exit 1
