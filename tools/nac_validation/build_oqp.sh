#!/bin/bash
# Incremental rebuild + install of liboqp.so after a Fortran change.
set -e
export PATH=/opt/soft/install/GCCcore/12.3.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/soft/install/GCCcore/12.3.0/lib64:$LD_LIBRARY_PATH
export CFLAGS="-g0"
cd /bighome/alireza/openqp-nac
ninja -C build install 2>&1 | tail -25
echo "BUILD_EXIT=${PIPESTATUS[0]}"
ls -la lib/liboqp.so
