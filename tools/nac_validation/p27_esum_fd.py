"""FD self-test of the closed-form esum: 2e (fock_jk ref) and XC (dftexcor ref).
MRSF is ROHF+functional -- run on the production BHHLYP input only."""
import os, sys
import numpy as np, oqp, oqp.library
from oqp.pyoqp import Runner
INP = sys.argv[1]
os.environ['NAC_ESUM_FDTEST'] = '1'
r = Runner(input_file=INP, log='/bighome/alireza/openqp-nac/tools/nac_validation/p27.log')
r.run()
mol = r.mol
for (i, j) in [(1, 2), (1, 3), (2, 3)]:
    oqp.mrsf_nac_esum(mol, i, j)
    print(f'=== pair ({i},{j}) ===')
    for f in ('/tmp/nac_esum_fdtest.out', '/tmp/nac_esum_xcfd.out', '/tmp/nac_esum_1efd.out'):
        try:
            print(open(f).read())
        except Exception as e:
            print(f'  [{f}: {e}]')
