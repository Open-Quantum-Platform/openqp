"""Dump canonical numerical-NAC reference vectors (post sign-fix) to npz."""
import hashlib
import os
from pathlib import Path
import sys


def main():
    import oqp
    from oqp.pyoqp import Runner
    from oqp.library.single_point import NAC
    import numpy as np
    inp = sys.argv[1]
    project = os.path.basename(inp).removesuffix('.inp')
    r = Runner(
        project=project,
        input_file=inp,
        log=os.path.abspath(inp.replace('.inp', '_frz.log')),
    )
    r.run()
    nac = NAC(r.mol)
    nacv, dcv, flags = nac.numerical_nac()
    E = np.array(r.mol.energies)
    input_path = Path(inp).resolve()
    np.savez(
        sys.argv[2],
        nacv=nacv,
        dcv=dcv,
        energies=E,
        flags=np.array(flags),
        input_path=str(input_path),
        input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
        scf_conv=float(r.mol.config['scf']['conv']),
        tdhf_conv=float(r.mol.config['tdhf']['conv']),
        displacement=float(r.mol.config['nac']['dx']),
        tlf=int(r.mol.config['tdhf']['tlf']),
    )
    print('saved', sys.argv[2], 'flags:', set(flags))


if __name__ == '__main__':
    main()
