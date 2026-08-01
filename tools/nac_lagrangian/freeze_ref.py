"""Dump canonical numerical-NAC reference vectors (post sign-fix) to npz."""
import os
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
    np.savez(sys.argv[2], nacv=nacv, dcv=dcv, energies=E, flags=np.array(flags))
    print('saved', sys.argv[2], 'flags:', set(flags))


if __name__ == '__main__':
    main()
