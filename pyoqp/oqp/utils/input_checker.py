"""Input value checker"""

import os
import multiprocessing

def check_input_values(config):
    runtype = config['input']['runtype']

    check_func = {
        'energy': check_energy_input,
        'grad': check_grad_input,
        'nac': check_nac_input,
        'bp': check_nac_input,
        'soc': check_soc_input,
        'optimize': check_optimize_input,
        'meci': check_optimize_input,
        'mecp': check_optimize_input,
        'mep': check_optimize_input,
        'ts': check_optimize_input,
        'neb': check_optimize_input,
        'hess': check_hess_input,
        'nacme': check_nacme_input,
        'prop': skip_check,
    }

    info = f'\nPyOQP checking input\n\n[input] runtype={runtype}'
    check_func[runtype](config, info)

def skip_check(config, info):
    pass

def not_available(func_name, info):
    exit(f'{info}\nPyOQP: runtype {func_name} is not available yet\n')

def check_scf_input(config, info):
    scf_type = config['scf']['type']
    scf_mult = config['scf']['multiplicity']

    info += f'\n[scf] type={scf_type}'
    info += f'\n[scf] multiplicity={scf_mult}'

    if scf_mult < 1:
        exit(f'{info}\nPyOQP:scf multiplicity must be greater than 0\n')

    if scf_mult > 1 and scf_type == 'rhf':
        exit(f'{info}\nPyOQP: scf multiplicity {scf_mult} does not match rhf, choose uhf or rohf\n')

def check_tdhf_input(config, info):
    check_scf_input(config, info)

    scf_type = config['scf']['type']
    scf_mult = config['scf']['multiplicity']
    td_type = config['tdhf']['type']
    td_mult = config['tdhf']['multiplicity']

    info += f'\n[scf] type={scf_type}'
    info += f'\n[scf] multiplicity={scf_mult}'
    info += f'\n[tdhf] type={td_type}'
    info += f'\n[tdhf] multiplicity={td_mult}'

    if td_type in ['rpa', 'tda'] and scf_mult != td_mult:
        print(f'{info}\nPyOQP: Caution! tdhf type {td_type} multiplicity {td_mult} is not equal to scf multiplicity {scf_mult}\n')

    if td_type in ['sf', 'mrsf'] and scf_mult == td_mult:
        print(f'{info}\nPyOQP: Caution! tdhf type {td_type} multiplicity {td_mult} is equal to scf multiplicity {scf_mult}\n')

    if td_type in ['mrsf', 'sf'] and scf_type != 'rohf':
        exit(f'{info}\nPyOQP: tdhf type {td_type} cannot use {scf_type} orbitals, choose rohf in scf type\n')

def check_energy_input(config, info):
    method = config['input']['method']

    check_func = {
        'hf': check_scf_input,
        'tdhf': check_tdhf_input,
    }

    info += f'\n[input] method={method}'
    check_func[method](config, info)

def check_grad_input(config, info):

    check_energy_input(config, info)

    method = config['input']['method']
    td_type = config['tdhf']['type']
    grad = config['properties']['grad']

    info += f'\n[input] method={method}'
    info += f'\n[tdhf] type={td_type}'
    info += f'\n[properties] grad={grad}'

    if method == 'hf' and sum(grad) > 0:
        exit(f'{info}\nPyOQP: scf cannot compute gradient for state > 0\n')

    if method == 'tdhf' and 0 in grad:
        exit(f'{info}\nPyOQP: tdhf type {td_type} cannot compute gradient for state 0 \n')

def check_optimize_input(config, info):
    runtype = config['input']['runtype']
    check_energy_input(config, info)
    method = config['input']['method']
    istate = config['optimize']['istate']
    jstate = config['optimize']['jstate']
    imult = config['optimize']['imult']
    jmult = config['optimize']['jmult']

    info += f'\n[input] method={method}'
    info += f'\n[optimize] istate={istate}'

    if method == 'hf' and istate > 0:
        exit(f'{info}\nPyOQP: scf cannot compute gradient for istate > 0\n')

    if method == 'tdhf' and istate == 0:
        exit(f'{info}\nPyOQP: tdhf cannot compute gradient for istate = 0 \n')

    info += f'\n[optimize] jstate={jstate}'
    if runtype == 'meci' and jstate - istate < 1:
        exit(f'{info}\nPyOQP: meci optimization require jstate >= istate + 1 \n')

    if runtype == 'mecp' and imult == jmult:
        exit(f'{info}\nPyOQP: mecp optimization require imult != jmult \n')

    lib = config['optimize']['lib']
    info += f'\n[optimize] lib={lib}'

    if lib not in ['scipy', 'dlfind']:
        exit(f'{info}\nPyOQP: cannot recognize the geometry optimization library {lib} \n')

    if lib == 'dlfind':
        check_dlfind_input(config, info)

def check_dlfind_input(config, info):
    runtype = config['input']['runtype']
    icoord = config['dlfind']['icoord']
    iopt = config['dlfind']['iopt']
    ims = config['dlfind']['ims']

    info += f'\n[dlfind] icoord={icoord}'
    info += f'\n[dlfind] iopt={iopt}'
    info += f'\n[dlfind] ims={ims}'

    check_func = {
        'optimize': check_dlfind_min_input,
        'meci': check_dlfind_meci_input,
        'ts': check_dlfind_ts_input,
    }

    check_func[runtype](icoord, iopt, ims, info)

def check_dlfind_min_input(icoord, iopt, ims, info):
    if icoord not in [0, 1, 2, 3, 4]:
        exit(f'{info}\nPyOQP: dl-find single state optimization only support icoord 0–4, but found {icoord}')

    if iopt not in [0, 1, 2, 3]:
        exit(f'{info}\nPyOQP: dl-find single state optimization only support iopt 0–3, but found {iopt}')

    if ims not in [0]:
        exit(f'{info}\nPyOQP: dl-find single state optimization only support ims 0, but found {ims}')

def check_dlfind_meci_input(icoord, iopt, ims, info):
    if iopt not in [0, 1, 2, 3]:
        exit(f'{info}\nPyOQP: dl-find meci optimization only support iopt 0–3, but found {iopt}')

    if ims not in [1, 2, 3]:
        exit(f'{info}\nPyOQP: dl-find meci optimization only support ims 1-3, but found {ims}')

    if ims == 3 and icoord not in [10, 11, 12, 13, 14]:
        exit(f'{info}\nPyOQP: dl-find meci optimization {ims} only support icoord 10–14, but found {icoord}')

    if ims in [1, 2] and icoord not in [0, 1, 2, 3, 4]:
        exit(f'{info}\nPyOQP: dl-find meci optimization {ims} only support icoord 0–4, but found {icoord}')

def check_dlfind_ts_input(icoord, iopt, ims, info):
    if iopt not in [0, 1, 2, 3]:
        exit(f'{info}\nPyOQP: dl-find meci optimization only support iopt 0–3, but found {iopt}')

    if ims not in [1, 2, 3]:
        exit(f'{info}\nPyOQP: dl-find meci optimization only support ims 1-3, but found {ims}')

    if ims == 3 and icoord not in [10, 11, 12, 13, 14]:
        exit(f'{info}\nPyOQP: dl-find meci optimization {ims} only support icoord 10–14, but found {icoord}')

    if ims in [1, 2] and icoord not in [0, 1, 2, 3, 4]:
        exit(f'{info}\nPyOQP: dl-find meci optimization {ims} only support icoord 0–4, but found {icoord}')


def check_nac_input(config, info):
    check_energy_input(config, info)
    check_energy_input(config, info)
    method = config['input']['method']
    td = config['tdhf']['type']
    nproc = config['nac']['nproc']
    ncpu = multiprocessing.cpu_count()

    try:
        omp = int(os.environ['OMP_NUM_THREADS'])
    except KeyError:
        omp = 0

    info += f'\n[input] method={method}'

    if method == 'hf':
        exit(f'{info}\nPyOQP: scf cannot compute nac for hf calculations\n')

    if method == 'tdhf' and td != 'mrsf':
        exit(f'{info}\nPyOQP: tdhf cannot compute nac for non-mrsf calculations\n')

        info += f'\n[hess] nproc={nproc}'

    if omp > 0 and omp*nproc > ncpu:
        info += f'\nOMP_NUM_THREADS={omp}'
        exit(f'{info}\nPyOQP: nac requested {nproc * omp} cpu, exceeding available {ncpu} cpus')

    if omp == 0 and nproc > 1:
        info += f'\nOMP_NUM_THREADS not set, using all threads'
        exit(f'{info}\nPyOQP: nac requested {nproc * ncpu}, exceeding available {ncpu} cpus')

    if nproc < 1:
        exit(f'{info}\nPyOQP: nac requested {nproc} process is illegal ')

def check_soc_input(config, info):
    not_available('soc', info)


def check_neb_input(config, info):
    not_available('neb', info)


def check_hess_input(config, info):
    check_energy_input(config, info)
    method = config['input']['method']
    state = config['hess']['state']
    nproc = config['hess']['nproc']
    restart = config['hess']['restart']
    ncpu = multiprocessing.cpu_count()

    try:
        omp = int(os.environ['OMP_NUM_THREADS'])
    except KeyError:
        omp = 0

    info += f'\n[input] method={method}'
    info += f'\n[hess] state={state}'

    if method == 'hf' and state > 0:
        exit(f'{info}\nPyOQP: scf cannot compute hessian for state > 0\n')

    if method == 'tdhf' and state == 0:
        exit(f'{info}\nPyOQP: tdhf cannot compute hessian for state = 0 \n')

    if restart != 'read':
        info += f'\n[hess] nproc={nproc}'

        if omp > 0 and omp*nproc > ncpu:
            info += f'\nOMP_NUM_THREADS={omp}'
            exit(f'{info}\nPyOQP: hessian requested {nproc * omp} cpu, exceeding available {ncpu} cpus')

        if omp == 0 and nproc > 1:
            info += f'\nOMP_NUM_THREADS not set, using all threads'
            exit(f'{info}\nPyOQP: hessian requested {nproc * ncpu}, exceeding available {ncpu} cpus')

        if nproc < 1:
            exit(f'{info}\nPyOQP: hessian requested {nproc} process is illegal ')


def check_nacme_input(config, info):
    check_energy_input(config, info)
