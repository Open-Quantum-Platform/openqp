"""Public model-curvature configuration without an electronic-structure backend."""
import ast
import importlib.util
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


INPUT = load('_curvature_input', 'pyoqp/oqp/utils/oqp_input.py')
CHECK = load('_curvature_checker', 'pyoqp/oqp/utils/input_checker.py')
PARSER = load('_curvature_parser', 'pyoqp/oqp/utils/input_parser.py')
OPTIONS = dict(model_hessian='lindh', hessian_update='gpr',
               gpr_history=12, gpr_length_scale=0.3)


@pytest.mark.parametrize('driver', ['opt', 'ts'])
def test_concise_options_lower_and_roundtrip(driver):
    text = (f'rhf/sto-3g\n{driver}(S0,model_hessian=lindh,hessian_update=gpr,'
            'gpr_history=12,gpr_length_scale=0.3)\n'
            'geom="""H 0 0 0\nH 0 0 0.7"""')
    result = INPUT.resolve_oqp_text(text)
    assert result.legacy_config['oqp'] == {k: str(v) for k, v in OPTIONS.items()}
    assert INPUT.resolve_oqp_text(result.canonical_text).legacy_config == result.legacy_config


def test_sectioned_input_uses_real_schema_types_and_defaults():
    tree = ast.parse((ROOT / 'pyoqp/oqp/molecule/oqpdata.py').read_text())
    schema_node = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                       and any(isinstance(t, ast.Name) and t.id == 'OQP_CONFIG_SCHEMA'
                               for t in n.targets))
    oqp_node = next(v for k, v in zip(schema_node.keys, schema_node.values)
                    if ast.literal_eval(k) == 'oqp')
    schema = {'input': {'system': {'type': str, 'default': 'h2.xyz'}},
              'oqp': eval(compile(ast.Expression(oqp_node), '<schema>', 'eval'))}
    parser = PARSER.OQPConfigParser(schema=schema)
    defaults = parser.validate()['oqp']
    assert {k: defaults[k] for k in OPTIONS} == dict(
        model_hessian='auto', hessian_update='auto', gpr_history=8, gpr_length_scale=0.5)
    parser.read_string('[oqp]\n' + '\n'.join(f'{k}={v}' for k, v in OPTIONS.items()))
    assert {k: parser.validate()['oqp'][k] for k in OPTIONS} == OPTIONS


def report_for(options, runtype='ts', lib='oqp', qmmm=False):
    config = {'input': {'runtype': runtype, 'method': 'hf', 'qmmm_flag': qmmm},
              'optimize': {'lib': lib}, 'oqp': options}
    report = CHECK.CheckReport()
    CHECK._check_model_curvature(config, report)
    return report


@pytest.mark.parametrize('key,value', [
    ('model_hessian', 'unknown'), ('hessian_update', 'bfgs'),
    ('gpr_history', 2), ('gpr_history', 21), ('gpr_history', 3.5),
    ('gpr_history', True), ('gpr_history', 'nan'),
    ('gpr_length_scale', 0), ('gpr_length_scale', 10.1),
    ('gpr_length_scale', float('nan')), ('gpr_length_scale', float('inf')),
    ('gpr_length_scale', True),
])
def test_invalid_curvature_controls_are_rejected(key, value):
    report = report_for({**OPTIONS, key: value})
    assert not report.ok
    assert f'oqp.{key}' in report.to_text()


@pytest.mark.parametrize('kwargs', [
    {'runtype': 'energy'}, {'runtype': 'neb'}, {'runtype': 'meci'},
    {'runtype': 'mecp'}, {'runtype': 'tci'}, {'runtype': 'irc'},
    {'runtype': 'mep'}, {'lib': 'scipy'}, {'lib': 'geometric'}, {'qmmm': True},
])
def test_nondefault_controls_reject_unsupported_drivers(kwargs):
    assert not report_for(OPTIONS, **kwargs).ok
    assert report_for({}, **kwargs).ok


@pytest.mark.parametrize('extra', [dict(freeze='distance(1,2)'),
                                   dict(init_hessian='analytical'),
                                   dict(init_hessian='numerical')])
def test_nondefault_controls_reject_constraints_and_molecular_hessian(extra):
    assert not report_for({**OPTIONS, **extra}).ok


@pytest.mark.parametrize('options', [dict(gpr_history=12), dict(gpr_length_scale=0.2)])
def test_nondefault_gpr_parameters_require_gpr(options):
    assert not report_for(options).ok


@pytest.mark.parametrize('history,scale', [(3, 0.001), (20, 10), (8, 0.5)])
def test_boundaries_are_valid(history, scale):
    assert report_for({**OPTIONS, 'gpr_history': history, 'gpr_length_scale': scale}).ok


def test_small_example_resolves_endpoint_files():
    path = ROOT / 'examples/OPT/NH3_RHF-HF_QST3_LINDH_GPR_OQP.oqp'
    cfg = INPUT.resolve_oqp_text(path.read_text(), source_path=path).legacy_config
    assert cfg['oqp']['model_hessian'] == 'lindh'
    assert cfg['oqp']['hessian_update'] == 'gpr'
    assert Path(cfg['oqp']['ts_guess']).is_file()
    assert Path(cfg['oqp']['ts_product']).is_file()


def test_native_runner_forwards_curvature_options_and_logger():
    # Stop at engine construction: this tests the real runner without E/G calls.
    captured = {}
    class EngineReached(Exception):
        pass
    def engine(*args, **kwargs):
        captured.update(kwargs)
        raise EngineReached
    modules = {}
    def stub(name, **attrs):
        module = ModuleType(name)
        module.__dict__.update(attrs)
        modules[name] = module
    base = type('Optimizer', (), {})
    stub('oqp')
    stub('oqp.library')
    stub('oqp.library.libscipy', Optimizer=base, StateSpecificOpt=base, MECIOpt=base, MECPOpt=base)
    stub('oqp.library.baeka', BaekAState=base)
    stub('oqp.library.oqp_engine', OQPEngine=engine, parse_frozen_distance_spec=lambda value: [])
    stub('oqp.library.oqp_neb', NEB=base)
    stub('oqp.library.neb_utils', _read_xyz=None, kabsch_align=None, write_neb_xyz=None)
    stub('oqp.periodic_table', ELEMENTS_NAME={}, SYMBOL_MAP={})
    messages = []
    stub('oqp.utils')
    stub('oqp.utils.file_utils', dump_log=lambda mol, title: messages.append(title), dump_data=None)
    with patch.dict(sys.modules, modules):
        runtime = load('_curvature_runner', 'pyoqp/oqp/library/liboqp.py')
    runner = runtime.OQPOpt.__new__(runtime.OQPOpt)
    runner.mol = SimpleNamespace(config={'input': {'runtype': 'optimize'},
                                       'oqp': {**OPTIONS, 'coordsys': 'cart'}},
                                 get_atoms=lambda: [1, 1], get_mass=lambda: [1., 1.])
    runner.pre_coord = np.array([0., 0., 0., 0., 0., 1.4])
    runner.maxit = 3
    with pytest.raises(EngineReached):
        runner.optimize()
    assert {k: captured[k] for k in OPTIONS} == OPTIONS
    captured['logger']('GP diagnostic')
    assert messages[-1] == 'GP diagnostic'


@pytest.mark.parametrize('runtype', ['energy', 'tci', 'optimize', 'ts'])
def test_general_runtype_preflight_checks_curvature_controls(runtype):
    config = {'input': {'runtype': runtype, 'method': 'hf'},
              'optimize': {'lib': 'oqp'}, 'oqp': {'gpr_history': 21}}
    report = CHECK.CheckReport()
    CHECK._check_runtype(config, report)
    assert 'oqp.gpr_history' in report.to_text()


@pytest.mark.parametrize('controls', [dict(gpr_history=3.5), dict(gpr_history=True),
    dict(gpr_history=float('inf')), dict(gpr_history=float('nan')), dict(gpr_history=21),
    dict(gpr_length_scale=True), dict(gpr_length_scale=float('nan')),
    dict(gpr_length_scale=0), dict(gpr_length_scale=11)])
def test_raw_runner_controls_are_not_silently_coerced(controls):
    tree = ast.parse((ROOT / 'pyoqp/oqp/library/liboqp.py').read_text())
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == '_OQPRunner')
    method = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == '_oqp_config')
    namespace = {'np': np, 'parse_frozen_distance_spec': lambda value: []}
    exec(compile(ast.Module(body=[method], type_ignores=[]), '<runner>', 'exec'), namespace)
    runner = SimpleNamespace(mol=SimpleNamespace(config={'oqp': controls}))
    with pytest.raises(ValueError):
        namespace['_oqp_config'](runner)


@pytest.mark.parametrize('kwargs', [dict(runtype='energy'), dict(runtype='meci'),
    dict(runtype='mecp'), dict(runtype='neb'), dict(runtype='irc'), dict(runtype='mep'),
    dict(lib='scipy'), dict(lib='geometric'), dict(qmmm=True)])
def test_automatic_model_default_is_valid_outside_supported_lindh_searches(kwargs):
    assert report_for(dict(model_hessian='auto'), **kwargs).ok


@pytest.mark.parametrize('extra', [dict(freeze='distance(1,2)'), dict(init_hessian='analytical'),
                                   dict(init_hessian='numerical')])
def test_automatic_model_preserves_constraints_and_calculated_hessians(extra):
    assert report_for(dict(model_hessian='auto', **extra)).ok


@pytest.mark.parametrize('runtype,state,expected', [
    ('optimize', 0, 'auto'), ('optimize', 1, 'auto'), ('optimize', 3, 'auto'),
    ('ts', 0, 'auto'), ('ts', 1, 'auto'),
    ('meci', 1, 'constant'), ('mecp', 1, 'constant'), ('tci', 1, 'constant'),
    ('neb', 0, 'constant'), ('irc', 0, 'constant'), ('mep', 1, 'constant')])
def test_runner_default_distinguishes_state_specific_and_crossing_objectives(runtype, state, expected):
    tree = ast.parse((ROOT / 'pyoqp/oqp/library/liboqp.py').read_text())
    cls = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == '_OQPRunner')
    method = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == '_oqp_config')
    namespace = {'np': np, 'parse_frozen_distance_spec': lambda value: []}
    exec(compile(ast.Module(body=[method], type_ignores=[]), '<runner>', 'exec'), namespace)
    config = {'input': {'runtype': runtype}, 'optimize': {'istate': state}}
    # All these runner classes can use mode=min; state index does not select
    # the initial model, while the physical optimization objective must.
    runner = SimpleNamespace(mode='min', mol=SimpleNamespace(config=config))
    options = namespace['_oqp_config'](runner)
    assert options['model_hessian'] == expected
    assert options['hessian_update'] == 'auto'


@pytest.mark.parametrize('name', ['H2O_RHF-DFT_OPTIMIZE_OQP',
    'HCN_RHF-DFT_TS_OQP', 'HCN_BHHLYP-MRSFTDDFT_TS_OQP'])
def test_trajectory_references_pin_old_model_while_auto_examples_keep_new_default(name):
    folder = ROOT / 'examples/OPT'
    path = folder / (name + '.oqp')
    original = INPUT.resolve_oqp_text(path.read_text(), source_path=path).legacy_config
    auto_path = folder / (name + '_AUTO.oqp')
    automatic = INPUT.resolve_oqp_text(auto_path.read_text(), source_path=auto_path).legacy_config
    assert (folder / (name + '.json')).is_file()
    assert original['oqp'].pop('model_hessian') == 'constant'
    assert 'model_hessian' not in automatic['oqp']
    assert original == automatic  # Same method/state/geometry and iteration controls.
    assert report_for(automatic['oqp'], runtype=automatic['input']['runtype']).ok
