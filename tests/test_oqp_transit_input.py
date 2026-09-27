"""Endpoint-option lowering, validation and canonical round trips."""
import importlib.util
from pathlib import Path
import sys
import pytest
ROOT = Path(__file__).resolve().parents[1]
def load(name,path):
    spec = importlib.util.spec_from_file_location(name,ROOT/path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module
INPUT=load('_transit_input','pyoqp/oqp/utils/oqp_input.py')
CHECK=load('_transit_checker','pyoqp/oqp/utils/input_checker.py')

@pytest.mark.parametrize('name',['NH3_RHF-HF_QST2_OQP','NH3_RHF-HF_QST3_OQP','NH3_RHF-HF_IDPP_NEB_OQP'])
def test_examples_lower_and_roundtrip(name):
    folder=ROOT/'examples/OPT'
    text=(folder/(name+'.oqp')).read_text()
    result=INPUT.resolve_oqp_text(text,source_path=folder/(name+".oqp"))
    config=result.legacy_config
    assert config['optimize'].get('lib','oqp')=='oqp'
    roundtrip = INPUT.resolve_oqp_text(result.canonical_text,source_path=folder/(name+".oqp")).legacy_config
    # Canonical rendering intentionally omits explicit schema defaults.
    for c in (config, roundtrip):
        c.setdefault("oqp", {}).setdefault("init_hessian", "model")
        if "neb" in c:
            c["neb"].setdefault("nimage", "5")
            c["oqp"].setdefault("align", "True")
    assert roundtrip == config
    if 'QST' in name:
        assert Path(config['oqp']['ts_product']).is_file()
        assert config['oqp']['init_hessian']=='model'
    else:
        assert config['oqp']['neb_interpolation']=='idpp'

@pytest.mark.parametrize('options,fragment',[
    ({'ts_search':'qst2'},'ts_product'),
    ({'ts_search':'qst3','ts_product':'p.xyz'},'ts_guess'),
    ({'ts_search':'qst2','ts_product':'p.xyz','init_hessian':'numerical'},'init_hessian'),
    ({'ts_search':'qst2','ts_product':'p.xyz','freeze':'distance(1,2)'},'freeze'),
    ({'ts_search':'wrong'},'ts_search'),
    ({'ts_guess':'g.xyz'},'ts_search'),
])
def test_invalid_search_settings_rejected(options,fragment):
    report=CHECK.CheckReport()
    CHECK._check_optimize({'input':{'runtype':'ts','method':'hf'},'optimize':{'lib':'oqp'},'oqp':options},report)
    assert not report.ok
    assert fragment in report.to_text()
