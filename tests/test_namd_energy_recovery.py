"""Exercise the production recovery block without a native SCF calculation."""
import ast
import copy
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock
import numpy as np
import pytest

SOURCE = Path(__file__).resolve().parents[1] / 'pyoqp/oqp/library/namd.py'


def make_recovery(residuals, converged=True, previous=True):
    source = SOURCE.read_text()
    tree = ast.parse(source)
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'NAMD')
    names = {'_energy_refinement_counts', '_energy_retry_state',
             '_restore_energy_retry_state',
             '_scale_analytic_velocity_contractions'}
    cls.body = [n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name in names]
    ns = dict(np=np, copy=copy, FS_TO_AU=1.0, dump_log=Mock())
    exec(compile(ast.Module(body=[cls], type_ignores=[]), str(SOURCE), 'exec'), ns)
    start = source.index('            refinement_attempted = False')
    end = source.index('            energy_before_transition = (', start)
    body = textwrap.indent(textwrap.dedent(source[start:end]), '    ')
    exec('def recover(self, mol, r, r_start, vel_start, accel_start, retry_state, istep=1):\n'
         + body + '\n    return r\n', ns)
    d = ns['NAMD']()
    d.disc_tol = .002
    d.disc_substeps = 10
    d.nve_gate = 'off'
    d._etot_prev = 0.0
    d._conservative_restraint_energy = 0.0
    d.prev_data = {'saved': [1]} if previous else None
    d.vel = np.ones((1, 1))
    d.mass = np.ones(1)
    d.active = 0
    d.dt = 1.0
    d.natom = 1
    d.rescale_provider = 'isotropic'
    d._somo_switch_step = False
    d._window_leak_step = False
    d._somo_switch_count = 0
    d.ref_switch_rescale = False
    d.disc_rescale = True
    d._disc_substep_events = d._disc_event_count = 0
    d._disc_energy_absorbed = 0.0
    d._ba_energy_center = ['accepted']
    d._evaluate_odp = lambda r: None
    d.mol = SimpleNamespace(energies=np.array([-.4]), put_data=Mock(), update_system=Mock())
    trace = []
    def advance(istep, r, vel, accel, dt, last, continuation):
        n = round(1/dt)
        trace.append(n)
        d.mol.energies = np.array([-.5 + residuals[n]])
        if last:
            d._ba_energy_center.append('trial')
        return np.ones_like(vel), np.zeros_like(accel), False
    d._advance_electronic_and_kick = advance
    d._state_overlap = Mock()
    def require():
        if not converged:
            raise RuntimeError('SCF did not converge')
    d._require_converged_reference = require
    def run():
        z = np.zeros((1, 1))
        return ns['recover'](d, d.mol, z, z, d.vel.copy(), z, d._energy_retry_state())
    return d, run, trace, ns['dump_log']


def test_refinement_stops_before_rescaling_when_energy_converges():
    d, run, trace, log = make_recovery({2: .02, 4: .001})
    run()
    assert trace == [2]*2 + [4]*4
    assert d._disc_event_count == 0
    np.testing.assert_array_equal(d.vel, [[1.]])
    assert d._ba_energy_center == ['accepted', 'trial']
    assert d.prev_data == {'saved': [1]}


def test_rescaling_only_after_all_configured_subdivisions():
    d, run, trace, log = make_recovery({2: .04, 4: .03, 8: .02, 10: .01})
    run()
    assert trace == [2]*2 + [4]*4 + [8]*8 + [10]*10
    assert d._disc_event_count == 1
    assert d._disc_energy_absorbed == pytest.approx(.01)
    assert float(d.mol.energies[0] + .5*np.sum(d.vel**2)) == pytest.approx(0., abs=1e-15)
    assert any('last-resort' in c.kwargs.get('title', '') for c in log.call_args_list)


def test_last_resort_correction_rescales_the_cached_analytic_tdc():
    d, run, trace, log = make_recovery({2: .04, 4: .03, 8: .02, 10: .01})
    tdc = np.array([[0.0, 0.4], [-0.4, 0.0]])
    d._last_analytic_tdc = tdc.copy()
    d._analytic_tdc_previous = tdc.copy()
    d._last_analytic_step = 1
    d._last_analytic_pair = None
    run()
    factor = float(d.vel[0, 0])
    assert d._disc_event_count == 1 and factor != pytest.approx(1.0)
    np.testing.assert_allclose(d._last_analytic_tdc, factor*tdc)
    np.testing.assert_allclose(d._analytic_tdc_previous, factor*tdc)


def test_no_rescaling_without_a_completed_refinement():
    d, run, trace, log = make_recovery({}, previous=False)
    run()
    assert not trace
    assert d._disc_event_count == 0
    np.testing.assert_array_equal(d.vel, [[1.]])


def test_unconverged_reference_cannot_be_rescaled():
    d, run, trace, log = make_recovery({2: .04, 4: .03, 8: .02, 10: .01}, converged=False)
    with pytest.raises(RuntimeError, match='SCF did not converge'):
        run()
    assert d._disc_event_count == 0
    np.testing.assert_array_equal(d.vel, [[1.]])


@pytest.mark.parametrize('active_changed, expected', [(False, np.nan), (True, .1)])
def test_numerical_correction_is_separate_from_hop_energy(active_changed, expected):
    source = SOURCE.read_text()
    start = source.index('            active_changed = new_active != active_old',
                         source.index('            refinement_attempted = False'))
    end = source.index('            self._apply_thermostat(istep)', start)
    d = SimpleNamespace(active=0, vel=np.ones((1, 1)), mass=np.ones(1),
                        _ref_switch_jump=.4, _conservative_restraint_force=0.,
                        _conservative_restraint_energy=0.,
                        _active_gradient=lambda: np.zeros((1, 1)))
    ns = dict(np=np, self=d, mol=SimpleNamespace(energies=np.array([-.5, -.4])),
              new_active=int(active_changed), active_old=0, bias_energy=0.,
              energy_before_transition=0.)
    exec(textwrap.dedent(source[start:end]), ns)
    if np.isnan(expected):
        assert np.isnan(ns['transition_energy_jump'])
    else:
        assert ns['transition_energy_jump'] == pytest.approx(expected)


def test_insufficient_kinetic_energy_does_not_count_as_a_correction():
    d, run, trace, log = make_recovery({2: .8, 4: .7, 8: .6, 10: .55})
    run()
    assert d._disc_event_count == 0
    assert d._disc_energy_absorbed == 0.
    np.testing.assert_array_equal(d.vel, [[1.]])
