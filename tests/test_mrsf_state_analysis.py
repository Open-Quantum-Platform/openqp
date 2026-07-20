"""Focused tests for physical-root MRSF state labels."""
import importlib.util
import os
import sys
import types

import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
_ANALYSIS = os.path.join(_ROOT, "pyoqp", "oqp", "analysis")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        f"_state_analysis_test.{name}", os.path.join(_ANALYSIS, f"{name}.py"))
    module = importlib.util.module_from_spec(spec)
    module.__package__ = "_state_analysis_test"
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


package = types.ModuleType("_state_analysis_test")
package.__path__ = [_ANALYSIS]
sys.modules[package.__name__] = package
_load("descriptors")
_load("nto")
state_analysis = _load("state_analysis")


class _Mol:
    @staticmethod
    def get_atoms():
        return np.array([8, 6])


class _States:
    nstates = 2
    nbf = 2
    S = np.eye(2)
    C = np.eye(2)
    energies = np.array([-10.0, -9.8])
    mol = _Mol()

    @staticmethod
    def tdm_mo(i, j):
        assert (i, j) == (0, 1)
        return np.array([[0.0, 0.0], [1.0, 0.0]])

    @classmethod
    def tdm_ao(cls, i, j):
        return cls.tdm_mo(i, j)


class _AO:
    nbf = 2
    ao_atom = np.array([0, 1])
    ao_index = [(0, (1, 0, 0), 1.0), (1, (0, 0, 1), 1.0)]
    coords = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])


def test_n_pi_star_and_ct_are_reported_separately():
    result = state_analysis.analyze_mrsf_transition(
        _States(), _AO(), 1, [[0], [1]], plane_normal=[0, 0, 1])
    assert result["pair"] == (0, 1)
    assert result["ct_fraction"] == 1.0
    assert result["le_fraction"] == 0.0
    assert result["orbital_character"]["label"] == "nπ*"
    fractions = result["orbital_character"]["fractions"]
    assert fractions["n->pi*"] == 1.0
    assert abs(sum(fractions.values()) - 1.0) < 1.0e-12
    assert result["nto_participation_ratio"] == 1.0


def test_mixed_label_retains_channel_fractions():
    result = state_analysis.analyze_mrsf_transition(
        _States(), _AO(), 1, [[0], [1]], plane_normal=[1, 0, 1],
        label_threshold=0.75)
    assert result["orbital_character"]["label"] == "mixed"
    assert abs(sum(result["orbital_character"]["fractions"].values()) - 1.0) < 1.0e-12


def test_le_ct_survives_when_no_unique_plane_exists():
    result = state_analysis.analyze_mrsf_transition(
        _States(), _AO(), 1, [[0], [1]])
    assert result["ct_fraction"] == 1.0
    assert result["orbital_character"]["label"] == "unclassified"
    assert not result["orbital_character"]["classified"]
    assert "three atoms" in result["orbital_character"]["reason"]
