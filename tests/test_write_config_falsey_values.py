from oqp.utils.file_utils import write_config


def test_write_config_preserves_false_zero_and_unpruned_grid():
    text, lowered = write_config({
        "dftgrid": {
            "rad_npts": 96,
            "ang_npts": 302,
            "pruned": "",
            "grid_ao_pruned": False,
            "grid_ao_threshold": 0.0,
        },
        "tests": {"exception": False},
    })

    assert "pruned=none\n" in text
    assert "grid_ao_pruned=False\n" in text
    assert "grid_ao_threshold=0.0\n" in text
    assert "exception=False\n" in text
    assert lowered["dftgrid"]["pruned"] == "none"
    assert lowered["dftgrid"]["grid_ao_pruned"] == "False"
    assert lowered["dftgrid"]["grid_ao_threshold"] == "0.0"
    assert lowered["tests"]["exception"] == "False"


def test_write_config_still_omits_empty_nonsemantic_strings_and_lists():
    text, lowered = write_config({
        "input": {"system2": ""},
        "properties": {"grad": []},
    })

    assert "system2=" not in text
    assert "grad=" not in text
    assert lowered == {"input": {}, "properties": {}}
