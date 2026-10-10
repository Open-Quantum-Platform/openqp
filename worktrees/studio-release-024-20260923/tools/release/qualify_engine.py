"""Run the bundled engine's NMR/ACID example in an isolated clean runtime."""
import argparse
import json
import math
import os
import subprocess
from pathlib import Path


def qualify(directory: Path, output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=False)
    executable = directory / ("openqp.exe" if os.name == "nt" else "openqp")
    source = directory / "examples/NMR/H2O_RHF-GIAO-ACID.inp"
    (output / "acid.inp").write_text(source.read_text().replace(
        "[guess]", "[guess]\nsave_mol=True"))
    env = {k: v for k, v in os.environ.items() if k in (
        "HOME", "USERPROFILE", "SYSTEMROOT", "WINDIR", "TEMP", "TMP", "COMSPEC")}
    env.update(PATH=(r"C:\Windows\System32;C:\Windows" if os.name == "nt" else "/usr/bin:/bin"),
               OMP_NUM_THREADS="2", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (output / "execution.log").open("w") as log:
        subprocess.run([str(executable.resolve()), "acid.inp", "--nompi"], cwd=output,
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=600)
    data = json.loads((output / "acid.json").read_text())
    assert abs(data["energy"] - (-74.9609154815971)) < 1e-7, data["energy"]
    assert len(data["nmr_shielding"]) == 3
    cubes = []
    for suffix in ("acid", "jx", "jy", "jz"):
        path = output / f"acid_{suffix}.cube"
        lines = path.read_text().splitlines()
        natoms = abs(int(lines[2].split()[0]))
        shape = [abs(int(lines[i].split()[0])) for i in (3, 4, 5)]
        values = [float(x) for line in lines[6 + natoms:] for x in line.split()]
        assert len(values) == math.prod(shape)
        assert all(math.isfinite(value) for value in values)
        cubes.append(path.name)
    evidence = {"engine": json.loads((directory / "engine-identity.json").read_text()),
                "energy": data["energy"], "nmr_atoms": 3, "cubes": cubes,
                "clean_runtime": True, "passed": True}
    (output / "qualification.json").write_text(json.dumps(evidence, indent=2) + "\n")
    return evidence


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(qualify(args.engine, args.output), indent=2))
