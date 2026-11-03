"""Tests for the command line interface, `python -m mrsimulator`."""
import os
import subprocess
import sys

from mrsimulator import save
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecaySpectrum


def _write_mrsim(path):
    method = BlochDecaySpectrum(channels=["1H"], spectral_dimensions=[{"count": 32}])
    site = Site(isotope="1H", isotropic_chemical_shift=1)
    sim = Simulator(spin_systems=[SpinSystem(sites=[site])], methods=[method])
    save(path, sim)


def _run_cli(cwd, *args):
    return subprocess.run(
        [sys.executable, "-m", "mrsimulator", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def test_cli_output_file_names(tmp_path):
    """The output name keeps the directories and dots of the input path."""
    os.makedirs(tmp_path / "sub.dir")
    _write_mrsim(str(tmp_path / "sub.dir" / "sample.mrsim"))
    result = _run_cli(tmp_path, os.path.join("sub.dir", "sample.mrsim"))
    assert result.returncode == 0, result.stderr
    assert os.path.isfile(tmp_path / "sub.dir" / "sample.mrsim_0.csdf")

    _write_mrsim(str(tmp_path / "sample.mrsim"))
    result = _run_cli(tmp_path, os.path.join(".", "sample.mrsim"))
    assert result.returncode == 0, result.stderr
    assert os.path.isfile(tmp_path / "sample.mrsim_0.csdf")

    result = _run_cli(tmp_path, "sample.mrsim", "-o", os.path.join(".", "out.v1.csdf"))
    assert result.returncode == 0, result.stderr
    assert os.path.isfile(tmp_path / "out.v1_0.csdf")


def test_cli_import_does_not_parse_arguments():
    """Importing the CLI module does not parse the command line arguments."""
    code = "import sys; sys.argv = ['x', '--unknown']; import mrsimulator.__main__"
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
