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


def test_cli_n_jobs_without_value():
    """--n_jobs without a value is one job."""
    from mrsimulator.__main__ import parser

    assert parser.parse_args(["--n_jobs"]).n_jobs == 1
    assert parser.parse_args(["--n_jobs", "3"]).n_jobs == 3
    assert parser.parse_args([]).n_jobs == 1


def test_cli_fewer_processors_than_methods(tmp_path, monkeypatch):
    """Every method is saved, also when the file holds fewer signal processors than
    methods. The methods without a processor are saved unprocessed."""
    import argparse

    from mrsimulator import signal_processor as sp
    from mrsimulator.__main__ import run_full_simulator

    methods = [
        BlochDecaySpectrum(channels=["1H"], spectral_dimensions=[{"count": 32}])
        for _ in range(2)
    ]
    site = Site(isotope="1H", isotropic_chemical_shift=1)
    sim = Simulator(spin_systems=[SpinSystem(sites=[site])], methods=methods)
    processors = [sp.SignalProcessor(operations=[sp.Scale(factor=2)])]
    save(str(tmp_path / "two.mrsim"), sim, processors)

    monkeypatch.chdir(tmp_path)
    args = argparse.Namespace(input="two.mrsim", output=None, n_jobs=1, plot=False)
    run_full_simulator(args)
    assert os.path.isfile(tmp_path / "two.mrsim_0.csdf")
    assert os.path.isfile(tmp_path / "two.mrsim_1.csdf")


def test_benchmark_output_without_a_terminal(capsys):
    """The benchmark report does not need a terminal."""
    from mrsimulator import benchmark

    benchmark.terminal_start_setup()
    benchmark.terminal_end_setup(1.0, 10, "description")
    assert "Computation method" in capsys.readouterr().out
