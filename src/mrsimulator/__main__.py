import argparse

import csdmpy as cp
import matplotlib.pyplot as plt

from . import __version__
from . import load
from .benchmark import Benchmark


def run_benchmark(args):
    """Run benchmark tests"""
    getattr(Benchmark, "prep")()
    getattr(Benchmark, args.benchmark)(
        n_jobs=int(args.n_jobs),
        interpolation=args.interpolation,
        simulation=args.simulation,
    )


def run_full_simulator(args):
    """Read, run, and save output of .mrsim file"""
    sim, signal_processor, _ = load(args.input)
    sim.run(n_jobs=args.n_jobs)
    proc_data = []
    for proc, mth in zip(signal_processor, sim.methods):
        proc_data.append(proc.apply_operations(mth.simulation))

    outfile = args.output or args.input
    segments = outfile.split(".")
    outname, ext = "".join(segments[:-1]), segments[-1]
    if ext not in ["csdf"]:
        outname = f"{outname}.{ext}"
        ext = "csdf"

    for i, f in enumerate(proc_data):
        f.save(f"{outname}_{i}.{ext}")

        if args.plot:
            cp.plot(f.real)
            plt.show()


def run():
    if args.input is not None:
        run_full_simulator(args)

    if args.version:
        print(f"mrsimulator {__version__}")

    if args.benchmark is not None:
        run_benchmark(args)


parser = argparse.ArgumentParser(description="Mrsimulator CLI")
parser.add_argument(
    "input", metavar="i", type=str, nargs="?", help="path to the .mrsim input file."
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    nargs="?",
    help="CSDM output file name. Extension is '.csdf'.",
)
parser.add_argument("--plot", action="store_true", help="plot the csdf output.")
parser.add_argument(
    "--n_jobs",
    type=int,
    nargs="?",
    default=1,
    help="number of processors for parallel computation.",
)
parser.add_argument("-v", "--version", action="store_true", help="mrsimulator version.")
parser.add_argument(
    "--benchmark", choices=["l0", "l1", "l2"], help="set benchmark level."
)
parser.add_argument(
    "--interpolation",
    action="store_true",
    help="run interpolation benchmark. Default is False.",
)
parser.add_argument(
    "--simulation",
    action="store_true",
    help="run simulation benchmark. Default is False.",
)
args = parser.parse_args()


if __name__ == "__main__":
    run()
