# -*- coding: utf-8 -*-
import getopt
import sys

from .benchmark import Benchmark

if __name__ == "__main__":
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hb:n:", ["benchmark=", "n_jobs="])
    except getopt.GetoptError:
        print("--benchmark=<option>")
        sys.exit(2)

    benchmark_level = None
    n_jobs = 1
    for opt, arg in opts:
        if opt == "-h":  # help
            print("--benchmark=<option>")
            break
        if opt in ("-b", "--benchmark"):  # benchmark
            benchmark_level = arg
        if opt in ("-n", "--n_jobs"):  # n_jobs
            n_jobs = arg

    if benchmark_level is not None:
        getattr(Benchmark, "prep")()
        getattr(Benchmark, benchmark_level)(n_jobs=int(n_jobs))
