# -*- coding: utf-8 -*-
import getopt
import sys

from .benchmark import Benchmark


class Main:
    def __init__(self):
        self.benchmark_level = None
        self.n_jobs = 1

    def get_args(self, opts):
        for opt, arg in opts:
            if opt == "-h":  # help
                print("--benchmark=<option=l0,l1,l2>")
                break
            if opt == "--benchmark":  # benchmark
                if int(arg[-1]) > 2:
                    allow = [f"l{i}" for i in range(3)]
                    print(f"Allowed levels are {', '.join(allow)}")
                    sys.exit(2)
                self.benchmark_level = arg
            if opt == "--n_jobs":  # n_jobs
                self.n_jobs = arg

    def benchmark(self):
        if self.benchmark_level is not None:
            getattr(Benchmark, "prep")()
            getattr(Benchmark, self.benchmark_level)(n_jobs=int(self.n_jobs))


if __name__ == "__main__":
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "h", ["benchmark=", "n_jobs="])
    except getopt.GetoptError:
        print("--benchmark=<option>")
        sys.exit(2)

    start = Main()
    start.get_args(opts)
    start.benchmark()
