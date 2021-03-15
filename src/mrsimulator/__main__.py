# -*- coding: utf-8 -*-
import getopt
import sys

from .benchmark import Benchmark

if __name__ == "__main__":
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, "hb:", ["benchmark="])
    except getopt.GetoptError:
        print("--benchmark=<option>")
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-h":  # help
            print("--benchmark=<option>")
        elif opt in ("-b", "--benchmark"):  # benchmark
            getattr(Benchmark, "prep")()
            getattr(Benchmark, arg)()
