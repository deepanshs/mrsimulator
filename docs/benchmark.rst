.. _benchmark:

=========
Benchmark
=========

One of the objectives in the design of the ``Mrsimulator`` library is to enable
fast NMR line-shape simulation.
For this, we have put a considerable effort into the optimization of the library.
The following benchmark shows the performance of the `Mrsimulator` library in computing
the solid-state NMR line-shapes of single-site spin systems for the shift and
quadrupolar tensor interaction at static and MAS conditions.

Computational specs
-------------------

A benchmark for the number of single-site spin systems computer per second.
Higher is better.

.. figure:: _static/benchmark_sites.*
    :figclass: figure

    The number of single-site spin systems computer per seconds.


A similar benchmark showing the execution time of a single-site spin system. Lower
is better.

.. figure:: _static/benchmark_time.*
    :figclass: figure

    The execution time (in ms) in computing a line-shape of a single-site spin system.

The benchmarks were performed on a 2.3 GHz Quad-Core Intel Core i5 Laptop using 8
GB 2133 MHz LPDDR3 memory. For consistent benchmarking, 1000 single-site
spin systems were constructed, where the tensor parameters of the sites (`zeta`
and `eta` for the shielding tensor, and `Cq` and `eta` for the quadrupolar
tensor) were randomly populated. The execution time for this setup was
recorded, and the process repeated 70 times, giving a mean and a standard
deviation.

All calculations were performed using the default Simulator
:attr:`~mrsimulator.Simulator.config` attribute values.
