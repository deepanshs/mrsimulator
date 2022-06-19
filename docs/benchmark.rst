.. _benchmark:

=====================
Performance benchmark
=====================

One of the objectives in the design of the ``mrsimulator`` library is to enable
fast NMR spectrum simulation.
For this, we have put considerable effort into optimizing the library.
The following benchmark shows the library's performance in computing the
solid-state NMR spectra from single-site spin systems for the shift and
quadrupolar tensor interactions at static and MAS conditions.



.. A benchmark for the number of single-site spin systems computer per second.

.. figure:: _static/benchmark.*
    :alt: benchmark

    (Left) The number of single-site spin systems computer per seconds. (Right)
    The execution time (in ms) in computing spectrum from a single-site spin system.


.. A similar benchmark showing the execution time of a single-site spin system. Lower
.. is better.

.. .. figure:: _static/benchmark_time.*
..     :figclass: figure

..     The execution time (in ms) in computing spectrum from a single-site spin system.

**Benchmark specs**

The benchmarks were performed on a 2.3 GHz Quad-Core Intel Core i5 Laptop using 8
GB 2133 MHz LPDDR3 memory. For consistent benchmarking, 1000 single-site
spin systems were constructed, where the tensor parameters of the sites (*zeta*
and *eta* for the shielding tensor, and *Cq* and *eta* for the quadrupolar
tensor) were randomly populated. The execution time for this setup was recorded,
and the process was repeated 70 times. The reported value is the mean and the
standard deviation.

All calculations were performed using the default Simulator
:attr:`~mrsimulator.Simulator.config` attribute values.

Benchmark for the previous versions
-----------------------------------

.. figure:: _static/benchmark_previous.*
    :alt: previous benchmark

    (Left) The number of single-site spin systems computer per seconds. (Right)
    The execution time (in ms) in computing spectrum from a single-site spin system.
