.. _benchmark:

=========
Benchmark
=========

The primary objective in the design of ``Mrsimulator`` is to enable fast
simulation of NMR lineshapes arising for amorphous materials. For this, we have
put in a lot of effort into optimizing the library. The following benchmark
shows the number of sites simulated per second for the shift and quadrupolar
tensor interaction at static and MAS conditions.

.. figure:: _images/benchmark.*
    :figclass: figure-polaroid

.. .. raw:: html

..     <iframe src="_static/benchmark_result.html" height="475px" width="100%" frameBorder="0"></iframe>

The benchmark was performed on a 2.3 GHz Quad-Core Intel Core i5 Laptop using
8 GB 2133 MHz LPDDR3 memory. For consistent benchmarking, we created 1000
single-site isotopomers, where the tensor parameters of the sites (zeta and eta
for shielding tensor, and Cq and eta for quadrupolar tensor) were randomly
populated. The execution time of simulating 1000 isotopomers was then recorded.
This process was repeated 70 times, and the execution time from the runs was
used to calculate the mean and the standard deviation. To obtain the number of
sites computed per second, we performed the following calculation,
``number_of_sites_computed_per_second = 1000/mean_time``. All calculations were
performed using the default Simulator
:attr:`~mrsimulator.simulator.Simulator.config` attribute values.
