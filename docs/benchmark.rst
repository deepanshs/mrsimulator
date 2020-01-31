.. _benchmark:

=========
Benchmark
=========

One of the objectives in the design of the ``Mrsimulator`` library is to enable
fast NMR lineshape simulation for amorphous materials. For this, we have
put a considerable effort into optimizing the library. The following benchmark
shows the number of sites simulated per second for the shift and quadrupolar
tensor interaction at static and MAS conditions.

.. figure:: _static/benchmark.*
    :figclass: figure-polaroid

The benchmark was performed on a 2.3 GHz Quad-Core Intel Core i5 Laptop using 8
GB 2133 MHz LPDDR3 memory. For consistent benchmarking, 1000 single-site
isotopomers were constructed, where the tensor parameters of the sites (`zeta`
and `eta` for the shielding tensor, and `Cq` and `eta` for the quadrupolar
tensor) were randomly populated. The execution time for this setup was
recorded, and the process repeated 70 times, giving a mean and a standard
deviation. The number of sites computed per second is then given by

``number_of_sites_computed_per_second = 1000/mean_time``.

All calculations were performed using the default Simulator
:attr:`~mrsimulator.simulator.Simulator.config` attribute values.
