
.. _simplify_spin_systems:

=====================
Simplify Spin Systems
=====================

In ``mrsimulator``, it is possible to define a spin system and then simplify it to
irreducible spin systems using :ref:`.simplify` method. Let's start by importing and
intializing the necessary objects.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator import Site, Coupling, SpinSystem
    >>> import pprint
    >>> pp = pprint.PrettyPrinter()


Now, let's define an uncoupled spin system with six sites.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> A = Site(isotope="1H", isotropic_chemical_shift=0, name="a")
    >>> B = Site(isotope="1H", isotropic_chemical_shift=2, name="b")
    >>> C = Site(isotope="1H", isotropic_chemical_shift=4, name="c")
    >>> D = Site(isotope="1H", isotropic_chemical_shift=6, name="d")
    >>> E = Site(isotope="1H", isotropic_chemical_shift=8, name="e")
    >>> F = Site(isotope="1H", isotropic_chemical_shift=10, name="f")
    >>> sites = [A, B, C, D, E, F]
    >>> uncoupled_sys = SpinSystem(sites=sites, abundance=10)


When we run the :ref:`simplify` method, we expect to see a list of six spin systems
containing one site each.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> simplified_uncoupled = uncoupled_sys.simplify()
    >>> simplified_uncoupled = [i.json() for i in simple_sys]
    >>> pp.pprint(simple_sys)
    [{'abundance': '10.0 %',
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '0.0 ppm',
                 'name': 'a'}]},
     {'abundance': '10.0 %',
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '2.0 ppm',
                 'name': 'b'}]},
     {'abundance': '10.0 %',
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '4.0 ppm',
                 'name': 'c'}]},
     {'abundance': '10.0 %',
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '6.0 ppm',
                 'name': 'd'}]},
     {'abundance': '10.0 %',
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '8.0 ppm',
                 'name': 'e'}]},
     {'abundance': '10.0 %',
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '10.0 ppm',
                 'name': 'f'}]}]


Now, we have a list of six systems, just as we expected.  Let's define another reducible
spin system that has some couplings between sites.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> AB_couple = Coupling(site_index=[0, 1], isotropic_j=10, name="AB")
    >>> BC_couple = Coupling(site_index=[1, 2], isotropic_j=10, name="BC")
    >>> DF_couple = Coupling(site_index=[3, 5], isotropic_j=30, name="DF")
    >>> couplings = [AB_couple, BC_couple, DF_couple]
    >>> coupled_sys_1 = SpinSystem(sites=sites, couplings=couplings, abundance=30)

We expect the simplified system to be a list of three spin systems, one containing sites
A, B, and C with their couplings, another containing sites D and F and their coupling,
and another containing only site E.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> simplified_coupled_1 = coupled_sys_1.simplify()
    >>> simplified_coupled_1 = [i.json() for i in simple_sys]
    >>> pp.pprint(simple_sys)
    [{'abundance': '30.0 %',
      'couplings': [{'isotropic_j': '10.0 Hz', 'name': 'AB', 'site_index': [0, 1]},
                    {'isotropic_j': '10.0 Hz', 'name': 'BC', 'site_index': [1, 2]}],
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '0.0 ppm',
                 'name': 'a'},
                {'isotope': '1H',
                 'isotropic_chemical_shift': '2.0 ppm',
                 'name': 'b'},
                {'isotope': '1H',
                 'isotropic_chemical_shift': '4.0 ppm',
                 'name': 'c'}]},
     {'abundance': '30.0 %',
      'couplings': [{'isotropic_j': '30.0 Hz', 'name': 'DF', 'site_index': [0, 1]}],
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '6.0 ppm',
                 'name': 'd'},
                {'isotope': '1H',
                 'isotropic_chemical_shift': '10.0 ppm',
                 'name': 'f'}]},
     {'abundance': '30.0 %',
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '8.0 ppm',
                 'name': 'e'}]}]
