.. _single_site_system_generator_documentation:

============================
Single Site System Generator
============================

Custom Site and SpinSystem parameters can be passed to the
:py:meth:`~mrsimulator.utils.collection.single_site_system_generator` method to create a list
of uncoupled spin systems. Each :ref:`spin_sys_api` in the returned list holds a single
:ref:`site_api` object since the backend simulation is more efficient for single site spin
systems. Import the method as below

.. code-block:: python

    from mrsimulator.utils.collection import single_site_system_generator

The arguments passed to the function, defined in :numref:`single_site_sys_gen_table`,
can either be a scalar quantity (``float`` or ``str``, where applicable) or a
``list``/``np.array`` of those quantities. All lists passed must
have the same length, otherwise an error will be thrown. For example,

.. skip: start

.. code-block:: python

    single_site_system_generator(
        isotope=["1H", "1H", "13C", "17O"],
        isotropic_chemical_shift=[1.3, 3.7, 65.0],
    )

.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    Traceback (most recent call last):
    ...
    ValueError: An array or list was either too short or too long. All arguments must be the
    same size. If one attribute is a type list of length n, then all attributes with list types
    must also be of length n, and all remaining attributes must be scalar (singular float, int,
    or str).

.. skip: end

The attributes of each returned spin system at a certain index correspond to the attribute passed
at that index. For example,

.. code-block:: python

    single_site_system_generator(
        isotope=["1H", "1H", "13C"],
        isotropic_chemical_shift=[1.3, 3.7, 65.0],
    )

returns a list of 3 spin systems. The first two spin systems represent proton sites with isotropic
chemical shifts of 1.3 ppm and 3.7 ppm, respectively. The third spin system is a
:math:`^{13}\text{C}` site with a chemical shift of 65.0 ppm.

Broadcasting Length of List
---------------------------

Arguments passed as a single value will be broadcast to a list of that value with the same
length as other lists passed. For example

.. code-block:: python

    single_site_system_generator(
        isotope=["1H", "1H", "1H"],
        isotropic_chemical_shift=2.0,
    )

is equivalent to calling

.. code-block:: python

    single_site_system_generator(
        isotope=["1H", "1H", "1H"],
        isotropic_chemical_shift=[2.0, 2.0, 2.0],
    )

Passing lists of Tensor Parameters
----------------------------------

Tensor parameters for sites are passed as dictionaries where the keywords represent the tensor
attribute and the values are single values or a ``list``/``np.array`` of values. Again, these
lists must have the same length of all other lists passed. Single values will be broadcast to a
list of that value with the same length as other lists passed. For example

.. code-block:: python

    single_site_system_generator(
        isotope="13C",
        shielding_symmetric={
            "zeta": [5, 10, 15, 20, 25],
            "eta": 0.3,
        },
    )

returns a list of five :math:`^{13}\text{C}` spin systems with different ``shielding_symmetric.zeta``
values but the same ``shielding_symmetric.eta`` value.

If you need to intermix sites with and without tensor parameters, simply put ``None`` at the index
of the site without the tensor parameter.

.. code-block:: python

    single_site_system_generator(
        isotope=["1H", "17O"],
        quadrupolar={
            "Cq": [None, 3.2e6],
            "eta": [None, 0.5],
        },
    )

.. minigallery:: mrsimulator.utils.collection.single_site_system_generator
  :add-heading: Examples using ``single_site_system_generator()``
  :heading-level: -

--------------------------------------------------------------------------------

.. cssclass:: table-bordered table-striped centered
.. _single_site_sys_gen_table:
.. list-table:: Arguments for ``single_site_system_generator``
    :widths: 30 20 50
    :header-rows: 1

    * - Name
      - Type
      - Description

    * - ``isotope``
      - ``str`` or list of ``str``
      - A *required* string or list of strings representing the label of the ``isotope`` attribute
        of the :ref:`site_api` (e.g. ``"1H"`` or ``["29Si", "17O"]``).

    * - ``isotropic_chemical_shift``
      - ``float``, list of ``float``, or numpy array
      - An *optional* number or list of numbers representing the ``isotropic_chemical_shift``
        attribute of the :ref:`site_api` (e.g. ``17.3`` or ``[2.4, 19.5]``) in ppm.
        The default value is ``0``.

    * - ``shielding_symmetric``
      - ``dict``
      - An *optional* dictionary representing the ``shielding_symmetric`` attribute of the :ref:`site_api`
        where the keys are valid :ref:`sy_api` attributes and the values are floats or lists/numpy
        arrays of floats. The default is ``None``.

    * - ``shielding_antisymmetric``
      - ``dict``
      - An *optional* dictionary representing the ``shielding_antisymmetric`` attribute of the
        :ref:`site_api` where the keys are valid :ref:`asy_api` attributes and the values are floats
        or lists/numpy arrays of floats. The default is ``None``.

    * - ``quadrupolar``
      - ``dict``
      - An *optional* dictionary representing the ``quadrupolar`` attribute of the
        :ref:`site_api` where the keys are valid :ref:`sy_api` attributes and the values are floats
        or lists/numpy arrays of floats. The default is ``None``.

    * - ``abundance``
      - ``float``, list of ``float``, or numpy array
      - An *optional* number or list of numbers representing the ``abundance`` attribute of
        the SpinSystem (e.g. ``0.182`` or ``[85, 7.3]``. By default, the abundance
        of each spin system will be set to ``1 / n_sys`` where ``n_sys`` is the number of spin
        systems generated.

    * - ``site_name``
      - ``str`` or list of ``str``
      - An *optional* string or list of strings representing the ``name`` attribute of each
        :ref:`site_api`. By default, each :ref:`site_api` will take the default name of ``None``

    * - ``site_label``
      - ``str`` or list of ``str``
      - An *optional* string or list of strings representing the ``label`` attribute of each
        :ref:`site_api`. By default, each :ref:`site_api` will take the default label of ``None``

    * - ``site_description``
      - ``str`` or list of ``str``
      - An *optional* string or list of strings representing the ``description`` attribute of each
        :ref:`site_api`. By default, each :ref:`site_api` will take the default description of ``None``
