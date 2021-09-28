.. _custom_method:

Writing custom methods
======================

Out of the box, ``mrsimulator`` comes with a variety of build-in methods for simulating NMR spectra such as the Bloch decay spectrum, MQMAS, STMAS, spinning sidebands correlation spectrum, and many more. These build-in methods provide a simple and quick simulation set up for improved user workflow. 

The core library, however, is designed such that the end-users can write their own simulation methods. At present, we provide two generic method objects, *Method1D* and *Method2D*, for writing one and two-dimensional simulation methods, respectively. 

Every method object has attributes: channels, magnetic_flux_density, rotor_angle, rotor_frequency, and spectral_dimensions. The value of the channels attribute is a list of isotopes that are involved in the method. The value of magnetic_flux_density, rotor_angle, and rotor_frequency are the static external magnetic flux density, sample rotation angle, and sample rotation frequency, respectively, that are defined globally. We will get to what globally means in a moment. The value of the spectral_dimensions attribute is a list of **SpectralDimension** objects. Note, by definition, *Method1D* and *Method2D* hold a list of one and two SpectralDimension objects respectively. Each SpectralDimension object has attributes: count, spectral_width, and reference_offset that collectively defines a coordinate grid over which the frequencies are sampled.

The attributes described above are the same as used with the build-in methods, and should be familiar to you. Besides this, a SpectralDimension object also has an *events* attribute, whose value is a list of Event objects. For build-in methods, we have pre-defined the events, and these cannot be modified by the users. In the following sec

In mrsimulator, event objects play a key role in how frequencies are evaluated along the spectral dimension. There are two types of Event objects: SpectralEvent and MixingEvent. A SpectralEvent has attributes: magnetic_flux_density, rotor_angle, rotor_frequency, and transition_query. The magnetic_flux_density, rotor_angle, and rotor_frequency are defined similar to the SpectralDimension, except these values are local to the event. If the attribute is not defined with in the event, the global value of the respective attribute is used instead. is what contributes to the frequencies. The MixingEvent defines how the transitions are connection in a multi-event method.

In a build-in method, 


You may have seen these with the build-in methods as well. Aside from these grid generating p

In Method1D and Method2D, the length of these objects is one and two, respectively.

core objects that you will frequently use when writing a method are the *SpectralDimension* and *Event* objects, where the attribute *spectra_dimensions* and *events* attribute can have a *SpectralDimension* objects are 


using the core objects. All build-in mrsimulator methods use these core objects.

When writing the API for building mrsimulator methods, our focuses attention is versatile to allow advanced users to write their own custom methods. In this section, we describe how one can 

The ``mrsimulator`` library allows the users to write their own custom methods. Before describing 
