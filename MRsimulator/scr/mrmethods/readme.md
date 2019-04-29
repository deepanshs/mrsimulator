# Fortran file compilation

The file `powder.f90` contains the code for the powder interpolation scheme by
Alderman, Solum and Grant, J. Chem. Phys, 84, 1985.
DOI: [10.1063/1.450211](https://aip.scitation.org/doi/10.1063/1.450211)

To compile ths fortran code, execute the following in the command line. A fortran compiler is required for compilation.

`f2py -m powder -c powder.f90`
