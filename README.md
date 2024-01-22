This repository contains the modifications I did to the double-precision MPI implementation of OpenDDA, namely:
- Implemented a [Corrected Lattice Dispersion Relation](https://doi.org/10.48550/arXiv.astro-ph/0403082) in `dipole_polarisabilities_MPI.c`
- Added a new function `scattered_field_MPI.c` that calculates the scattered field in terms of 4-by-4 complex-valued matrix, according to the conventions outlined in [this article](https://doi.org/10.1086/166795)
- Miscellaneous updates to deprecated MPI functions, and modifications to the standard output in the main `opendda_MPI.c` file.

To use this code, you need to install the two main library dependencies: 
- [FFTW](https://fftw.org/) 
- [Open MPI](https://www.open-mpi.org/).

To learn more about Open DDA, see the following file: [[About OpenDDA]]
