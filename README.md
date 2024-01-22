This repository contains the modifications I did to the double-precision MPI implementation of OpenDDA, namely:
- Implemented a [Corrected Lattice Dispersion Relation](https://doi.org/10.48550/arXiv.astro-ph/0403082) in `dipole_polarisabilities_MPI.c`
- Added a new function `scattered_field_MPI.c` that calculates the scattered field in terms of 4-by-4 complex-valued matrix, according to the conventions outlined in [this article](https://doi.org/10.1086/166795)
- Miscellaneous updates to deprecated MPI functions, and modifications to the standard output in the main `opendda_MPI.c` file.

To compile this code:
1. Install the two main library dependencies: 
- [FFTW](https://fftw.org/) 
- [Open MPI](https://www.open-mpi.org/)
2. Update the Makefile with path to `include` and `lib` folders. The specific location will be installation-dependent.
3. Update Makefile with compiler (`CC`) and relevant flags (`CFLAGS`)
4. Type `make` in the terminal to compile the program. The executable is `opendda_MPI`

To learn more about Open DDA, see the [About OpenDDA](https://github.com/bcrodrigo/OpenDDA_MPI_modifications/blob/main/About%20OpenDDA.md) file.
