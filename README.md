# Background

The Discrete-Dipole Approximation (DDA) is a numerical method that represents an arbitrarily-shaped particle as electric dipoles on a cubic lattice. This representation yields a system of linear equations for the induced dipole moments (light-matter interaction). If there are $N$ dipoles representing a particle, this results in a system of linear equations with $3N$ complex-valued unknowns. Typical values of $N$ range between 10,000 and 100,000 so the equations are solved with iterative techniques (i.e. Conjugate Gradient) in high-performance computing clusters. The solutions are used to calculate the theoretical scattering response for any observation direction.

There are many implementations, but the one I used was [OpenDDA](https://github.com/drjmcdonald/OpenDDA). In particular, I made some modifications to use the parallelized version using the Message-Passing Interface (MPI) for distributed memory clusters. See below for more details.

# Repository Description

This repository contains the modifications I did to the double-precision MPI implementation of OpenDDA, namely:
- Implemented a [Corrected Lattice Dispersion Relation](https://doi.org/10.48550/arXiv.astro-ph/0403082) in `dipole_polarisabilities_MPI.c`
- Added a new function `scattered_field_MPI.c` that calculates the scattered field in terms of 4-by-4 complex-valued matrix, according to the conventions outlined in [this article](https://doi.org/10.1086/166795)
- Miscellaneous updates to deprecated MPI functions, and modifications to the standard output in the main `opendda_MPI.c` file

## To Compile

1. Download the source code
2. Install the two main library dependencies: 
- [FFTW](https://fftw.org/) 
- [Open MPI](https://www.open-mpi.org/)
2. Update the Makefile with path to `include` and `lib` folders. The specific location will be installation-dependent
3. Update Makefile with compiler (`CC`) and relevant flags (`CFLAGS`)
4. Type `make` in the terminal to compile the program. The executable is `opendda_MPI`

## To Run

1. Open your terminal and change to the folder containing `opendda_MPI`
2. Type
```bash
mpirun -np X opendda_MPI
```

where X is the number of processors available

## Further Reading
To learn more about Open DDA, see the [About OpenDDA](https://github.com/bcrodrigo/OpenDDA_MPI_modifications/blob/main/About%20OpenDDA.md) file.
