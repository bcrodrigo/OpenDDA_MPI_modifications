
OpenDDA is a highly optimised computational framework, written in the C language, for the Discrete Dipole Approximation, a numerical method for calculating the optical properties associated with a target of arbitrary geometry that is widely used in atmospheric, astrophysical and industrial simulations.

Core optimisations include the bit-fielding of integer data and iterative methods that complement a new Discrete Fourier Transform (DFT) kernel, which efficiently calculates the matrix-vector products required by these iterative solution schemes. The new kernel performs the requisite 3D DFTs as ensembles of 1D transforms, and by doing so, is able to reduce the number of constituent 1D transforms by 60% and the memory by over 80%. The optimisations also facilitate the use of parallel techniques to further enhance the performance. Complete OpenMP-based shared-memory and MPI-based distributed-memory implementations have been created to take full advantage of the various architectures.

OpenDDA is free to use for everyone. The only request pertaining to its use is that any subsequent work, use, modification and / or augmentation of the available software, should cite the associated publication, see below for details.

The are two associated documents that are useful when working with the OpenDDA framework:

1. The associated publication:
OpenDDA: A Novel High-Performance Computational Framework for the Discrete Dipole Approximation, Mc Donald, J., Golden, A., Jennings, S. G., International Journal of High Performance Computing Applications (IJHPCA), Volume 23, No. 1, 42-61, 2009.
http://hpc.sagepub.com/cgi/content/abstract/23/1/42

2. The associated Ph.D. Thesis:
OpenDDA: A Novel High-Performance Computational Framework for the Discrete Dipole Approximation, Mc Donald, J., Ph.D. Thesis, School of Physics, National University of Ireland, Galway, Ireland, September 2007. https://github.com/drjmcdonald/OpenDDA/blob/main/thesis_phd_OpenDDA_2007.pdf

