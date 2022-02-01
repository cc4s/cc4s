module load gcc/9.1.0-gcc-4.8.5-mj7s6dg
module load openmpi/4.1.1-gcc-9.1.0-fwmcvon
module load openblas/0.3.6-gcc-9.1.0-36c6y5u
module load netlib-lapack/3.8.0-gcc-9.1.0-z5qzrho

module load cmake/3.15.1-gcc-9.1.0-4jglial
module load python/3.6.6-gcc-9.1.0-mkqdc62

module load openblas/0.3.6-gcc-9.1.0-36c6y5u

export BLAS_PATH=$( module show openblas/0.3.6-gcc-9.1.0-36c6y5u |
                    awk -F'LD_LIBRARY_PATH' '{print $2}' |
                    xargs -n1 dirname
                    )

