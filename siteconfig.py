# Build GPAW against the lr8 stack. -march=znver4 matches the lr8 EPYC nodes.
# Do not build this on a login node with -march=native.
mpi = True
compiler = 'mpic++'
fftw = True
scalapack = True
libxc = '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-11.4.0/libxc-6.2.2-2kiqi2c7pmmdyziffcsfci5izv4zny3o'

library_dirs = [
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/fftw-3.3.10-dgsve3sn52epzjn5cdbzu6lasshydnva/lib',
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-11.4.0/libxc-6.2.2-2kiqi2c7pmmdyziffcsfci5izv4zny3o/lib',
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/netlib-scalapack-2.2.0-qdbdb3gae45vbovb4vuvilccb45awxtv/lib',
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/intel-oneapi-mkl-2023.2.0-clworlsmp7fmaribct2hwnkslvjayg2k/mkl/2023.2.0/lib/intel64',
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/openmpi-4.1.6-fvizfdj24f3gxalmikf2bszpbjgzba76/lib',
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/gcc-13.2.0-utdynkvgasyrbya2p7eug6257talnpgc/lib64',
]
include_dirs = [
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/fftw-3.3.10-dgsve3sn52epzjn5cdbzu6lasshydnva/include',
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-11.4.0/libxc-6.2.2-2kiqi2c7pmmdyziffcsfci5izv4zny3o/include',
    '/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/openmpi-4.1.6-fvizfdj24f3gxalmikf2bszpbjgzba76/include',
]
# Leave runtime_library_dirs empty. GPAW 26.7's Makefile generator crashes
# when rpath flags come back as a list. Jobs pick the libs up from
# LD_LIBRARY_PATH in env.sh.

# C BLAS interface. Sequential MKL because the job is pure MPI (OMP=1).
# Replace the default ['xc'] list. MKL provides the Fortran BLAS names GPAW calls.
libraries = [
    'xc',
    'fftw3',
    'scalapack',
    'mkl_intel_lp64',
    'mkl_sequential',
    'mkl_core',
    'pthread',
    'dl',
    'm',
    'gfortran',
]

extra_compile_args += ['-O3', '-march=znver4']
