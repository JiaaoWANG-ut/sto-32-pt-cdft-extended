#!/bin/bash
# Same MPI/MKL/FFTW/libxc stack as env.sh, but the GPAW extension is the
# Cascade Lake build in ~/software/gpaw-venv-lr6.
GCC=/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-8.5.0/gcc-13.2.0-utdynkvgasyrbya2p7eug6257talnpgc
OMPI=/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/openmpi-4.1.6-fvizfdj24f3gxalmikf2bszpbjgzba76
MKLROOT=/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/intel-oneapi-mkl-2023.2.0-clworlsmp7fmaribct2hwnkslvjayg2k/mkl/2023.2.0
SCALAPACK=/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/netlib-scalapack-2.2.0-qdbdb3gae45vbovb4vuvilccb45awxtv
FFTW=/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-13.2.0/fftw-3.3.10-dgsve3sn52epzjn5cdbzu6lasshydnva
LIBXC=/global/software/rocky-8.x86_64/gcc/linux-rocky8-x86_64/gcc-11.4.0/libxc-6.2.2-2kiqi2c7pmmdyziffcsfci5izv4zny3o
PY=/global/software/rocky-8.x86_64/python-3.11.6/python-3.11.6/python/3.11.6-gcc
VENV=/global/home/users/jiaao/software/gpaw-venv-lr6

export PATH=$VENV/bin:$OMPI/bin:$GCC/bin:$PY/bin:$PATH
export LD_LIBRARY_PATH=$SCALAPACK/lib:$FFTW/lib:$LIBXC/lib:$MKLROOT/lib/intel64:$OMPI/lib:$GCC/lib64:${LD_LIBRARY_PATH:-}
export PKG_CONFIG_PATH=$LIBXC/lib/pkgconfig:$FFTW/lib/pkgconfig:${PKG_CONFIG_PATH:-}
export GPAW_CONFIG=/global/home/users/jiaao/workspace/STO/32-Pt-DEDL-CDFT-extended/siteconfig-lr6.py
export GPAW_MPI_BACKEND=cgpaw
export PMIX_MCA_gds=hash
export UCX_TLS=shm,self
export UCX_SHM_TLS=posix
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -l unlimited 2>/dev/null || true
ulimit -s unlimited 2>/dev/null || true
