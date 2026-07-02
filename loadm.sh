#!/bin/bash

module --force purge
module load anaconda3
module load craype-x86-spr
module load PrgEnv-intel
module load perftools-base/22.12.0
module load netcdf/4.9.2-intel-oneapi-mpi-intel
module load hdf5/1.14.1-2-intel-oneapi-mpi-intel
module load pnetcdf/1.12.2-intel
module load intel-oneapi-mkl/2023.2.0-intel
export OMP_NUM_THREADS=1
export NETCDF=/nopt/nrel/apps/libraries/08-23/spack/opt/spack/linux-rhel8-icelake/intel-2021.10.0/netcdf
export PATH="/nopt/nrel/apps/software/wrf/intel/installs/wps_dep_libs/bin":$PATH
export WRF_EM_CORE=1

module list
