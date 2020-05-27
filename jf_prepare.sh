#module load cmake/3.7.1
#module load intel/18.0.2
#module load gsl
#module load boost
#module load hdf5
#module load impi
#module load eigen

#export EIGEN_DOWNLOAD_DIR=$HOME/Software
export EIGEN_INSTALL_DIR=~/runs/sims/run_single_event/eigen/
export EIGEN3_ROOT=~/runs/sims/run_single_event/eigen/
export GSL=$(gsl-config --prefix)
export GSL_HOME=$(gsl-config --prefix)
export GSL_ROOT_DIR=$(gsl-config --prefix)
export JETSCAPE_DIR=`readlink -f .`

export SMASH_DIR=${JETSCAPE_DIR}/external_packages/smash/smash_code
export number_of_cores=`nproc --all`

export PYTHIAINSTALLDIR=~/runs/sims/run_single_event/

export PYTHIA8DIR=${PYTHIAINSTALLDIR}/pythia8235
export PYTHIA8_ROOT_DIR=${PYTHIAINSTALLDIR}/pythia8235

export CC=gcc
export CXX=g++
export OpenMP_CXX=g++

#cd external_packages; bash get_music.sh ; bash get_freestream-milne.sh ; bash get_iSS.sh ; bash get_smash.sh
#mkdir  build; cd build; cmake cmake -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_SMASH=ON ..; make -j
