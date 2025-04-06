source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh

dirName=v10_04_06_01
cd ${dirName}
source localProducts_larsoft_*/setup

mrbsetenv
cd ..

kx509
httokensh -a htvaultprod.fnal.gov -i sbnd -- /bin/bash

export TERM=screen
