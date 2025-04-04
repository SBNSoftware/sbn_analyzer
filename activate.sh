source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh

dirName=v10_04_06_01
cd ${dirName}
source localProducts_larsoft_*/setup

mrbsetenv
cd ..

kx509
voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/sbnd/Role=Analysis'
export TERM=screen
