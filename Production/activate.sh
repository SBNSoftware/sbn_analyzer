source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh

cd my_larsoft
source localProducts_larsoft_v09_88_00_04_e26_prof/setup
mrbslp
cd ..

export XDG_CACHE_HOME=/tmp/${USER}
kx509
voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/sbnd/Role=Analysis'