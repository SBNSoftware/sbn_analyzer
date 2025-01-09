source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh

mkdir my_larsoft; cd my_larsoft
setup sbndcode v09_88_00_04 -q e26:prof

mrb newDev -v v09_88_00_04 -q e26:prof
source localProducts_larsoft_v09_88_00_04_e26_prof/setup

cd srcs/
mrb g -t v09_88_00_04 sbndcode

cd $MRB_BUILDDIR
mrbsetenv
mrb i -j4
cd ../
mrbslp

export XDG_CACHE_HOME=/tmp/${USER}
kx509
voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/sbnd/Role=Analysis'