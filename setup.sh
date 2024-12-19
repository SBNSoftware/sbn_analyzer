source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh

#dirName=v09_75_03
#qualifier=e20:prof

dirName=v09_93_01_01
qualifier=e20:prof

mkdir ${dirName}; cd ${dirName}
mrb newDev -v ${dirName} -q ${qualifier}
source localProducts_larsoft_*/setup
cd srcs

# this needs to be fixed when Jaesung & Bruce merge things
mrb g -t v09_75_03 sbnana
cd sbnana
# add your remote; skip if you already have it
git remote add upstream https://github.com/epelaaez/sbnana.git
# pull your feature branch
git pull upstream feature/epelaez_TrueEnsembleSpectrum

# setup
mrbsetenv

# compile
mrb i -j4

# go back to main directory
cd ../../..

kx509
voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/sbnd/Role=Analysis'
export TERM=screen
