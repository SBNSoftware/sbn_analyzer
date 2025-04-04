source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh

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

#####################################

# creating the relevant directories

me="$(whoami)"
fig_dir="/exp/sbnd/data/users/${me}"

mkdir ${fig_dir}/cc2p
mkdir ${fig_dir}/cc2p/Overlay
mkdir ${fig_dir}/cc2p/Overlay/PreFSI
mkdir ${fig_dir}/cc2p/Overlay/PostFSI
mkdir ${fig_dir}/cc2p/Figs
mkdir ${fig_dir}/cc2p/Figs/Overlay
mkdir ${fig_dir}/cc2p/Figs/Overlay/Serial
mkdir ${fig_dir}/cc2p/Figs/Overlay/MEC
mkdir ${fig_dir}/cc2p/Figs/Overlay/MEC/Serial
mkdir ${fig_dir}/cc2p/Figs/InteBreakDown
mkdir ${fig_dir}/cc2p/Figs/InteBreakDown/PostFSI
mkdir ${fig_dir}/cc2p/Figs/InteBreakDown/PreFSI
mkdir ${fig_dir}/cc2p/Figs/CAFAna
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Cuts
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Smear
mkdir ${fig_dir}/cc2p/Figs/CAFAna/InteBreakdown
mkdir ${fig_dir}/cc2p/Figs/CAFAna/TopologyBreakdown
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Serial
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Efficiency
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Matrices
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Unfolded
mkdir ${fig_dir}/cc2p/Figs/CAFAna/WienerSVDOverlay
mkdir ${fig_dir}/cc2p/Figs/CAFAna/BinUncertainties
mkdir ${fig_dir}/cc2p/Figs/CAFAna/UnfBinUncertainties
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/NTargets
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/POT
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/Reinteraction
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/Detector
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/MCStat
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/Statistical
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/BinUncertainties
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/Unfolded
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/WienerSVDOverlay
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Uncertainties/FakeDataStudies
mkdir ${fig_dir}/cc2p/Figs/CAFAna/FakeDataStudies
mkdir ${fig_dir}/cc2p/Figs/CAFAna/FakeDataStudies/TwiceMEC
mkdir ${fig_dir}/cc2p/Figs/CAFAna/FakeDataStudies/Combined
mkdir ${fig_dir}/cc2p/Figs/CAFAna/DataCounts
mkdir ${fig_dir}/cc2p/Figs/CAFAna/DataCounts/NoCut
mkdir ${fig_dir}/cc2p/Figs/CAFAna/Resolution
mkdir ${fig_dir}/cc2p/Figs/CAFAna/CutScan