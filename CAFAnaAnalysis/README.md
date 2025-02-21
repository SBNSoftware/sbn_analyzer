### CAFAna analysis

###CAFAna is used to perform event selection on flat `caf` files. The specific event selection is one muon, two protons, and no pions. 

###Setting up

###To setup your workspace to run the scripts in this directory, you have to run `source setup.sh/activate.sh` as indicated [here](https://github.com/afropapp13/sbn_analyzer/blob/main/README.md).

### Running event selection

###To perform the event selection without any systematics included, you have to run

###```bash
#cafe -bq Scripts/Selection.C
##```

###This will generate the cuts defined in `Definitions.h` and will produce the corresponding figures under `Figs/`. 
###For scripts that do not process the events direclty but rather load the already created histogram and create new plots, you have to run, for example

###```bash
#root -b -q Scripts/SerialPlotGenerator.cpp
###```

###For the systematics study, in which we produce the covariance matrices and error bands for each systematic and plotted variable, we split up the selection code from the script that runs it. 
###This way, we can get results in a stream instead of having to wait many hours for all the results to come at once. 
###The selection logic is stored in `Scripts/SelectionSytematics.C`, and it can be run with all the systematics defined in `Scripts/Definitions.h` by doing

###```bash
#cafe -bq Scripts/RunAllSystematics.C
###```

###There is a `bool ConstructSpectra` flag defined in `Scripts/SelectionSystematics.C`. 
###This is intended to speed up the script after the first time running it. 
###If it is your first time running the systematics study, you will have to set the flag to `true`, and this will store all the histograms needed for the script to run in `.root` files. 
###Then, any time you run the script again (as long as no changes to how the spectra are constructed, e.g., the variables and their binning, are made) you can set the flag to `false` and the runtime of the script will be greatly reduced.

### Available scripts

The scripts availabe to run are:
- `Selection.C`: creates plots for all variables, serialized for double differential variables. 
- `SelectionInteBreakdown.C`: creates plots with event interaction type breakdown.
- `SelectionTopologyBreakdown.C`: creats plots with event topology breakdown.
- `SelectionCutPlots.C`: creates plots for the variables used to perform cuts in our signal definition.
- `SerialPlotGenerator.cpp`: creates sliced plots for the serial double differential histograms.
- `SelectionRunData.C`: saves data of events that pass our reconstructed signal definition into a `.csv` file.
- `SelectionEfficiency.C`: creates signal efficiency plots.
- `SelectionMigrationMatrix.C`: creates migration and response matrices plots.
- `RunAllSystematics.C`: creates covariance matrices for cross-section and flux systematics.
- `SelectionNTargetSystematics.C`: creates covariance matrices for number of target systematics.
- `SelectionPOTSystematics.C`: creates covariance matrices for protons on target systematics.
- `SelectionReinteractionSystematics.C`: creates covariance matrices for reinteraction systematics.
- `SelectionDetectorSystematics.C`: creates covariance matrices for detector systematics.
- `SelectionMCStatSystematics.C`: creates covariance matrices for MC statistical systematics.
- `StatSystematics.cpp`: creates covariance matrices for statistical systematics.
- `TotalCovMatrices.cpp`: creates total covariance matrices using all the systematics above.
- `Unfold.cpp`: unfolds data using the Wiener-SVD method to get cross-sections.
- `WienerSVDOverlay.cpp`: overlays unfolded cross-sections with smeared generator cross-sections from `GeneratorAnalysis/`.
- `SelectionFakeData.C`: creates fake signals to perform fake data studies.
- `UnfoldFakeData.cpp`: performs fake data studies using the Wiener-SVD method.
- `SelectionData.C`: area normalized data/MC plots

###All scripts with name starting with `Selection` are ran using `cafe -bq`, the rest can be ran using `root -b -q` (loading only ROOT is faster). 
###Further, all constant variables that require `sbnana` to be initialized are put in `Definitions.h`, while all other constants that do not require `sbnana` or are usd in pure ROOT scripts, are put in `../Utils/Constants.h`. 
###`Helpers.cpp` contains some helpers functions specific to the CAFAna analysis. 

