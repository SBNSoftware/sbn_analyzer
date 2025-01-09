// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"

// std includes.
#include <vector>
#include <memory>
#include <fstream>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionRunData() {
    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Open csv file to store data
    TString FilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/EventData.csv";
    fstream file; file.open(FilePath, fstream::out); 
    file << "fno,run,subrun,evt,subevt" << std::endl;
    file.close();

    Spectrum sRecoSignals(
        "RecoSignals",
        bEventCount, // does not really matter
        NuLoader,
        kSpillData,
        kNoSpillCut
    );

    NuLoader.Go();
}