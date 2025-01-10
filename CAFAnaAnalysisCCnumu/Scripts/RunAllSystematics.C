// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/ISyst.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TROOT.h"

// std includes.
#include <filesystem>
#include <vector>
#include <memory>
#include <string>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void RunAllSystematics() {
    gROOT->ProcessLine(".L ./Scripts/SelectionSystematics.C");

    std::vector<std::tuple<std::string, int>> SystsVector(XSecSystsVector);
    SystsVector.insert(SystsVector.end(), FluxSystsVector.begin(), FluxSystsVector.end());

    int nSysts = SystsVector.size();
    for (int SystIndex = 0; SystIndex < nSysts; SystIndex++) {
        std::cout << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Starting systematic " << SystIndex << std::endl;
        std::cout << "========================================" << std::endl;

        std::string SystName = std::get<0>(SystsVector.at(SystIndex));
        int SystNUniv = std::get<1>(SystsVector.at(SystIndex));
        if (SystIndex < (int) XSecSystsVector.size()) {
            gROOT->ProcessLine(("SelectionSystematics(\"" + SystName + "\", " + std::to_string(SystNUniv) + ", true)").c_str()); 
        } else {
            gROOT->ProcessLine(("SelectionSystematics(\"" + SystName + "\", " + std::to_string(SystNUniv) + ", false)").c_str()); 
        }

        std::cout << "========================================" << std::endl;
        std::cout << "Done with systematic " << SystIndex << std::endl;
        std::cout << ((SystIndex + 1.) / nSysts) * 100 << "\% done" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << std::endl;
    }
}
