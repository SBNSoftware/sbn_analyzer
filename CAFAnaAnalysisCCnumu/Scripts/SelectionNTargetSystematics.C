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
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"

// std includes.
#include <filesystem>
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionNTargetSystematics() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir_figs+"/Figs/CAFAna/Uncertainties/NTargets");

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsNTargets.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    const int NVars = Vars.size();

    // Construct all spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >> Spectra;
    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsSignal); 
        auto RecoTrueSignals = std::make_unique<Spectrum> (VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsTrueReco); 
        auto RecoBkgSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsBackground); 
        Spectra.push_back({std::move(RecoSignals), std::move(RecoTrueSignals), std::move(RecoBkgSignals)});
    }
    NuLoader.Go();

    // Modified NTargets by 1%
    double ModifiedNTargets = NTargets * 1.01;

    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto& [RecoSpectra, RecoTrueSpectra, RecoBkgSpectra] = Spectra.at(iVar);

        // Get nominal plots
        TH1D* RecoHisto = RecoSpectra->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSpectra->ToTH1(TargetPOT);
        TH1D* RecoBkgHisto = RecoBkgSpectra->ToTH1(TargetPOT);

        // Get plots with modified POT
        TH1D* ModifiedRecoHisto = RecoSpectra->ToTH1(TargetPOT);
        TH1D* ModifiedRecoTrueHisto = RecoTrueSpectra->ToTH1(TargetPOT);
        TH1D* ModifiedRecoBkgHisto = RecoBkgSpectra->ToTH1(TargetPOT);

        // Scale histograms for cov matrices
        RecoHisto->Scale(Units / (IntegratedFlux * NTargets));
        RecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
        RecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

        // Scale modified histograms for cov matrices
        ModifiedRecoHisto->Scale(Units / (IntegratedFlux * ModifiedNTargets));
        ModifiedRecoTrueHisto->Scale(Units / (IntegratedFlux * ModifiedNTargets));
        ModifiedRecoBkgHisto->Scale(Units / (IntegratedFlux * ModifiedNTargets));

        // Plot cov, frac cov, and corr matrices
        SelectionHelpers::DrawMatrices(
            RecoHisto,
            RecoTrueHisto,
            RecoBkgHisto,
            {ModifiedRecoHisto},
            {ModifiedRecoTrueHisto},
            {ModifiedRecoBkgHisto},
            dir_figs,
            "NTargets",
            VarLabels[iVar],
            PlotNames[iVar],
            SaveFile
        );
    }
    // Close file
    SaveFile->Close();
}
