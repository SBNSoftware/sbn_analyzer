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
#include <string>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Tools.cxx"
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionDetectorSystematics() {
    Tools tools;

    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // Target files for each sample
    std::string CVPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04/cv";
    std::vector<std::string> CVInputFiles = tools.GetInputFiles(CVPath,false,true);
    auto CVLoader = std::make_unique<SpectrumLoader>(CVInputFiles);

    std::string NoDifPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04/no_diffusion";
    std::vector<std::string> NoDifInputFiles = tools.GetInputFiles(NoDifPath,false,true);
    auto NoDifLoader = std::make_unique<SpectrumLoader>(NoDifInputFiles);

    std::string NoLonDifPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04/no_longitudinal_diffusion";
    std::vector<std::string> NoLonDifInputFiles = tools.GetInputFiles(NoLonDifPath,false, true);
    auto NoLonDifLoader = std::make_unique<SpectrumLoader>(NoLonDifInputFiles);

    std::string NoTraDifPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04/no_transverse_diffusion";
    std::vector<std::string> NoTraDifInputFiles = tools.GetInputFiles(NoTraDifPath,false,true);
    auto NoTraDifLoader = std::make_unique<SpectrumLoader>(NoTraDifInputFiles);

    std::string SCENoDifPath = "/pnfs/sbnd/persistent/users/epelaez/sbnd_detector_variations/v09_88_00_04/with_sce_no_diffusion";
    std::vector<std::string> SCENoDifInputFiles = tools.GetInputFiles(SCENoDifPath,false,true);
    auto SCENoDifLoader = std::make_unique<SpectrumLoader>(SCENoDifInputFiles);

    // Vector with loader for each universe
    std::vector<std::unique_ptr<SpectrumLoader>> UnivLoaders;
    UnivLoaders.push_back(std::move(CVLoader));
    UnivLoaders.push_back(std::move(NoDifLoader));
    UnivLoaders.push_back(std::move(NoLonDifLoader));
    UnivLoaders.push_back(std::move(NoTraDifLoader));
    UnivLoaders.push_back(std::move(SCENoDifLoader));

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir_figs+"/Figs/CAFAna/Uncertainties/Detector");

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsDetector.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    const int NVars = Vars.size();

    // Construct all spectra
    std::vector<std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >>> Spectra;
    for (int iVar = 0; iVar < NVars; ++iVar) {
        std::vector<std::tuple<std::unique_ptr<Spectrum>, std::unique_ptr<Spectrum>, std::unique_ptr<Spectrum>>> InnerSpectra;
        for (int iUniv = 0; iUniv < (int) UnivLoaders.size(); iUniv++) {
            auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), *UnivLoaders.at(iUniv), std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsSignal); 
            auto RecoTrueSignals = std::make_unique<Spectrum> (VarLabels.at(iVar), VarBins.at(iVar), *UnivLoaders.at(iUniv), std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsTrueReco); 
            auto RecoBkgSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), *UnivLoaders.at(iUniv), std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsBackground); 
            InnerSpectra.push_back({std::move(RecoSignals), std::move(RecoTrueSignals), std::move(RecoBkgSignals)});
        }
        Spectra.push_back(std::move(InnerSpectra));
    }
    for (int iUniv = 0; iUniv < (int) UnivLoaders.size(); iUniv++) {
        UnivLoaders.at(iUniv)->Go();
    }

    for (int iVar = 0; iVar < NVars; ++iVar) {
        TH1D* RecoHisto = std::get<0>(Spectra.at(iVar).at(0))->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = std::get<1>(Spectra.at(iVar).at(0))->ToTH1(TargetPOT);
        TH1D* RecoBkgHisto = std::get<2>(Spectra.at(iVar).at(0))->ToTH1(TargetPOT);

        std::vector<TH1D*> UnivRecoHisto;
        std::vector<TH1D*> UnivRecoTrueHisto;
        std::vector<TH1D*> UnivRecoBkgHisto;
        for (int iUniv = 1; iUniv < (int) UnivLoaders.size(); iUniv++) {
            UnivRecoHisto.push_back(std::get<0>(Spectra.at(iVar).at(iUniv))->ToTH1(TargetPOT));
            UnivRecoTrueHisto.push_back(std::get<1>(Spectra.at(iVar).at(iUniv))->ToTH1(TargetPOT));
            UnivRecoBkgHisto.push_back(std::get<2>(Spectra.at(iVar).at(iUniv))->ToTH1(TargetPOT));
        }

        SelectionHelpers::DrawHistosWithErrorBands(
            RecoHisto,
            RecoTrueHisto,
            RecoBkgHisto,
            UnivRecoHisto,
            UnivRecoTrueHisto,
            UnivRecoBkgHisto,
            dir_figs,
            "Detector",
            PlotNames[iVar]
        );

        // Scale histograms for cov matrices
        RecoHisto->Scale(Units / (IntegratedFlux * NTargets));
        RecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
        RecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

        // Scale modified histograms for cov matrices
        for (int iUniv = 0; iUniv < (int) UnivRecoHisto.size(); iUniv++) {
            UnivRecoHisto.at(iUniv)->Scale(Units / (IntegratedFlux * NTargets));
            UnivRecoTrueHisto.at(iUniv)->Scale(Units / (IntegratedFlux * NTargets));
            UnivRecoBkgHisto.at(iUniv)->Scale(Units / (IntegratedFlux * NTargets));
        }

        // Plot cov, frac cov, and corr matrices
        SelectionHelpers::DrawMatrices(
            RecoHisto,
            RecoTrueHisto,
            RecoBkgHisto,
            UnivRecoHisto,
            UnivRecoTrueHisto,
            UnivRecoBkgHisto,
            dir_figs,
            "Detector",
            VarLabels[iVar],
            PlotNames[iVar],
            SaveFile
        );
    }
    // Close file
    SaveFile->Close();
}
