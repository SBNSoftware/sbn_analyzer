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

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

// Helpers includes.
#include "Helpers.cpp"

using namespace std;
using namespace ana;
using namespace Constants;

void PowHisto(TH1* hist, double pow, bool abs) {
    for (int i = 0; i < hist->GetNbinsX() + 1; ++i) {
        double content = hist->GetBinContent(i);
        if (abs == true) content = TMath::Abs(content);
        double powContent = TMath::Power(content, pow);
        if ((abs == true) && (hist->GetBinContent(i) < 0)) powContent = -powContent; 
        hist->SetBinContent(i, powContent);
    }
}

void SelectionFakeData() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionFake.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    // Construct all spectra
    std::vector<std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >>> Spectra;
    for (std::size_t iData = 0; iData < FakeWeights.size(); ++iData) {
        std::vector<std::tuple<std::unique_ptr<Spectrum>, std::unique_ptr<Spectrum>>> InnerSpectra;
        for (std::size_t iVar = 0; iVar < Vars.size(); ++iVar) {
            auto FakeReco = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<0>(Vars.at(iVar)), kRecoIsSignal, kNoShift, std::get<0>(FakeWeights.at(iData))); 
            auto FakeTrue = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<2>(Vars.at(iVar)), kTruthIsSignal, kNoSpillCut, kNoShift, std::get<1>(FakeWeights.at(iData)));
            InnerSpectra.push_back({std::move(FakeReco), std::move(FakeTrue)});
        }
        Spectra.push_back(std::move(InnerSpectra));
    }

    NuLoader.Go();

    // Loop over all fake data 
    for (std::size_t iData = 0; iData < FakeWeights.size(); ++iData) {
        for (std::size_t iVar = 0; iVar < Vars.size(); ++iVar) {
            TH1D* FakeRecoHisto = std::get<0>(Spectra[iData][iVar])->ToTH1(TargetPOT);
            TH1D* FakeTrueHisto = std::get<1>(Spectra[iData][iVar])->ToTH1(TargetPOT);

            // Manage under/overflow bins
            FakeRecoHisto->SetBinContent(FakeRecoHisto->GetNbinsX(), FakeRecoHisto->GetBinContent(FakeRecoHisto->GetNbinsX()) + FakeRecoHisto->GetBinContent(FakeRecoHisto->GetNbinsX() + 1));
            FakeTrueHisto->SetBinContent(FakeTrueHisto->GetNbinsX(), FakeTrueHisto->GetBinContent(FakeTrueHisto->GetNbinsX()) + FakeTrueHisto->GetBinContent(FakeTrueHisto->GetNbinsX() + 1));

            FakeRecoHisto->SetBinContent(1, FakeRecoHisto->GetBinContent(0) + FakeRecoHisto->GetBinContent(1));
            FakeTrueHisto->SetBinContent(1, FakeTrueHisto->GetBinContent(0) + FakeTrueHisto->GetBinContent(1));

            // Save histograms to file
            SaveFile->WriteObject(FakeRecoHisto, PlotNames[iVar] + FakeDataNames[iData] + "_reco");
            SaveFile->WriteObject(FakeTrueHisto, PlotNames[iVar] + FakeDataNames[iData] + "_true");

            int n = FakeRecoHisto->GetXaxis()->GetNbins();
            double edges[n+1];
            for (int i = 0; i < n+1; i++) { edges[i] = FakeRecoHisto->GetBinLowEdge(i+1); }

            ////////////////////////
            // Get stat cov matrices
            ////////////////////////

            // Make modified plots
            TH1D* PowRecoHist = (TH1D*) FakeRecoHisto->Clone("hnew"); PowHisto(PowRecoHist, 0.5, true);
            TH1D* UnivRecoHist = (TH1D*) FakeRecoHisto->Clone("hnew"); UnivRecoHist->Add(PowRecoHist, -1);

            // Create covariance matrix
            std::string CovName = "Cov" + (std::string) FakeDataNames[iData];
            TH2* CovMatrix = new TH2D(
                (CovName + (std::string)PlotNames[iVar]).c_str(),
                CovName.c_str(),
                n, edges,
                n, edges
            );

            // Scale histograms
            FakeRecoHisto->Scale(Units / (IntegratedFlux * NTargets));
            UnivRecoHist->Scale(Units / (IntegratedFlux * NTargets));

            // Populate cov matrix
            SelectionHelpers::GetCovMatrix(
                FakeRecoHisto,
                {UnivRecoHist},
                CovMatrix,
                n
            );

            // Save cov matrix
            SaveFile->WriteObject(CovMatrix, PlotNames[iVar] + FakeDataNames[iData] + "_stat_cov");
        }
    }
}
