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

// Utils includes.
#include "../../Utils/Constants.h"

// Helpers includes.
#include "Helpers.cpp"

using namespace std;
using namespace Constants;

void PowHisto(TH1* hist, double pow) {
    for (int i = 0; i < hist->GetNbinsX() + 1; ++i) {
        double content = hist->GetBinContent(i);
        double powContent = TMath::Power(content, pow);
        hist->SetBinContent(i, powContent);
    }
}

void StatSystematics() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir+"/Figs/CAFAna/Uncertainties/Statistical");

    // Load root file with histograms
    TString HistoFile = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> File(TFile::Open(HistoFile));

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsStats.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    const int NVars = PlotNames.size();

    TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);

    for (int iVar = 0; iVar < NVars; iVar++) {
        // Load true plots
        TH1D* RecoHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_reco"));
        TH1D* RecoTrueHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_reco_true"));
        TH1D* RecoBkgHist = (TH1D*)(File->Get<TH1D>(PlotNames[iVar] + (TString) "_bkg"));

        // Make modified plots
        TH1D* PowRecoHist = (TH1D*) RecoHist->Clone("hnew"); PowHisto(PowRecoHist, 0.5);
        TH1D* PowRecoTrueHist = (TH1D*) RecoTrueHist->Clone("hnew"); PowHisto(PowRecoTrueHist, 0.5);
        TH1D* PowRecoBkgHist = (TH1D*) RecoBkgHist->Clone("hnew"); PowHisto(PowRecoBkgHist, 0.5);

        TH1D* UnivRecoHist = (TH1D*) RecoHist->Clone("hnew"); UnivRecoHist->Add(PowRecoHist, -1);
        TH1D* UnivRecoTrueHist = (TH1D*) RecoTrueHist->Clone("hnew"); UnivRecoTrueHist->Add(PowRecoTrueHist, -1);
        TH1D* UnivRecoBkgHist = (TH1D*) RecoBkgHist->Clone("hnew"); UnivRecoBkgHist->Add(PowRecoBkgHist, -1);

        // Plot histograms with error band
        SelectionHelpers::DrawHistosWithErrorBands(
            RecoHist,
            RecoTrueHist,
            RecoBkgHist,
            {UnivRecoHist},
            {UnivRecoTrueHist},
            {UnivRecoBkgHist},
            dir, 
            "Statistical",
            PlotNames[iVar]
        );

        // Scale histograms for cov matrices
        RecoHist->Scale(Units / (IntegratedFlux * NTargets));
        RecoTrueHist->Scale(Units / (IntegratedFlux * NTargets));
        RecoBkgHist->Scale(Units / (IntegratedFlux * NTargets));

        // Scale modified histograms for cov matrices
        UnivRecoHist->Scale(Units / (IntegratedFlux * NTargets));
        UnivRecoTrueHist->Scale(Units / (IntegratedFlux * NTargets));
        UnivRecoBkgHist->Scale(Units / (IntegratedFlux * NTargets));

        // Plot cov, frac cov, and corr matrices
        std::string RecoAxisTitle = (std::string)RecoHist->GetXaxis()->GetTitle();
        std::string AxisTitle = RecoAxisTitle.substr(5, RecoAxisTitle.size() - 1);

        SelectionHelpers::DrawMatrices(
            RecoHist,
            RecoTrueHist,
            RecoBkgHist,
            {UnivRecoHist},
            {UnivRecoTrueHist},
            {UnivRecoBkgHist},
            dir,
            "Statistical",
            AxisTitle,
            PlotNames[iVar],
            SaveFile
        );
    }
    // Close file
    SaveFile->Close();
}
