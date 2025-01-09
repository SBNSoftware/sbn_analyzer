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

// Definitions for CutVars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionCutPlots() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionCutPlots.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    // Colors to use
    std::vector<int> Colors;
    Colors.push_back(kBlue+2);
    Colors.push_back(kRed+1);
    Colors.push_back(kOrange+7);
    Colors.push_back(kMagenta+1);

    // Cuts and CutVars to plot
    std::vector<std::vector<Var>> CutVars; std::vector<Binning> CutVarBins; std::vector<Cut> SignalCuts;
    std::vector<TString> PlotNames; std::vector<std::string> CutVarLabels;
    std::vector<std::vector<Cut>> HistosCuts; std::vector<std::vector<TString>> HistosLabels;

    // Nu score before any cut
    CutVars.push_back({kNuScore}); CutVarBins.push_back(bNuScore); SignalCuts.push_back(kNoCut);
    PlotNames.push_back("NoCutNuScore"); CutVarLabels.push_back("#nu score");
    HistosCuts.push_back({kIsNotCosmic, kIsCosmic}); HistosLabels.push_back({"Not cosmic", "Cosmic"});

    // Fmatch score before any cut
    CutVars.push_back({kFMatchScore}); CutVarBins.push_back(bFMatchScore); SignalCuts.push_back(kNoCut);
    PlotNames.push_back("NoCutFMatchScore"); CutVarLabels.push_back("flash matching score");
    HistosCuts.push_back({kIsNotCosmic, kIsCosmic}); HistosLabels.push_back({"Not cosmic", "Cosmic"});

    // Fmatch time before any cut
    CutVars.push_back({kFMatchTime}); CutVarBins.push_back(bFMatchTime); SignalCuts.push_back(kNoCut);
    PlotNames.push_back("NoCutFMatchTime"); CutVarLabels.push_back("flash matching time");
    HistosCuts.push_back({kIsNotCosmic, kIsCosmic}); HistosLabels.push_back({"Not cosmic", "Cosmic"});

    // Muon chi square for muon cut
    CutVars.push_back({kMuMuChi2, kProtonMuChi2, kSecondProtonMuChi2, kPionMuChi2}); CutVarBins.push_back(bMuChi2); SignalCuts.push_back(kCosmicCut);
    PlotNames.push_back("ParticleCutChi2Muon"); CutVarLabels.push_back("#chi^{2}_{#mu}"); 
    HistosCuts.push_back({kHasMuon, kHasProton, kHasSecondProton, kHasPion}); HistosLabels.push_back({"#mu", "p1", "p2", "#pi"});

    // Proton chi square for muon cut
    CutVars.push_back({kMuProtonChi2, kProtonProtonChi2, kSecondProtonProtonChi2, kPionProtonChi2}); CutVarBins.push_back(bProtonChi2); SignalCuts.push_back(kCosmicCut);
    PlotNames.push_back("ParticleCutChi2Proton"); CutVarLabels.push_back("#chi^{2}_{p}"); 
    HistosCuts.push_back({kHasMuon, kHasProton, kHasSecondProton, kHasPion}); HistosLabels.push_back({"#mu", "p1", "p2", "#pi"});

    // Construct all spectra
    std::vector<std::vector<std::unique_ptr<Spectrum>>> Spectra;
    for (std::size_t i = 0; i < CutVars.size(); i++) {
        std::vector<std::unique_ptr<Spectrum>> InnerSpectra;
        for (std::size_t j = 0; j < HistosCuts[i].size(); j++) {
            std::unique_ptr<Spectrum> Signals;
            if (CutVars[i].size() == 1) {
                Signals = std::make_unique<Spectrum>(CutVarLabels[i], CutVarBins[i], NuLoader, CutVars[i][0], kNoSpillCut, SignalCuts[i] && HistosCuts[i][j]);
            } else {
                Signals = std::make_unique<Spectrum>(CutVarLabels[i], CutVarBins[i], NuLoader, CutVars[i][j], kNoSpillCut, SignalCuts[i] && HistosCuts[i][j]);
            }
            InnerSpectra.push_back(std::move(Signals));
        }
        Spectra.push_back(std::move(InnerSpectra));
    }

    NuLoader.Go();

    for (std::size_t i = 0; i < CutVars.size(); i++) {
        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(4);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        std::vector<TH1D*> Histos; Histos.resize(HistosCuts[i].size());
        for (std::size_t j = 0; j < HistosCuts[i].size(); j++) {
            auto& Signals = Spectra[i][j];
            Histos[j] = Signals->ToTH1(TargetPOT);

            TLegendEntry* legReco = leg->AddEntry(Histos[j],HistosLabels[i][j],"l");
            Histos[j]->SetLineColor(Colors[j]);
            Histos[j]->SetLineWidth(4);

            Histos[j]->GetXaxis()->SetTitleFont(FontStyle);
            Histos[j]->GetXaxis()->SetLabelFont(FontStyle);
            Histos[j]->GetXaxis()->SetNdivisions(8);
            Histos[j]->GetXaxis()->SetLabelSize(TextSize);
            Histos[j]->GetXaxis()->SetTitleSize(TextSize);
            Histos[j]->GetXaxis()->SetTitleOffset(1.1);
            Histos[j]->GetXaxis()->CenterTitle();
            Histos[j]->GetXaxis()->SetTitle((CutVarLabels[i]).c_str());

            Histos[j]->GetYaxis()->SetTitleFont(FontStyle);
            Histos[j]->GetYaxis()->SetLabelFont(FontStyle);
            Histos[j]->GetYaxis()->SetNdivisions(6);
            Histos[j]->GetYaxis()->SetLabelSize(TextSize);
            Histos[j]->GetYaxis()->SetTitleSize(TextSize);
            Histos[j]->GetYaxis()->SetTitleOffset(1.3);
            Histos[j]->GetYaxis()->SetTickSize(0);
            Histos[j]->GetYaxis()->CenterTitle();

            double imax = TMath::Max(Histos[j]->GetMaximum(), Histos[0]->GetMaximum());
            double YAxisRange;
            if (PlotNames[i] == "NoCutNuScore") {
                YAxisRange = 1.3*imax;
            } else {
                YAxisRange = 1.15*imax;
            }
            Histos[0]->GetYaxis()->SetRangeUser(0.,YAxisRange);
            Histos[j]->GetYaxis()->SetRangeUser(0.,YAxisRange);

            PlotCanvas->cd();
            Histos[0]->Draw("hist same");
            Histos[j]->Draw("hist same");

            // Save to root file
            SaveFile->WriteObject(Histos[j], PlotNames[i]+HistosLabels[i][j]+"_reco");
        }
        leg->Draw();

        // Save as png
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Cuts/"+PlotNames[i]+".png");

        delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();

}
