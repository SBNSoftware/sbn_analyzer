// Root includes.
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TMatrixD.h>

// std includes.
#include <vector>

// Helpers includes.
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Util.C"
#include "../../Utils/Constants.h"

using namespace std;
using namespace Constants;

void TotalCovMatrices() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.1f");

    // Vector with all systematic files
    std::vector<std::unique_ptr<TFile>> CovFiles;

    // Add xsec systematics
    for (int iSyst = 0; iSyst < (int) XSecSystsVector.size(); iSyst++) {
        std::string SystName = std::get<0>(XSecSystsVector.at(iSyst));
        TString FilePath =  "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematics"+TString(SystName)+".root";
        std::unique_ptr<TFile> File(TFile::Open(FilePath));
        CovFiles.push_back(std::move(File));
    }

    // Add flux systematics
    for (int iSyst = 0; iSyst < (int) FluxSystsVector.size(); iSyst++) {
        std::string SystName = std::get<0>(FluxSystsVector.at(iSyst));
        TString FilePath =  "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematics"+TString(SystName)+".root";
        std::unique_ptr<TFile> File(TFile::Open(FilePath));
        CovFiles.push_back(std::move(File));
    }

    // Add stat systematics
    std::unique_ptr<TFile> StatsFile(TFile::Open("/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsStats.root"));
    CovFiles.push_back(std::move(StatsFile));

    // Add POT systematics
    std::unique_ptr<TFile> POTFile(TFile::Open("/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsPOT.root"));
    CovFiles.push_back(std::move(POTFile));

    // Add NTargets systematics
    std::unique_ptr<TFile> NTargetsFile(TFile::Open("/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsNTargets.root"));
    CovFiles.push_back(std::move(NTargetsFile));

    // Add Detector systematics
    std::unique_ptr<TFile> DetectorFile(TFile::Open("/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsDetector.root"));
    CovFiles.push_back(std::move(DetectorFile));

    // Add Reinteraction systematics
    std::unique_ptr<TFile> ReinteractionFile(TFile::Open("/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsReinteraction.root"));
    CovFiles.push_back(std::move(ReinteractionFile));

    // Add MCStat systematics
    std::unique_ptr<TFile> MCStatFile(TFile::Open("/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsMCStat.root"));
    CovFiles.push_back(std::move(MCStatFile));

    // File with reco histograms
    TString HistoFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> HistoFile(TFile::Open(HistoFilePath));

    // File to store total cov matrices
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/TotalCovMatrices.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    TCanvas* PlotCanvas = new TCanvas("Cov","Cov",205,34,1124,768);

    // Loop over all vars
    for (int iVar = 0; iVar < (int) PlotNames.size(); iVar++) {
        // Margins for matrices
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        // Get first matrix to get correct dimensions
        TH2D* FirstCovHist = (TH2D*)(CovFiles[0]->Get<TH2D>(PlotNames[iVar]+"_cov"));
        int n = FirstCovHist->GetXaxis()->GetNbins();
        double edges[n+1];
        for (int i = 0; i < n+1; i++) { edges[i] = FirstCovHist->GetXaxis()->GetBinLowEdge(i+1); }
        TMatrixD TotalCovMatrix(n, n); H2M(FirstCovHist, TotalCovMatrix, kTRUE);

        // Add all uncertainties
        for (int iCov = 1; iCov < (int) CovFiles.size(); ++iCov) {
            TH2D* CovHist = (TH2D*)(CovFiles[iCov]->Get<TH2D>(PlotNames[iVar]+"_cov"));
            TMatrixD CovMatrix(n, n); H2M(CovHist, CovMatrix, kTRUE);
            TotalCovMatrix += CovMatrix;
        }
        
        TH2D* TotalCovHist = new TH2D(
            "TotalCov" + PlotNames[iVar],
            "TotalCov" + PlotNames[iVar],
            n, edges, 
            n, edges
        );
        M2H(TotalCovMatrix, TotalCovHist);

        double CovMin = TotalCovHist->GetMinimum();
        double CovMax = TotalCovHist->GetMaximum();
        TotalCovHist->GetZaxis()->SetRangeUser(CovMin,CovMax);
        TotalCovHist->GetZaxis()->CenterTitle();
        TotalCovHist->GetZaxis()->SetTitleFont(FontStyle);
        TotalCovHist->GetZaxis()->SetTitleSize(TextSize);
        TotalCovHist->GetZaxis()->SetLabelFont(FontStyle);
        TotalCovHist->GetZaxis()->SetLabelSize(TextSize);
        TotalCovHist->GetZaxis()->SetNdivisions(5);
        TotalCovHist->GetXaxis()->SetTitle(FirstCovHist->GetXaxis()->GetTitle());
        TotalCovHist->GetYaxis()->SetTitle(FirstCovHist->GetYaxis()->GetTitle());

        PlotCanvas->cd();
        TotalCovHist->Draw("colz");
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Uncertainties/TotalCov"+PlotNames[iVar]+".pdf");

        SaveFile->WriteObject(TotalCovHist, PlotNames[iVar]);

        // Get xsec cov matrix
        TMatrixD XSecFracCov(n, n);
        for (int iSyst = 0; iSyst < (int) XSecSystsVector.size(); iSyst++) {
            TH2D* FracCovHist = (TH2D*)(CovFiles[iSyst]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
            TMatrixD FracCovMatrix(n, n); H2M(FracCovHist, FracCovMatrix, kTRUE);
            XSecFracCov += FracCovMatrix;
        }

        // Get flux cov matrix
        TMatrixD FluxFracCov(n, n);
        for (int iSyst = 0; iSyst < (int) FluxSystsVector.size(); iSyst++) {
            TH2D* FracCovHist = (TH2D*)(CovFiles[iSyst + (int)XSecSystsVector.size()]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
            TMatrixD FracCovMatrix(n, n); H2M(FracCovHist, FracCovMatrix, kTRUE);
            FluxFracCov += FracCovMatrix;
        }

        int offset = XSecSystsVector.size() + FluxSystsVector.size();

        // Get stat cov matrix
        TMatrixD StatFracCov(n, n);
        TH2D* StatFracCovHist = (TH2D*)(CovFiles[offset]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
        H2M(StatFracCovHist, StatFracCov, kTRUE);

        // Get POT cov matrix
        TMatrixD POTFracCov(n, n);
        TH2D* POTFracCovHist = (TH2D*)(CovFiles[offset + 1]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
        H2M(POTFracCovHist, POTFracCov, kTRUE);

        // Get NTargets cov matrix
        TMatrixD NTargetsFracCov(n, n);
        TH2D* NTargetsFracCovHist = (TH2D*)(CovFiles[offset + 2]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
        H2M(NTargetsFracCovHist, NTargetsFracCov, kTRUE);

        // Get Detector cov matrix
        TMatrixD DetectorFracCov(n, n);
        TH2D* DetectorFracCovHist = (TH2D*)(CovFiles[offset + 3]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
        H2M(DetectorFracCovHist, DetectorFracCov, kTRUE);

        // Get Reinteraction cov matrix
        TMatrixD ReinteractionFracCov(n, n);
        TH2D* ReinteractionFracCovHist = (TH2D*)(CovFiles[offset + 4]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
        H2M(ReinteractionFracCovHist, ReinteractionFracCov, kTRUE);

        // Get MCstat cov matrix
        TMatrixD MCStatFracCov(n, n);
        TH2D* MCStatFracCovHist = (TH2D*)(CovFiles[offset + 5]->Get<TH2D>(PlotNames[iVar]+"_fraccov"));
        H2M(MCStatFracCovHist, MCStatFracCov, kTRUE);

        // Histograms for uncertainties
        TH1D* XSecHisto = new TH1D("XSec"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* FluxHisto = new TH1D("Flux"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* StatsHisto = new TH1D("Stats"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* POTHisto = new TH1D("POT"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* NTargetsHisto = new TH1D("NTargets"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* DetectorHisto = new TH1D("Detector"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* ReinteractionHisto = new TH1D("Reinteraction"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* MCStatHisto = new TH1D("MCStat"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* TotalHisto = new TH1D("Total"+PlotNames[iVar], ";;Uncertainty [%]", n, edges);
        TH1D* RecoHist = (TH1D*)(HistoFile->Get<TH1D>(PlotNames[iVar] + (TString) "_reco"));

        // Histograms to store total frac cov and corr
        TH2D* TotalFracCovHisto = new TH2D("TotalFracCov"+PlotNames[iVar],"TotalFracCov" + PlotNames[iVar],n, edges, n, edges);
        TH2D* TotalCorrHisto = new TH2D("TotalCorr"+PlotNames[iVar],"TotalCorr" + PlotNames[iVar],n, edges, n, edges);
        RecoHist->Scale(Units / (IntegratedFlux * NTargets));
        SelectionHelpers::GetFracCovAndCorrMatrix(RecoHist, TotalCovHist, TotalFracCovHisto, TotalCorrHisto, n);
        TMatrixD TotalFracCov(n, n); H2M(TotalFracCovHisto, TotalFracCov, kTRUE);

        // Loop over each bin
        for (int iBin = 0; iBin < n; iBin++) {
            double Total = 0.;

            // Add xsec uncertainties
            XSecHisto->SetBinContent(iBin + 1, TMath::Sqrt(XSecFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(XSecFracCov(iBin, iBin)) * 100, 2);

            // Add flux uncertainties
            FluxHisto->SetBinContent(iBin + 1, TMath::Sqrt(FluxFracCov(iBin, iBin)) * 100);
            Total +=  TMath::Power(TMath::Sqrt(FluxFracCov(iBin, iBin)) * 100, 2);

            // Add Stat systematics
            StatsHisto->SetBinContent(iBin + 1, TMath::Sqrt(StatFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(StatFracCov(iBin, iBin)) * 100, 2);

            // Add POT systematics
            POTHisto->SetBinContent(iBin + 1, TMath::Sqrt(POTFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(POTFracCov(iBin, iBin)) * 100, 2);

            // Add NTargets systematics
            NTargetsHisto->SetBinContent(iBin + 1, TMath::Sqrt(NTargetsFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(NTargetsFracCov(iBin, iBin)) * 100, 2);

            // Add Detector systematics
            DetectorHisto->SetBinContent(iBin + 1, TMath::Sqrt(DetectorFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(DetectorFracCov(iBin, iBin)) * 100, 2);

            // Add Reinteraction systematics
            ReinteractionHisto->SetBinContent(iBin + 1, TMath::Sqrt(ReinteractionFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(ReinteractionFracCov(iBin, iBin)) * 100, 2);

            // Add MCStat systematics
            MCStatHisto->SetBinContent(iBin + 1, TMath::Sqrt(MCStatFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(MCStatFracCov(iBin, iBin)) * 100, 2);

            // Total 
            TotalHisto->SetBinContent(iBin + 1, TMath::Sqrt(TotalFracCov(iBin, iBin)) * 100);

            // These should be the same
            std::cout << TMath::Sqrt(TotalFracCov(iBin, iBin)) * 100 << " ==? ";
            std::cout << TMath::Sqrt(Total) << std::endl;;
        }
        // Margins for bin by bin uncertainties plot
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legXSec = leg->AddEntry(XSecHisto,"XSec","l");
        XSecHisto->SetLineColor(kBlue+2);
        XSecHisto->SetLineWidth(4);
        XSecHisto->SetMarkerSize(1.5);
        XSecHisto->SetMarkerColor(kBlue+2);

        // Style 
        XSecHisto->GetXaxis()->SetNdivisions(5);
        if (PlotNames[iVar] == "EventCount") {
            XSecHisto->GetXaxis()->SetLabelSize(0);
            XSecHisto->GetXaxis()->SetTitleSize(0);
        } else {
            XSecHisto->GetXaxis()->SetTitleFont(FontStyle);
            XSecHisto->GetXaxis()->SetLabelFont(FontStyle);
            XSecHisto->GetXaxis()->SetLabelSize(TextSize);
            XSecHisto->GetXaxis()->SetTitleSize(TextSize);
            XSecHisto->GetXaxis()->SetTitleOffset(1.);
            XSecHisto->GetXaxis()->CenterTitle();
            XSecHisto->GetXaxis()->SetTitle(RecoHist->GetXaxis()->GetTitle());
        }

        XSecHisto->GetYaxis()->SetTitleFont(FontStyle);
        XSecHisto->GetYaxis()->SetLabelFont(FontStyle);
        XSecHisto->GetYaxis()->SetLabelSize(TextSize);
        XSecHisto->GetYaxis()->SetTitleSize(TextSize);
        XSecHisto->GetYaxis()->SetNdivisions(6);
        XSecHisto->GetYaxis()->SetTitleOffset(1.);
        XSecHisto->GetYaxis()->SetTickSize(0);
        XSecHisto->GetYaxis()->CenterTitle();
        XSecHisto->GetYaxis()->SetRangeUser(0.,TotalHisto->GetMaximum()*1.35);
        XSecHisto->Draw("hist text0");

        TLegendEntry* legFlux = leg->AddEntry(FluxHisto,"Flux","l");
        FluxHisto->SetLineColor(kRed+1);
        FluxHisto->SetLineWidth(4);
        FluxHisto->SetMarkerSize(1.5);
        FluxHisto->SetMarkerColor(kRed+1);
        FluxHisto->Draw("hist text0 same");

        TLegendEntry* legStats= leg->AddEntry(StatsHisto,"Stat","l");
        StatsHisto->SetLineColor(kMagenta+1);
        StatsHisto->SetLineWidth(4);
        StatsHisto->SetMarkerSize(1.5);
        StatsHisto->SetMarkerColor(kMagenta+1);
        StatsHisto->Draw("hist text0 same");

        TLegendEntry* legPOT= leg->AddEntry(POTHisto,"POT","l");
        POTHisto->SetLineColor(kGreen);
        POTHisto->SetLineWidth(4);
        POTHisto->SetMarkerSize(1.5);
        POTHisto->SetMarkerColor(kGreen);
        POTHisto->Draw("hist text0 same");

        TLegendEntry* legNTargets= leg->AddEntry(NTargetsHisto,"NTargets","l");
        NTargetsHisto->SetLineColor(kYellow+1);
        NTargetsHisto->SetLineWidth(4);
        NTargetsHisto->SetMarkerSize(1.5);
        NTargetsHisto->SetMarkerColor(kYellow+1);
        NTargetsHisto->Draw("hist text0 same");

        TLegendEntry* legDetector= leg->AddEntry(DetectorHisto,"Detector","l");
        DetectorHisto->SetLineColor(kCyan-3);
        DetectorHisto->SetLineWidth(4);
        DetectorHisto->SetMarkerSize(1.5);
        DetectorHisto->SetMarkerColor(kCyan-3);
        DetectorHisto->Draw("hist text0 same");

        TLegendEntry* legReinteraction= leg->AddEntry(ReinteractionHisto,"Reinteraction","l");
        ReinteractionHisto->SetLineColor(kTeal + 3);
        ReinteractionHisto->SetLineWidth(4);
        ReinteractionHisto->SetMarkerSize(1.5);
        ReinteractionHisto->SetMarkerColor(kTeal+3);
        ReinteractionHisto->Draw("hist text0 same");

        TLegendEntry* legMCStat= leg->AddEntry(MCStatHisto,"MCStat","l");
        MCStatHisto->SetLineColor(kSpring+9);
        MCStatHisto->SetLineWidth(4);
        MCStatHisto->SetMarkerSize(1.5);
        MCStatHisto->SetMarkerColor(kSpring+9);
        MCStatHisto->Draw("hist text0 same");

        TLegendEntry* legTotal= leg->AddEntry(TotalHisto,"Total","l");
        TotalHisto->SetLineColor(kBlack);
        TotalHisto->SetLineWidth(4);
        TotalHisto->SetMarkerSize(1.5);
        TotalHisto->SetMarkerColor(kBlack);
        TotalHisto->Draw("hist text0 same");

        // Save plot
        leg->Draw();
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/BinUncertainties/"+PlotNames[iVar]+".pdf");
    }
}
