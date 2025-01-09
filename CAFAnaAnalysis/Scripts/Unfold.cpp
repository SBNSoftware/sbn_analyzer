// ROOT includes.
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TGraphAsymmErrors.h"

// std includes.
#include <vector>

// Helpers includes.
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Tools.cxx"
#include "../../Utils/Util.C"
#include "../../Utils/WienerSVD.C"
#include "../../Utils/Constants.h"

using namespace std;
using namespace Constants;

void Unfold() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.1f");

    Tools tools;

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Unfolded.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    // Load root file with histograms
    TString SelectionRootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> SelectionFile(TFile::Open(SelectionRootFilePath));

    // Load root file with response matrices
    TString MatrixRootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Matrix.root";
    std::unique_ptr<TFile> MatrixFile(TFile::Open(MatrixRootFilePath));

    // Load root file with total covariance matricex
    TString CovRootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/TotalCovMatrices.root";
    std::unique_ptr<TFile> CovFile(TFile::Open(CovRootFilePath));

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

    const int NPlots = PlotNames.size();

    for (int iPlot = 0; iPlot < NPlots; iPlot++) {
        // Load necessary plots
        TH2D* ResponseHist = (TH2D*)(MatrixFile->Get<TH2D>(PlotNames[iPlot]+"_response")); // response matrix
        TH2D* TotalCovHist = (TH2D*)(CovFile->Get<TH2D>(PlotNames[iPlot])); // total cov matrix
        TH1D* TruePlot = (TH1D*)(MatrixFile->Get<TH1D>(PlotNames[iPlot]+"_true")); // all true generated events
        TH1D* RecoPlot = (TH1D*)(SelectionFile->Get<TH1D>(PlotNames[iPlot]+"_reco")); // reco events
        TH1D* BkgPlot = (TH1D*)(SelectionFile->Get<TH1D>(PlotNames[iPlot]+"_bkg")); // bkg events
        RecoPlot->Add(BkgPlot, -1); // subtract background from reco events

        // Get dimensions for matrices
        int n = TruePlot->GetNbinsX();
        int m = RecoPlot->GetNbinsX();
        double edges[n+1];
        for (int i = 0; i < n+1; i++) { edges[i] = TruePlot->GetBinLowEdge(i+1); }

        // Create objects to store matrices/vectors from Wiener SVD
        TMatrixD AddSmear(n,n);
        TVectorD WF(n);
        TMatrixD UnfoldCov(n,n);
        TMatrixD CovRotation(n,n);

        // Scale histograms
        TruePlot->Scale(Units / (IntegratedFlux * NTargets));
        RecoPlot->Scale(Units / (IntegratedFlux * NTargets));

        // Convert histograms to matrices/vectors
        TVectorD SignalVector(n); H2V(TruePlot, SignalVector);
        TVectorD MeasureVector(m); H2V(RecoPlot, MeasureVector);
        TMatrixD ResponseMatrix(m, n); H2M(ResponseHist, ResponseMatrix, kFALSE);
        TMatrixD CovarianceMatrix(m, m); H2M(TotalCovHist, CovarianceMatrix, kTRUE);

        TVectorD unfold = WienerSVD(
            ResponseMatrix,
            SignalVector,
            MeasureVector,
            CovarianceMatrix,
            2,
            0.5,
            AddSmear,
            WF,
            UnfoldCov,
            CovRotation
        );
        // Reweight unfolded covariance matrix
        TH2D* UnfTotalCovHisto = new TH2D("UnfTotalCov"+PlotNames[iPlot],"UnfTotalCov" + PlotNames[iPlot],n, edges, n, edges);
        M2H(UnfoldCov, UnfTotalCovHisto); tools.Reweight2D(UnfTotalCovHisto);
        SaveFile->WriteObject(UnfTotalCovHisto, PlotNames[iPlot]+"_unf_cov");

        // Get transpose cov rotation matrix
        TMatrixD CovRotationT (TMatrixD::kTransposed, CovRotation);

        // Get unfolded cross-section
        TH1D* UnfoldedSpectrum = new TH1D("Unfolded"+PlotNames[iPlot],";"+(TString)VarLabels[iPlot]+";"+(TString)YLabels[iPlot], n, edges);
        V2H(unfold, UnfoldedSpectrum); tools.Reweight(UnfoldedSpectrum);
        SaveFile->WriteObject(UnfoldedSpectrum, PlotNames[iPlot]+"_unf_spectrum");

        // Add smear to signal
        TH1D* SmearedSignal = new TH1D("SmearedTrue"+PlotNames[iPlot],";"+(TString)VarLabels[iPlot]+";"+(TString)YLabels[iPlot], n, edges);
        TVectorD SmearedVector = AddSmear * SignalVector;
        V2H(SmearedVector, SmearedSignal); tools.Reweight(SmearedSignal);

        // Declare canvas
        TCanvas* PlotCanvas = new TCanvas(PlotNames[iPlot],PlotNames[iPlot],205,34,1124,768);

        //////////////////
        // Smearing matrix
        //////////////////

        TH2D* SmearMatrixHisto = new TH2D("Smearing"+PlotNames[iPlot], "Smearing"+PlotNames[iPlot], n, edges, n, edges);
        M2H(AddSmear, SmearMatrixHisto);

        // Margins for matrix
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        double SmearMin = SmearMatrixHisto->GetMinimum();
        double SmearMax = SmearMatrixHisto->GetMaximum();
        SmearMatrixHisto->GetZaxis()->SetRangeUser(SmearMin,SmearMax);
        SmearMatrixHisto->GetZaxis()->CenterTitle();
        SmearMatrixHisto->GetZaxis()->SetTitleFont(FontStyle);
        SmearMatrixHisto->GetZaxis()->SetTitleSize(TextSize);
        SmearMatrixHisto->GetZaxis()->SetLabelFont(FontStyle);
        SmearMatrixHisto->GetZaxis()->SetLabelSize(TextSize);
        SmearMatrixHisto->GetZaxis()->SetNdivisions(5);

        SmearMatrixHisto->GetXaxis()->SetTitle("True " + (TString)VarLabels.at(iPlot));
        SmearMatrixHisto->GetXaxis()->CenterTitle();
        SmearMatrixHisto->GetXaxis()->SetTitleOffset(1.1);
        SmearMatrixHisto->GetXaxis()->SetTitleFont(FontStyle);
        SmearMatrixHisto->GetXaxis()->SetTitleSize(TextSize);
        SmearMatrixHisto->GetXaxis()->SetLabelFont(FontStyle);
        SmearMatrixHisto->GetXaxis()->SetLabelSize(TextSize);
        SmearMatrixHisto->GetXaxis()->SetNdivisions(5);

        SmearMatrixHisto->GetYaxis()->SetTitle("Reco " + (TString)VarLabels.at(iPlot));
        SmearMatrixHisto->GetYaxis()->CenterTitle();
        SmearMatrixHisto->GetYaxis()->SetTitleOffset(1.1);
        SmearMatrixHisto->GetYaxis()->SetTitleFont(FontStyle);
        SmearMatrixHisto->GetYaxis()->SetTitleSize(TextSize);
        SmearMatrixHisto->GetYaxis()->SetLabelFont(FontStyle);
        SmearMatrixHisto->GetYaxis()->SetLabelSize(TextSize);
        SmearMatrixHisto->GetYaxis()->SetNdivisions(5);

        PlotCanvas->cd();
        SmearMatrixHisto->Draw("colz");
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Smear/"+PlotNames[iPlot]+".pdf");
        SaveFile->WriteObject(SmearMatrixHisto, PlotNames[iPlot]+"_smear");

        ///////////////////////////
        // Bin by bin uncertainties
        ///////////////////////////

        // Get xsec cov matrix
        TMatrixD XSecCov(n, n);
        for (int iSyst = 0; iSyst < (int) XSecSystsVector.size(); iSyst++) {
            TH2D* CovHist = (TH2D*)(CovFiles[iSyst]->Get<TH2D>(PlotNames[iPlot]+"_cov"));
            TMatrixD CovMatrix(n, n); H2M(CovHist, CovMatrix, kTRUE);
            XSecCov += CovMatrix;
        }
        TMatrixD UnfXSecCov = CovRotation * XSecCov * CovRotationT;
        TH2D* UnfXSecCovHist = new TH2D("UnfCovXSec"+PlotNames[iPlot], "UnfCovXSec"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfXSecFracCovHist = new TH2D("UnfFracCovXSec"+PlotNames[iPlot], "UnfFracCovXSec"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfXSecCorrHist = new TH2D("UnfCorrXSec"+PlotNames[iPlot], "UnfCorrXSec"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfXSecCov, UnfXSecCovHist); tools.Reweight2D(UnfXSecCovHist);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfXSecCovHist, UnfXSecFracCovHist, UnfXSecCorrHist, n);
        TMatrixD UnfXSecFracCov(n, n); H2M(UnfXSecFracCovHist, UnfXSecFracCov, kTRUE);

        // Get flux cov matrix
        TMatrixD FluxCov(n, n);
        for (int iSyst = 0; iSyst < (int) FluxSystsVector.size(); iSyst++) {
            TH2D* CovHist = (TH2D*)(CovFiles[iSyst + (int)XSecSystsVector.size()]->Get<TH2D>(PlotNames[iPlot]+"_cov"));
            TMatrixD CovMatrix(n, n); H2M(CovHist, CovMatrix, kTRUE);
            FluxCov += CovMatrix;
        }
        TMatrixD UnfFluxCov = CovRotation * FluxCov * CovRotationT;
        TH2D* UnfFluxCovHist = new TH2D("UnfCovFlux"+PlotNames[iPlot], "UnfCovFlux"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfFluxFracCovHist = new TH2D("UnfFracCovFlux"+PlotNames[iPlot], "UnfFracCovFlux"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfFluxCorrHist = new TH2D("UnfCorrFlux"+PlotNames[iPlot], "UnfCorrFlux"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfFluxCov, UnfFluxCovHist); tools.Reweight2D(UnfFluxCovHist);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfFluxCovHist, UnfFluxFracCovHist, UnfFluxCorrHist, n);
        TMatrixD UnfFluxFracCov(n, n); H2M(UnfFluxFracCovHist, UnfFluxFracCov, kTRUE);

        int offset = XSecSystsVector.size() + FluxSystsVector.size();

        // Get stat cov matrix
        TH2D* StatCovHist = (TH2D*)(CovFiles[offset]->Get<TH2D>(PlotNames[iPlot]+"_cov"));
        TMatrixD StatCov(n, n); H2M(StatCovHist, StatCov, kTRUE); TMatrixD UnfStatCov = CovRotation * StatCov * CovRotationT;
        TH2D* UnfStatCovHist = new TH2D("UnfCovStat"+PlotNames[iPlot], "UnfCovStat"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfStatFracCovHist = new TH2D("UnfFracCovStat"+PlotNames[iPlot], "UnfFracCovStat"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfStatCorrHist = new TH2D("UnfCorrStat"+PlotNames[iPlot], "UnfCorrStat"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfStatCov, UnfStatCovHist); tools.Reweight2D(UnfStatCovHist);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfStatCovHist, UnfStatFracCovHist, UnfStatCorrHist, n);
        TMatrixD UnfStatFracCov(n, n); H2M(UnfStatFracCovHist, UnfStatFracCov, kTRUE);

        // Get POT cov matrix
        TH2D* POTCovHist = (TH2D*)(CovFiles[offset + 1]->Get<TH2D>(PlotNames[iPlot]+"_cov"));
        TMatrixD POTCov(n, n); H2M(POTCovHist, POTCov, kTRUE); TMatrixD UnfPOTCov = CovRotation * POTCov * CovRotationT;
        TH2D* UnfPOTCovHist = new TH2D("UnfCovPOT"+PlotNames[iPlot], "UnfCovPOT"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfPOTFracCovHist = new TH2D("UnfFracCovPOT"+PlotNames[iPlot], "UnfFracCovPOT"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfPOTCorrHist = new TH2D("UnfCorrPOT"+PlotNames[iPlot], "UnfCorrPOT"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfPOTCov, UnfPOTCovHist); tools.Reweight2D(UnfPOTCovHist);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfPOTCovHist, UnfPOTFracCovHist, UnfPOTCorrHist, n);
        TMatrixD UnfPOTFracCov(n, n); H2M(UnfPOTFracCovHist, UnfPOTFracCov, kTRUE);

        // Get NTargets cov matrix
        TH2D* NTargetsCovHist = (TH2D*)(CovFiles[offset + 2]->Get<TH2D>(PlotNames[iPlot]+"_cov"));
        TMatrixD NTargetsCov(n, n); H2M(NTargetsCovHist, NTargetsCov, kTRUE); TMatrixD UnfNTargetsCov = CovRotation * NTargetsCov * CovRotationT;
        TH2D* UnfNTargetsCovHist = new TH2D("UnfCovNTargets"+PlotNames[iPlot], "UnfCovNTargets"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfNTargetsFracCovHist = new TH2D("UnfFracCovNTargets"+PlotNames[iPlot], "UnfFracCovNTargets"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfNTargetsCorrHist = new TH2D("UnfCorrNTargets"+PlotNames[iPlot], "UnfCorrNTargets"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfNTargetsCov, UnfNTargetsCovHist); tools.Reweight2D(UnfNTargetsCovHist);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfNTargetsCovHist, UnfNTargetsFracCovHist,UnfNTargetsCorrHist, n);
        TMatrixD UnfNTargetsFracCov(n, n); H2M(UnfNTargetsFracCovHist, UnfNTargetsFracCov, kTRUE);

        // Get Detector cov matrix
        TH2D* DetectorCovHist = (TH2D*)(CovFiles[offset + 3]->Get<TH2D>(PlotNames[iPlot]+"_cov"));
        TMatrixD DetectorCov(n, n); H2M(DetectorCovHist, DetectorCov, kTRUE); TMatrixD UnfDetectorCov = CovRotation * DetectorCov * CovRotationT;
        TH2D* UnfDetectorCovHist = new TH2D("UnfCovDetector"+PlotNames[iPlot], "UnfCovDetector"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfDetectorFracCovHist = new TH2D("UnfFracCovDetector"+PlotNames[iPlot], "UnfFracCovDetector"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfDetectorCorrHist = new TH2D("UnfCorrDetector"+PlotNames[iPlot], "UnfCorrDetector"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfDetectorCov, UnfDetectorCovHist); tools.Reweight2D(UnfDetectorCovHist);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfDetectorCovHist, UnfDetectorFracCovHist, UnfDetectorCorrHist, n);
        TMatrixD UnfDetectorFracCov(n, n); H2M(UnfDetectorFracCovHist, UnfDetectorFracCov, kTRUE);

        // Get Reinteraction cov matrix
        TH2D* ReinteractionCovHist = (TH2D*)(CovFiles[offset + 4]->Get<TH2D>(PlotNames[iPlot]+"_cov"));
        TMatrixD ReinteractionCov(n, n); H2M(ReinteractionCovHist, ReinteractionCov, kTRUE); TMatrixD UnfReinteractionCov = CovRotation * ReinteractionCov * CovRotationT;
        TH2D* UnfReinteractionCovHist = new TH2D("UnfCovReinteraction"+PlotNames[iPlot], "UnfCovReinteraction"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfReinteractionFracCovHist = new TH2D("UnfFracCovReinteraction"+PlotNames[iPlot], "UnfFracCovReinteraction"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfReinteractionCorrHist = new TH2D("UnfCorrReinteraction"+PlotNames[iPlot], "UnfCorrReinteraction"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfReinteractionCov, UnfReinteractionCovHist); tools.Reweight2D(UnfReinteractionCovHist);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfReinteractionCovHist, UnfReinteractionFracCovHist, UnfReinteractionCorrHist, n);
        TMatrixD UnfReinteractionFracCov(n, n); H2M(UnfReinteractionFracCovHist, UnfReinteractionFracCov, kTRUE);

        // Get MCStat cov matrix
        TH2D* MCStatCovHist = (TH2D*)(CovFiles[offset + 5]->Get<TH2D>(PlotNames[iPlot]+"_cov"));
        TMatrixD MCStatCov(n, n); H2M(MCStatCovHist, MCStatCov, kTRUE); TMatrixD UnfMCStatCov = CovRotation * MCStatCov * CovRotationT;
        TH2D* UnfMCStatCovHist = new TH2D("UnfCovMCStat"+PlotNames[iPlot], "UnfCovMCStat"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfMCStatFracCovHist = new TH2D("UnfFracCovMCStat"+PlotNames[iPlot], "UnfFracCovMCStat"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfMCStatCorrHist = new TH2D("UnfCorrMCStat"+PlotNames[iPlot], "UnfCorrMCStat"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfMCStatCov, UnfMCStatCovHist); tools.Reweight2D(UnfMCStatCovHist);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfMCStatCovHist, UnfMCStatFracCovHist, UnfMCStatCorrHist, n);
        TMatrixD UnfMCStatFracCov(n, n); H2M(UnfMCStatFracCovHist, UnfMCStatFracCov, kTRUE);

        // Histograms for uncertainties
        TH1D* XSecHisto = new TH1D("XSec"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);
        TH1D* FluxHisto = new TH1D("Flux"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);
        TH1D* StatsHisto = new TH1D("Stats"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);
        TH1D* POTHisto = new TH1D("POT"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);
        TH1D* NTargetsHisto = new TH1D("NTargets"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);
        TH1D* DetectorHisto = new TH1D("Detector"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);
        TH1D* ReinteractionHisto = new TH1D("Reinteraction"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);
        TH1D* MCStatHisto = new TH1D("MCStat"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);
        TH1D* TotalHisto = new TH1D("Total"+PlotNames[iPlot], ";;Uncertainty [%]", n, edges);

        // Histograms to store total frac cov and corr
        TH2D* UnfTotalFracCovHisto = new TH2D("UnfTotalFracCov"+PlotNames[iPlot],"UnfTotalFracCov" + PlotNames[iPlot],n, edges, n, edges);
        TH2D* UnfTotalCorrHisto = new TH2D("UnfTotalCorr"+PlotNames[iPlot],"UnfTotalCorr" + PlotNames[iPlot],n, edges, n, edges);
        SelectionHelpers::GetFracCovAndCorrMatrix(UnfoldedSpectrum, UnfTotalCovHisto, UnfTotalFracCovHisto, UnfTotalCorrHisto, n);
        TMatrixD UnfTotalFracCov(n, n); H2M(UnfTotalFracCovHisto, UnfTotalFracCov, kTRUE);

        // Loop over each bin
        for (int iBin = 0; iBin < n; iBin++) {
            double Total = 0.;

            // Add xsec uncertainties
            XSecHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfXSecFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(UnfXSecFracCov(iBin, iBin)) * 100, 2);

            // Add flux uncertainties
            FluxHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfFluxFracCov(iBin, iBin)) * 100);
            Total +=  TMath::Power(TMath::Sqrt(UnfFluxFracCov(iBin, iBin)) * 100, 2);

            // Add Stat systematics
            StatsHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfStatFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(UnfStatFracCov(iBin, iBin)) * 100, 2);

            // Add POT systematics
            POTHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfPOTFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(UnfPOTFracCov(iBin, iBin)) * 100, 2);

            // Add NTargets systematics
            NTargetsHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfNTargetsFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(UnfNTargetsFracCov(iBin, iBin)) * 100, 2);

            // Add Detector systematics
            DetectorHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfDetectorFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(UnfDetectorFracCov(iBin, iBin)) * 100, 2);

            // Add Reinteraction systematics
            ReinteractionHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfReinteractionFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(UnfReinteractionFracCov(iBin, iBin)) * 100, 2);

            // Add MCStat systematics
            MCStatHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfMCStatFracCov(iBin, iBin)) * 100);
            Total += TMath::Power(TMath::Sqrt(UnfMCStatFracCov(iBin, iBin)) * 100, 2);

            // Total 
            TotalHisto->SetBinContent(iBin + 1, TMath::Sqrt(UnfTotalFracCov(iBin, iBin)) * 100);

            // These should be the same
            double frac_unc = TMath::Sqrt(UnfTotalFracCov(iBin, iBin)) * 100;
	    double sum_unc = TMath::Sqrt(Total);

	    if ( frac_unc != sum_unc ) {

            	std::cout << frac_unc << " ==? ";
            	std::cout << sum_unc << std::endl;;

	    }

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
        if (PlotNames[iPlot] == "EventCount") {
            XSecHisto->GetXaxis()->SetLabelSize(0);
            XSecHisto->GetXaxis()->SetTitleSize(0);
        } else {
            XSecHisto->GetXaxis()->SetTitleFont(FontStyle);
            XSecHisto->GetXaxis()->SetLabelFont(FontStyle);
            XSecHisto->GetXaxis()->SetLabelSize(TextSize);
            XSecHisto->GetXaxis()->SetTitleSize(TextSize);
            XSecHisto->GetXaxis()->SetTitleOffset(1.);
            XSecHisto->GetXaxis()->CenterTitle();
            XSecHisto->GetXaxis()->SetTitle(UnfoldedSpectrum->GetXaxis()->GetTitle());
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
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/UnfBinUncertainties/"+PlotNames[iPlot]+".pdf");

        leg->Clear();

        //////////////////////////
        // Unfolded cross-sections
        //////////////////////////

        // Separate unfolded covariance matrices into: norm, shape, and stat
        TMatrixD ReweightedUnfTotalCov(n, n); H2M(UnfTotalCovHisto, ReweightedUnfTotalCov, kTRUE);
        TMatrixD ReweightedUnfStatCov(n, n); H2M(UnfStatCovHist, ReweightedUnfStatCov, kTRUE);
        TMatrixD UnfoldCovNoStat = ReweightedUnfTotalCov - ReweightedUnfStatCov; // Remove stat from total cov matrix
        std::vector<TMatrixD> DecompOutput = tools.MatrixDecomp(n, MeasureVector, UnfoldCovNoStat); // Get norm and shape covariances
        TMatrixD UnfNormCov = DecompOutput[0]; TMatrixD UnfShapeCov = DecompOutput[1];

        TH2D* UnfNormCovHist = new TH2D("UnfCovNorm"+PlotNames[iPlot], "UnfCovNorm"+PlotNames[iPlot], n, edges, n, edges);
        TH2D* UnfShapeCovHist = new TH2D("UnfCovShape"+PlotNames[iPlot], "UnfCovShape"+PlotNames[iPlot], n, edges, n, edges);
        M2H(UnfNormCov, UnfNormCovHist); M2H(UnfShapeCov, UnfShapeCovHist);

        // Margins for unfolded xsecs
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        // Deserialize double differential plots
        if (PlotNames[iPlot].Contains("Serial")) {
            auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator["True"+PlotNames[iPlot]+"Plot"];
            auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
            int StartIndex = 0; int MatrixIndex = 0;

            // Loop over slices
            for (int iSlice = 0; iSlice < NSlices; iSlice++) {
                // Slice name
                TString SlicePlotName = PlotNames[iPlot] + "_" + TString(std::to_string(iSlice));

                // Get slice width
                double SliceWidth = SliceDiscriminators[iSlice + 1] - SliceDiscriminators[iSlice]; 

                // Get number of bins
                int SliceNBins = SerialVectorBins.at(iSlice);
                std::vector<double> SerialSliceBinning;

                for (int iBin = 0; iBin < SliceNBins + 1; iBin++) {
                    double value = SerialVectorRanges.at(StartIndex + iBin);
                    SerialSliceBinning.push_back(value);
                } // End of the number of bins and the bin ranges declaration

                // Slice true and reco true histos
                TH1D* SlicedSmearedSignal = tools.GetHistoBins(
                    SmearedSignal,
                    SerialVectorLowBin.at(iSlice),
                    SerialVectorHighBin.at(iSlice),
                    SliceWidth,
                    SerialSliceBinning,
                    "SmearedSignal"
                );
                TH1D* SlicedUnfoldedSpectrum = tools.GetHistoBins(
                    UnfoldedSpectrum,
                    SerialVectorLowBin.at(iSlice),
                    SerialVectorHighBin.at(iSlice),
                    SliceWidth,
                    SerialSliceBinning,
                    "UnfoldedSpectrum"
                );

                // Create error band
                TGraphAsymmErrors* StatErrorBand = new TGraphAsymmErrors;
                TGraphAsymmErrors* ShapeErrorBand = new TGraphAsymmErrors;
                TH1D* NormErrorHisto = new TH1D("Norm"+SlicePlotName, "", SerialSliceBinning.size() - 1, SerialSliceBinning.data());

                for (int iBin = 1; iBin < SliceNBins + 1; ++iBin) {
                    const double xnom = SlicedUnfoldedSpectrum->GetXaxis()->GetBinCenter(iBin);
                    const double ynom = SlicedUnfoldedSpectrum->GetBinContent(iBin);
                    const double dx = SlicedUnfoldedSpectrum->GetXaxis()->GetBinWidth(iBin);

                    StatErrorBand->SetPoint(iBin, xnom, ynom);
                    ShapeErrorBand->SetPoint(iBin, xnom, ynom);

                    StatErrorBand->SetPointError(
                        iBin, 0, 0,
                        TMath::Sqrt(UnfStatCovHist->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx),
                        TMath::Sqrt(UnfStatCovHist->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx)
                    );
                    ShapeErrorBand->SetPointError(
                        iBin, 0, 0,
                        TMath::Sqrt(UnfShapeCovHist->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx) + 
                        TMath::Sqrt(UnfStatCovHist->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx),
                        TMath::Sqrt(UnfShapeCovHist->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx) + 
                        TMath::Sqrt(UnfStatCovHist->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx)
                    );
                    NormErrorHisto->SetBinContent(iBin, TMath::Sqrt(UnfNormCovHist->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx));
                }

                // Style 
                SlicedSmearedSignal->GetXaxis()->SetNdivisions(5);
                if (PlotNames[iPlot] == "EventCount") {
                    SlicedSmearedSignal->GetXaxis()->SetLabelSize(0);
                    SlicedSmearedSignal->GetXaxis()->SetTitleSize(0);
                } else {
                    SlicedSmearedSignal->GetXaxis()->SetTitleFont(FontStyle);
                    SlicedSmearedSignal->GetXaxis()->SetLabelFont(FontStyle);
                    SlicedSmearedSignal->GetXaxis()->SetLabelSize(TextSize);
                    SlicedSmearedSignal->GetXaxis()->SetTitleSize(TextSize);
                    SlicedSmearedSignal->GetXaxis()->SetTitleOffset(1.);
                    SlicedSmearedSignal->GetXaxis()->CenterTitle();

                    std::string VarLabel = (std::string) VarLabels.at(iPlot);
                    VarLabel.erase(VarLabel.end() - 7, VarLabel.end());
                    SlicedSmearedSignal->GetXaxis()->SetTitle(VarLabel.c_str());
                }
                SlicedSmearedSignal->GetYaxis()->SetTitleFont(FontStyle);
                SlicedSmearedSignal->GetYaxis()->SetLabelFont(FontStyle);
                SlicedSmearedSignal->GetYaxis()->SetLabelSize(TextSize);
                SlicedSmearedSignal->GetYaxis()->SetTitleSize(TextSize);
                SlicedSmearedSignal->GetYaxis()->SetNdivisions(6);
                SlicedSmearedSignal->GetYaxis()->SetTitleOffset(1.);
                SlicedSmearedSignal->GetYaxis()->SetTickSize(0);
                SlicedSmearedSignal->GetYaxis()->CenterTitle();

                // Create legend object
                TLegend* leg = new TLegend(0.2,0.73,0.55,0.83);
                leg->SetBorderSize(0);
                leg->SetNColumns(2);
                leg->SetTextSize(TextSize*0.8);
                leg->SetTextFont(FontStyle);

                TLegendEntry* legSlicedSmear = leg->AddEntry(SlicedSmearedSignal,"True","l");
                SlicedSmearedSignal->SetLineColor(kRed+1);
                SlicedSmearedSignal->SetLineWidth(4);

                TLegendEntry* legSlicedUnf = leg->AddEntry(SlicedUnfoldedSpectrum,"Unfolded","p");
                SlicedUnfoldedSpectrum->SetLineColor(kBlack);
                SlicedUnfoldedSpectrum->SetMarkerColor(kBlack);
                SlicedUnfoldedSpectrum->SetMarkerStyle(20);
                SlicedUnfoldedSpectrum->SetMarkerSize(0.5);

                TLegendEntry* legNorm = leg->AddEntry(NormErrorHisto,"Norm","f");
                NormErrorHisto->SetLineColor(kGray);
                NormErrorHisto->SetFillColorAlpha(kGray, 0.5);
                NormErrorHisto->SetFillStyle(1001);

                TLegendEntry* legStatShape = leg->AddEntry(ShapeErrorBand,"Stat#oplusShape","ep");	
                ShapeErrorBand->SetLineWidth(2);

                double imax = TMath::Max(SlicedUnfoldedSpectrum->GetMaximum(),SlicedSmearedSignal->GetMaximum());
                SlicedUnfoldedSpectrum->GetYaxis()->SetRangeUser(0.,1.35*imax);
                SlicedSmearedSignal->GetYaxis()->SetRangeUser(0.,1.35*imax);

                PlotCanvas->cd();
                SlicedSmearedSignal->Draw("hist");
                SlicedUnfoldedSpectrum->Draw("p0 hist same");
                StatErrorBand->Draw("e1 same");
                ShapeErrorBand->Draw("e1 same");
                NormErrorHisto->Draw("hist same");
                leg->Draw();

                // Slice label
                TLatex *textSlice = new TLatex();
                TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel["True"+PlotNames[iPlot]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
                textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

                // Save histogram
                PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Unfolded/"+SlicePlotName+".pdf");	

                // Save error band
                SaveFile->WriteObject(StatErrorBand, SlicePlotName+"_stat_band");
                SaveFile->WriteObject(ShapeErrorBand, SlicePlotName+"_shape_band");
                SaveFile->WriteObject(NormErrorHisto, SlicePlotName+"_norm_histo");

                StartIndex += (SliceNBins + 1); MatrixIndex += SliceNBins;
            }
        } else {
            // Create error band
            TGraphAsymmErrors* StatErrorBand = new TGraphAsymmErrors;
            TGraphAsymmErrors* ShapeErrorBand = new TGraphAsymmErrors;
            TH1D* NormErrorHisto = new TH1D("Norm"+PlotNames[iPlot], "", n, edges);

            for (int iBin = 1; iBin < UnfoldedSpectrum->GetNbinsX() + 1; ++iBin) {
                const double xnom = UnfoldedSpectrum->GetXaxis()->GetBinCenter(iBin);
                const double ynom = UnfoldedSpectrum->GetBinContent(iBin);

                StatErrorBand->SetPoint(iBin, xnom, ynom);
                ShapeErrorBand->SetPoint(iBin, xnom, ynom);

                StatErrorBand->SetPointError(
                    iBin, 0, 0,
                    TMath::Sqrt(UnfStatCovHist->GetBinContent(iBin, iBin)),
                    TMath::Sqrt(UnfStatCovHist->GetBinContent(iBin, iBin))
                );
                ShapeErrorBand->SetPointError(
                    iBin, 0, 0,
                    TMath::Sqrt(UnfShapeCovHist->GetBinContent(iBin, iBin)) +
                    TMath::Sqrt(UnfStatCovHist->GetBinContent(iBin, iBin)),
                    TMath::Sqrt(UnfShapeCovHist->GetBinContent(iBin, iBin)) + 
                    TMath::Sqrt(UnfStatCovHist->GetBinContent(iBin, iBin))
                );
                NormErrorHisto->SetBinContent(iBin, TMath::Sqrt(UnfNormCovHist->GetBinContent(iBin, iBin)));
            }

            // Style 
            SmearedSignal->GetXaxis()->SetNdivisions(5);
            if (PlotNames[iPlot] == "EventCount") {
                SmearedSignal->GetXaxis()->SetLabelSize(0);
                SmearedSignal->GetXaxis()->SetTitleSize(0);
            } else {
                SmearedSignal->GetXaxis()->SetTitleFont(FontStyle);
                SmearedSignal->GetXaxis()->SetLabelFont(FontStyle);
                SmearedSignal->GetXaxis()->SetLabelSize(TextSize);
                SmearedSignal->GetXaxis()->SetTitleSize(TextSize);
                SmearedSignal->GetXaxis()->SetTitleOffset(1.);
                SmearedSignal->GetXaxis()->CenterTitle();
                SmearedSignal->GetXaxis()->SetTitle(VarLabels.at(iPlot).c_str());
            }
            SmearedSignal->GetYaxis()->SetTitleFont(FontStyle);
            SmearedSignal->GetYaxis()->SetLabelFont(FontStyle);
            SmearedSignal->GetYaxis()->SetLabelSize(TextSize);
            SmearedSignal->GetYaxis()->SetTitleSize(TextSize);
            SmearedSignal->GetYaxis()->SetNdivisions(6);
            SmearedSignal->GetYaxis()->SetTitleOffset(1.);
            SmearedSignal->GetYaxis()->SetTickSize(0);
            SmearedSignal->GetYaxis()->CenterTitle();

            TLegend* leg = new TLegend(0.2,0.73,0.55,0.83);
            leg->SetBorderSize(0);
            leg->SetNColumns(2);
            leg->SetTextSize(TextSize*0.8);
            leg->SetTextFont(FontStyle);

            TLegendEntry* legSmeared = leg->AddEntry(SmearedSignal,"True","l");
            SmearedSignal->SetLineColor(kRed+1);
            SmearedSignal->SetLineWidth(4);

            TLegendEntry* legUnfSpectrum = leg->AddEntry(UnfoldedSpectrum,"Unfolded","p");
            UnfoldedSpectrum->SetLineColor(kBlack);
            UnfoldedSpectrum->SetLineWidth(4);
            UnfoldedSpectrum->SetMarkerColor(kBlack);
			UnfoldedSpectrum->SetMarkerStyle(20);
			UnfoldedSpectrum->SetMarkerSize(0.5);

            TLegendEntry* legNorm = leg->AddEntry(NormErrorHisto,"Norm","f");
            NormErrorHisto->SetLineColor(kGray);
            NormErrorHisto->SetFillColorAlpha(kGray, 0.5);
            NormErrorHisto->SetFillStyle(1001);

            TLegendEntry* legStatShape = leg->AddEntry(ShapeErrorBand,"Stat#oplusShape","ep");	
            ShapeErrorBand->SetLineWidth(2);

            double imax = TMath::Max(UnfoldedSpectrum->GetMaximum(),SmearedSignal->GetMaximum());
            for(int i = 1; i < ShapeErrorBand->GetN() - 1; ++i){
                imax = std::max(imax, ShapeErrorBand->GetY()[i] + ShapeErrorBand->GetErrorYhigh(i));
            }
            UnfoldedSpectrum->GetYaxis()->SetRangeUser(0.,1.3 * imax);
            SmearedSignal->GetYaxis()->SetRangeUser(0.,1.3 * imax);	

            PlotCanvas->cd();
            SmearedSignal->Draw("hist");
            UnfoldedSpectrum->Draw("p0 hist same");
            StatErrorBand->Draw("e1 same");
            ShapeErrorBand->Draw("e1 same");
            NormErrorHisto->Draw("hist same");
            leg->Draw();

            // Save histogram
            PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Unfolded/"+PlotNames[iPlot]+".pdf");

            // Save error band
            SaveFile->WriteObject(StatErrorBand, PlotNames[iPlot]+"_stat_band");
            SaveFile->WriteObject(ShapeErrorBand, PlotNames[iPlot]+"_shape_band");
            SaveFile->WriteObject(NormErrorHisto, PlotNames[iPlot]+"_norm_histo");
        }
        delete PlotCanvas;
    }
}
