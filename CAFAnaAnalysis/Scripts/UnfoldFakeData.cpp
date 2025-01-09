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
#include <filesystem>

// Helpers includes.
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Tools.cxx"
#include "../../Utils/Util.C"
#include "../../Utils/WienerSVD.C"
#include "../../Utils/Constants.h"

using namespace std;
using namespace Constants;

void UnfoldFakeData() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.1f");

    Tools tools;

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/UnfoldedFake.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    // Load root file with fake data
    TString FakeRootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionFake.root";
    std::unique_ptr<TFile> FakeFile(TFile::Open(FakeRootFilePath));

    // Load root file with nominal histograms
    TString SelectionRootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> SelectionFile(TFile::Open(SelectionRootFilePath));

    // Load root file with response matrices
    TString MatrixRootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Matrix.root";
    std::unique_ptr<TFile> MatrixFile(TFile::Open(MatrixRootFilePath));

    // Vector with all cov matrices to use
    std::vector<std::unique_ptr<TFile>> CovFiles;

    // Add xsec systematics
    for (int iSyst = 0; iSyst < (int) XSecSystsVector.size(); iSyst++) {
        std::string SystName = std::get<0>(XSecSystsVector.at(iSyst));
        TString FilePath =  "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematics"+TString(SystName)+".root";
        std::unique_ptr<TFile> File(TFile::Open(FilePath));
        CovFiles.push_back(std::move(File));
    }

    // Add MCStat systematics
    std::unique_ptr<TFile> MCStatFile(TFile::Open("/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsMCStat.root"));
    CovFiles.push_back(std::move(MCStatFile));

    // Dir to save plots
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    for (std::size_t iData = 0; iData < FakeDataNames.size(); ++iData) {
        // Create directory for this fake data if it does not exist yet
        std::filesystem::create_directory((std::string)dir+"/Figs/CAFAna/FakeDataStudies/"+(std::string)FakeDataNames[iData]);

        for (std::size_t iVar = 0; iVar < PlotNames.size(); ++iVar) {
            // Load nominal plots
            TH2D* ResponseHist = (TH2D*)(MatrixFile->Get<TH2D>(PlotNames[iVar]+"_response"));
            TH1D* TruePlot = (TH1D*)(MatrixFile->Get<TH1D>(PlotNames[iVar]+"_true"));
            TH1D* BkgPlot = (TH1D*)(SelectionFile->Get<TH1D>(PlotNames[iVar]+"_bkg"));
            
            // Get fake data reco histogram
            TH1D* FakeRecoHisto = (TH1D*)(FakeFile->Get<TH1D>(PlotNames[iVar] + FakeDataNames[iData] + "_reco"));
            TH1D* FakeTrueHisto = (TH1D*)(FakeFile->Get<TH1D>(PlotNames[iVar] + FakeDataNames[iData] + "_true"));

            // Subtract MC background from fake data
            FakeRecoHisto->Add(BkgPlot, -1);

            // Scale plots
            TruePlot->Scale(Units / (IntegratedFlux * NTargets));
            FakeRecoHisto->Scale(Units / (IntegratedFlux * NTargets));
            FakeTrueHisto->Scale(Units / (IntegratedFlux * NTargets));

            int n = TruePlot->GetXaxis()->GetNbins();
            int m = FakeRecoHisto->GetNbinsX();
            double edges[n+1];
            for (int i = 0; i < n+1; i++) { edges[i] = FakeRecoHisto->GetBinLowEdge(i+1); }

            // Get our total cov matrix
            TMatrixD TotalCov(n, n);
            for (int iSyst = 0; iSyst < (int) XSecSystsVector.size(); iSyst++) {
                TH2D* CovHist = (TH2D*)(CovFiles[iSyst]->Get<TH2D>(PlotNames[iVar]+"_cov"));
                TMatrixD CovMatrix(n, n); H2M(CovHist, CovMatrix, kTRUE);
                TotalCov += CovMatrix;
            }
            int offset = XSecSystsVector.size();

            TH2D* MCStatCovHist = (TH2D*)(CovFiles[offset]->Get<TH2D>(PlotNames[iVar]+"_cov"));
            TMatrixD MCStatCov(n, n); H2M(MCStatCovHist, MCStatCov, kTRUE); TotalCov += MCStatCov;

            TH2D* StatCovHist = (TH2D*)(FakeFile->Get<TH2D>(PlotNames[iVar] + FakeDataNames[iData] + "_stat_cov"));
            TMatrixD StatCov(n, n); H2M(StatCovHist, StatCov, kTRUE); TotalCov += StatCov;

            /////////
            // Unfold
            /////////

            // Matrices/vectors to store output of unfolding
            TMatrixD AddSmear(n,n);
            TVectorD WF(n);
            TMatrixD UnfoldCov(n,n);
            TMatrixD CovRotation(n,n);

            // Convert histograms to matrices/vectors
            TVectorD SignalVector(n); H2V(TruePlot, SignalVector);
            TVectorD MeasureVector(m); H2V(FakeRecoHisto, MeasureVector);
            TMatrixD ResponseMatrix(m, n); H2M(ResponseHist, ResponseMatrix, kFALSE);

            TVectorD unfold = WienerSVD(
                ResponseMatrix,
                SignalVector,
                MeasureVector,
                TotalCov,
                2,
                0.5,
                AddSmear,
                WF,
                UnfoldCov,
                CovRotation
            );

            // Reweight unfolded covariance matrix
            TH2D* UnfTotalCovHisto = new TH2D("UnfTotalCov"+FakeDataNames[iData]+PlotNames[iVar],"UnfTotalCov" + PlotNames[iVar],n, edges, n, edges);
            M2H(UnfoldCov, UnfTotalCovHisto); tools.Reweight2D(UnfTotalCovHisto);

            // Get fake data unfolded cross-section
            TH1D* UnfoldedSpectrum = new TH1D("Unf"+FakeDataNames[iData]+PlotNames[iVar],";"+(TString)VarLabels[iVar]+";"+(TString)YLabels[iVar], n, edges);
            V2H(unfold, UnfoldedSpectrum); tools.Reweight(UnfoldedSpectrum);
            
            // Add smear to fake true signal
            TH1D* SmearedFakeSignal = new TH1D("Smeared"+FakeDataNames[iData]+PlotNames[iVar],";"+(TString)VarLabels[iVar]+";"+(TString)YLabels[iVar], n, edges);
            TVectorD FakeSignalVector(n); H2V(FakeTrueHisto, FakeSignalVector);
            TVectorD SmearedFakeVector = AddSmear * FakeSignalVector;
            V2H(SmearedFakeVector, SmearedFakeSignal); tools.Reweight(SmearedFakeSignal);
           
            // Add smear to nominal signal
            TH1D* SmearedNomSignal = new TH1D("SmearedNom"+FakeDataNames[iData]+PlotNames[iVar],";"+(TString)VarLabels[iVar]+";"+(TString)YLabels[iVar], n, edges);
            TVectorD SmearedNomVector = AddSmear * SignalVector;
            V2H(SmearedNomVector, SmearedNomSignal); tools.Reweight(SmearedNomSignal);

            ///////
            // Plot
            ///////

            TCanvas* PlotCanvas = new TCanvas(PlotNames[iVar],PlotNames[iVar],205,34,1124,768);
            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);

            if (PlotNames[iVar].Contains("Serial")) {
                auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator["True"+PlotNames[iVar]+"Plot"];
                auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
                int StartIndex = 0; int MatrixIndex = 0;

                // Get total unfolded cov matrix
                TMatrixD TotalCovMatrix(n, n); H2M(UnfTotalCovHisto, TotalCovMatrix, kTRUE);

                // Loop over slices
                for (int iSlice = 0; iSlice < NSlices; ++iSlice) {
                    ////////////////////
                    // Slice information
                    ////////////////////

                    TString SlicePlotName = PlotNames[iVar] + "_" + TString(std::to_string(iSlice));
                    double SliceWidth = SliceDiscriminators[iSlice + 1] - SliceDiscriminators[iSlice]; 
                    int SliceNBins = SerialVectorBins.at(iSlice); std::vector<double> SerialSliceBinning;
                    for (int iBin = 0; iBin < SliceNBins + 1; iBin++) {
                        double value = SerialVectorRanges.at(StartIndex + iBin);
                        SerialSliceBinning.push_back(value);
                    }

                    // Get subcov matrix
                    TMatrixD SubCovMatrix(SliceNBins, SliceNBins);
                    for (int i = MatrixIndex; i < MatrixIndex + SliceNBins; ++i) {
                        for (int j = MatrixIndex; j < MatrixIndex + SliceNBins; ++j) {
                            SubCovMatrix(i - MatrixIndex, j - MatrixIndex) = UnfoldCov(i, j);
                        }
                    }

                    // Convert sub cov matrix to histo and reweight
                    TH2D* UnfSubCovHisto = new TH2D("UnfSubCov"+FakeDataNames[iData]+SlicePlotName,"UnfSubCov"+SlicePlotName, SliceNBins, SerialSliceBinning.data(), SliceNBins, SerialSliceBinning.data());
                    M2H(SubCovMatrix, UnfSubCovHisto); 
                    UnfSubCovHisto->Scale(1 / (TMath::Power(SliceWidth, 2)));
                    tools.Reweight2D(UnfSubCovHisto);

                    ////////////////
                    // Slice histos
                    ////////////////

                    TH1D* SlicedSmearedFakeSignal = tools.GetHistoBins(
                        SmearedFakeSignal,
                        SerialVectorLowBin.at(iSlice),
                        SerialVectorHighBin.at(iSlice),
                        SliceWidth,
                        SerialSliceBinning,
                        "SmearedFake" + SlicePlotName
                    );

                    TH1D* SlicedSmearedNomSignal = tools.GetHistoBins(
                        SmearedNomSignal,
                        SerialVectorLowBin.at(iSlice),
                        SerialVectorHighBin.at(iSlice),
                        SliceWidth,
                        SerialSliceBinning,
                        "SmearedNom" + SlicePlotName
                    );
                    
                    TH1D* SlicedUnfoldedSpectrum = tools.GetHistoBins(
                        UnfoldedSpectrum,
                        SerialVectorLowBin.at(iSlice),
                        SerialVectorHighBin.at(iSlice),
                        SliceWidth,
                        SerialSliceBinning,
                        "UnfoldedSpectrum" + SlicePlotName
                    );

                    ////////////////////////////////////////////
                    // Calculate chi squared, p value, and sigma
                    ////////////////////////////////////////////

                    double ChiFake; int NDofFake; double PValFake; double SigmaFake;
                    tools.CalcChiSquared(SlicedSmearedFakeSignal, SlicedUnfoldedSpectrum, UnfSubCovHisto, ChiFake, NDofFake, PValFake, SigmaFake);

                    double ChiNom; int NDofNom; double PValNom; double SigmaNom;
                    tools.CalcChiSquared(SlicedSmearedNomSignal, SlicedUnfoldedSpectrum, UnfSubCovHisto, ChiNom, NDofNom, PValNom, SigmaNom);

                    // Create error band
                    TGraphAsymmErrors* ErrorBand = new TGraphAsymmErrors;
                    for (int iBin = 1; iBin < SliceNBins + 1; ++iBin) {
                        const double xnom = SlicedUnfoldedSpectrum->GetXaxis()->GetBinCenter(iBin);
                        const double ynom = SlicedUnfoldedSpectrum->GetBinContent(iBin);
                        ErrorBand->SetPoint(iBin, xnom, ynom);
                        const double dx = SlicedUnfoldedSpectrum->GetXaxis()->GetBinWidth(iBin);
                        ErrorBand->SetPointError(
                            iBin, 0, 0,
                            TMath::Sqrt(UnfTotalCovHisto->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx),
                            TMath::Sqrt(UnfTotalCovHisto->GetBinContent(MatrixIndex + iBin, MatrixIndex + iBin)) / (SliceWidth * dx)
                        );
                    }

                    // Style 
                    SlicedSmearedFakeSignal->GetXaxis()->SetNdivisions(5);
                    if (PlotNames[iVar] == "EventCount") {
                        SlicedSmearedFakeSignal->GetXaxis()->SetLabelSize(0);
                        SlicedSmearedFakeSignal->GetXaxis()->SetTitleSize(0);
                    } else {
                        SlicedSmearedFakeSignal->GetXaxis()->SetTitleFont(FontStyle);
                        SlicedSmearedFakeSignal->GetXaxis()->SetLabelFont(FontStyle);
                        SlicedSmearedFakeSignal->GetXaxis()->SetLabelSize(TextSize);
                        SlicedSmearedFakeSignal->GetXaxis()->SetTitleSize(TextSize);
                        SlicedSmearedFakeSignal->GetXaxis()->SetTitleOffset(1.);
                        SlicedSmearedFakeSignal->GetXaxis()->CenterTitle();
                        std::string VarLabel = (std::string) VarLabels.at(iVar);
                        VarLabel.erase(VarLabel.end() - 7, VarLabel.end()); // get rid of (bin #)
                        SlicedSmearedFakeSignal->GetXaxis()->SetTitle((TString)VarLabel + SerialNameToUnit[PlotNames[iVar]]);
                    }
                    SlicedSmearedFakeSignal->GetYaxis()->SetTitleFont(FontStyle);
                    SlicedSmearedFakeSignal->GetYaxis()->SetLabelFont(FontStyle);
                    SlicedSmearedFakeSignal->GetYaxis()->SetLabelSize(TextSize);
                    SlicedSmearedFakeSignal->GetYaxis()->SetTitleSize(TextSize);
                    SlicedSmearedFakeSignal->GetYaxis()->SetNdivisions(6);
                    SlicedSmearedFakeSignal->GetYaxis()->SetTitleOffset(1.);
                    SlicedSmearedFakeSignal->GetYaxis()->SetTickSize(0);
                    SlicedSmearedFakeSignal->GetYaxis()->CenterTitle();

                    // Create legend object
                    TLegend* leg = new TLegend(0.5,0.73,0.9,0.83);
                    leg->SetBorderSize(0);
                    leg->SetNColumns(1);
                    leg->SetTextSize(TextSize*0.8);
                    leg->SetTextFont(FontStyle);

                    TString ChiFakeLabel = "(" + 
                        tools.to_string_with_precision(ChiFake, 1) + 
                        "/" + 
                        tools.to_string_with_precision(NDofFake, 0) + 
                        ", " + 
                        tools.to_string_with_precision(PValFake, 1) + 
                        ", " + 
                        tools.to_string_with_precision(SigmaFake, 1) + "#sigma" +
                        ")";
                    TLegendEntry* legSlicedFakeSmear = leg->AddEntry(SlicedSmearedFakeSignal,"Fake true" + ChiFakeLabel,"l");
                    SlicedSmearedFakeSignal->SetLineColor(kRed+1);
                    SlicedSmearedFakeSignal->SetLineWidth(4);

                    TString ChiNomLabel = "(" + 
                        tools.to_string_with_precision(ChiNom, 1) + 
                        "/" + 
                        tools.to_string_with_precision(NDofNom, 0) + 
                        ", " + 
                        tools.to_string_with_precision(PValNom, 1) + 
                        ", " + 
                        tools.to_string_with_precision(SigmaNom, 1) + "#sigma" +
                        ")";
                    TLegendEntry* legSlicedNomSmeared = leg->AddEntry(SlicedSmearedNomSignal,"Nom true " + ChiNomLabel,"l");
                    SlicedSmearedNomSignal->SetLineColor(kBlue+8);
                    SlicedSmearedNomSignal->SetLineWidth(4);

                    TLegendEntry* legSlicedUnf = leg->AddEntry(SlicedUnfoldedSpectrum,"Unfolded","ep");
                    SlicedUnfoldedSpectrum->SetLineColor(kBlack);
                    SlicedUnfoldedSpectrum->SetMarkerColor(kBlack);
                    SlicedUnfoldedSpectrum->SetMarkerStyle(20);
                    SlicedUnfoldedSpectrum->SetMarkerSize(1.);

                    double imax = TMath::Max(SlicedUnfoldedSpectrum->GetMaximum(), SlicedSmearedFakeSignal->GetMaximum());
                    for(int i = 1; i < ErrorBand->GetN() - 1; ++i){
                        imax = std::max(imax, ErrorBand->GetY()[i] + ErrorBand->GetErrorYhigh(i));
                    }
                    SlicedSmearedFakeSignal->GetYaxis()->SetRangeUser(0., 1.35*imax);

                    PlotCanvas->cd();
                    SlicedSmearedFakeSignal->Draw("hist");
                    SlicedSmearedNomSignal->Draw("hist same");
                    SlicedUnfoldedSpectrum->Draw("e1x0 same");
                    ErrorBand->Draw("e1 same");
                    leg->Draw();

                    // Slice label
                    TLatex *textSlice = new TLatex();
                    TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel["True"+PlotNames[iVar]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
                    textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

                    // Save histogram
                    PlotCanvas->SaveAs(dir+"/Figs/CAFAna/FakeDataStudies/"+FakeDataNames[iData]+"/"+SlicePlotName+".png");

                    StartIndex += (SliceNBins + 1); MatrixIndex += SliceNBins;
                }
            } 

            ////////////////////////////////////////////
            // Calculate chi squared, p value, and sigma
            ////////////////////////////////////////////

            double ChiFake; int NDofFake; double PValFake; double SigmaFake;
            tools.CalcChiSquared(SmearedFakeSignal, UnfoldedSpectrum, UnfTotalCovHisto, ChiFake, NDofFake, PValFake, SigmaFake);

            double ChiNom; int NDofNom; double PValNom; double SigmaNom;
            tools.CalcChiSquared(SmearedNomSignal, UnfoldedSpectrum, UnfTotalCovHisto, ChiNom, NDofNom, PValNom, SigmaNom);

            // Create error band
            TGraphAsymmErrors* ErrorBand = new TGraphAsymmErrors;
            for (int iBin = 1; iBin < UnfoldedSpectrum->GetNbinsX() + 1; ++iBin) {
                const double xnom = UnfoldedSpectrum->GetXaxis()->GetBinCenter(iBin);
                const double ynom = UnfoldedSpectrum->GetBinContent(iBin);
                ErrorBand->SetPoint(iBin, xnom, ynom);
                ErrorBand->SetPointError(
                    iBin, 0, 0,
                    TMath::Sqrt(UnfTotalCovHisto->GetBinContent(iBin, iBin)),
                    TMath::Sqrt(UnfTotalCovHisto->GetBinContent(iBin, iBin))
                );
            }

            // Style 
            SmearedFakeSignal->GetXaxis()->SetNdivisions(5);
            if (PlotNames[iVar] == "EventCount") {
                SmearedFakeSignal->GetXaxis()->SetLabelSize(0);
                SmearedFakeSignal->GetXaxis()->SetTitleSize(0);
            } else {
                SmearedFakeSignal->GetXaxis()->SetTitleFont(FontStyle);
                SmearedFakeSignal->GetXaxis()->SetLabelFont(FontStyle);
                SmearedFakeSignal->GetXaxis()->SetLabelSize(TextSize);
                SmearedFakeSignal->GetXaxis()->SetTitleSize(TextSize);
                SmearedFakeSignal->GetXaxis()->SetTitleOffset(1.);
                SmearedFakeSignal->GetXaxis()->CenterTitle();
                SmearedFakeSignal->GetXaxis()->SetTitle(VarLabels.at(iVar).c_str());
            }
            SmearedFakeSignal->GetYaxis()->SetTitleFont(FontStyle);
            SmearedFakeSignal->GetYaxis()->SetLabelFont(FontStyle);
            SmearedFakeSignal->GetYaxis()->SetLabelSize(TextSize);
            SmearedFakeSignal->GetYaxis()->SetTitleSize(TextSize);
            SmearedFakeSignal->GetYaxis()->SetNdivisions(6);
            SmearedFakeSignal->GetYaxis()->SetTitleOffset(1.);
            SmearedFakeSignal->GetYaxis()->SetTickSize(0);
            SmearedFakeSignal->GetYaxis()->CenterTitle();

            TLegend* leg = new TLegend(0.5,0.73,0.9,0.83);
            leg->SetBorderSize(0);
            leg->SetNColumns(1);
            leg->SetTextSize(TextSize*0.8);
            leg->SetTextFont(FontStyle);

            TString ChiFakeLabel = "(" + 
                tools.to_string_with_precision(ChiFake, 1) + 
                "/" + 
                tools.to_string_with_precision(NDofFake, 0) + 
                ", " + 
                tools.to_string_with_precision(PValFake, 1) + 
                ", " + 
                tools.to_string_with_precision(SigmaFake, 1) + "#sigma" +
                ")";
            TLegendEntry* legFakeSmeared = leg->AddEntry(SmearedFakeSignal,"Fake true " + ChiFakeLabel,"l");
            SmearedFakeSignal->SetLineColor(kRed+1);
            SmearedFakeSignal->SetLineWidth(4);

            TString ChiNomLabel = "(" + 
                tools.to_string_with_precision(ChiNom, 1) + 
                "/" + 
                tools.to_string_with_precision(NDofNom, 0) + 
                ", " + 
                tools.to_string_with_precision(PValNom, 1) + 
                ", " + 
                tools.to_string_with_precision(SigmaNom, 1) + "#sigma" +
                ")";
            TLegendEntry* legNomSmeared = leg->AddEntry(SmearedNomSignal,"Nom true " + ChiNomLabel,"l");
            SmearedNomSignal->SetLineColor(kBlue+8);
            SmearedNomSignal->SetLineWidth(4);

            TLegendEntry* legUnfSpectrum = leg->AddEntry(UnfoldedSpectrum,"Unfolded","ep");
            UnfoldedSpectrum->SetLineColor(kBlack);
            UnfoldedSpectrum->SetMarkerColor(kBlack);
            UnfoldedSpectrum->SetMarkerStyle(20);
            UnfoldedSpectrum->SetMarkerSize(1.);

            double imax = TMath::Max(UnfoldedSpectrum->GetMaximum(), SmearedFakeSignal->GetMaximum());
            for(int i = 1; i < ErrorBand->GetN() - 1; ++i){
                imax = std::max(imax, ErrorBand->GetY()[i] + ErrorBand->GetErrorYhigh(i));
            }
            SmearedFakeSignal->GetYaxis()->SetRangeUser(0., 1.35*imax);

            PlotCanvas->cd();
            SmearedFakeSignal->Draw("hist");
            SmearedNomSignal->Draw("hist same");
            UnfoldedSpectrum->Draw("e1x0 same");
            ErrorBand->Draw("e1 same");
            leg->Draw();

            // Save histogram
            PlotCanvas->SaveAs(dir+"/Figs/CAFAna/FakeDataStudies/"+FakeDataNames[iData]+"/"+PlotNames[iVar]+".png");
            delete PlotCanvas;
        }
    }

}
