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

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionCutScan() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Directory to store figs
    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";

    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >> Spectra;

    for (int i = 0; i < 40; ++i) {
        double CutPosition = -1. + (0.05 * i);
        std::vector<double> CutArray{-1.0,CutPosition,1.0};

        const Var kVar([=](const caf::SRSliceProxy* slc) -> double {
            float fCosOpeningAngleProtons = kCosOpeningAngleProtons(slc);
            int MuonCosThetaTwoDIndex = tools.ReturnIndex(kMuonCosTheta(slc), CutArray);
            int SerialCosOpeningAngleProtonsInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                fCosOpeningAngleProtons
            );
            return SerialCosOpeningAngleProtonsInMuonCosThetaIndex;
        });
        const TruthVar kTruthVar([=](const caf::SRTrueInteractionProxy* nu) -> double {
            float fCosOpeningAngleProtons = kTruthCosOpeningAngleProtons(nu);
            int MuonCosThetaTwoDIndex = tools.ReturnIndex(kTruthMuonCosTheta(nu), CutArray);
            int SerialCosOpeningAngleProtonsInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                fCosOpeningAngleProtons
            );
            return SerialCosOpeningAngleProtonsInMuonCosThetaIndex;
        });

        auto RecoSignals = std::make_unique<Spectrum>(
            "CosOpeningAngleProtonsInMuonCosTheta",
            bCosOpeningAngleProtonsInMuonCosTheta,
            NuLoader,
            kVar,
            kNoSpillCut,
            kRecoIsSignal
        );
        auto RecoTrueSignals = std::make_unique<Spectrum>(
            "CosOpeningAngleProtonsInMuonCosTheta",
            bCosOpeningAngleProtonsInMuonCosTheta,
            NuLoader,
            kVar,
            kNoSpillCut,
            kRecoIsTrueReco
        );
        auto TrueSignals = std::make_unique<Spectrum>(
            "CosOpeningAngleProtonsInMuonCosTheta",
            bCosOpeningAngleProtonsInMuonCosTheta,
            NuLoader,
            kTruthVar,
            kTruthIsSignal,
            kNoSpillCut
        );
        Spectra.push_back({std::move(RecoSignals), std::move(RecoTrueSignals), std::move(TrueSignals)});
    }
    NuLoader.Go();

    TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
    PlotCanvas->SetTopMargin(0.13);
    PlotCanvas->SetLeftMargin(0.17);
    PlotCanvas->SetRightMargin(0.05);
    PlotCanvas->SetBottomMargin(0.16);

    TH1* EfficiencyHisto = new TH1D("Efficiency",";Cut position",40,-1,1);
    TH1* PurityHisto = new TH1D("Purity",";Cut position",40,-1,1);
    TH1* EfficiencyPurityHisto = new TH1D("EfficiencyPurity",";Cut position",40,-1,1);

    for (int i = 0; i < 40; ++i) {
        auto& [RecoSignals, RecoTrueSignals, TrueSignals] = Spectra[i];

        TH1D* RecoHisto = RecoSignals->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSignals->ToTH1(TargetPOT);
        TH1D* TrueHisto = TrueSignals->ToTH1(TargetPOT);

        // Slice the plots
        double CutPosition = -1. + (0.05 * i);
        std::vector<double> CutArray{-1.0,CutPosition,1.0};
        auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(CutArray, TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices);
        int StartIndex = 0;

        // Create vector to store deserialize plots
        std::vector<std::vector<TH1D*>> Histos;
        Histos.resize(NSlices);
        for (int iSlice = 0; iSlice < NSlices; iSlice++) {
            Histos[iSlice].resize(3);
        }

        // Loop over slices
        for (int iSlice = 0; iSlice < NSlices; iSlice++) {
            double SliceWidth = CutArray[iSlice + 1] - CutArray[iSlice];

            // Get number of bins
            int SliceNBins = SerialVectorBins.at(iSlice);
            std::vector<double> SerialSliceBinning;
            for (int iBin = 0; iBin < SliceNBins + 1; iBin++) {
                double value = SerialVectorRanges.at(StartIndex + iBin);
                SerialSliceBinning.push_back(value);
            }

            Histos[iSlice][0]= tools.GetHistoBinsNoScale(
                RecoHisto,
                SerialVectorLowBin.at(iSlice),
                SerialVectorHighBin.at(iSlice),
                SerialSliceBinning,
                "Reco"
            );
            Histos[iSlice][1]= tools.GetHistoBinsNoScale(
                RecoTrueHisto,
                SerialVectorLowBin.at(iSlice),
                SerialVectorHighBin.at(iSlice),
                SerialSliceBinning,
                "RecoTrue"
            );
            Histos[iSlice][2]= tools.GetHistoBinsNoScale(
                TrueHisto,
                SerialVectorLowBin.at(iSlice),
                SerialVectorHighBin.at(iSlice),
                SerialSliceBinning,
                "True"
            );
            StartIndex += (SliceNBins + 1);
        }

        // Compute efficiency and purity for forward slice
        double ForwardRecoInt = Histos[1][0]->Integral();
        double ForwardRecoTrueInt = Histos[1][1]->Integral();
        double ForwardTrueInt = Histos[1][2]->Integral();

        std::cout << "Cut position: " << CutPosition << std::endl;
        std::cout << "Reco: " << ForwardRecoInt << std::endl;
        std::cout << "Reco true signal: " << ForwardRecoTrueInt << std::endl;
        std::cout << "True signal: " << ForwardTrueInt << std::endl;

        double Efficiency = ForwardRecoTrueInt / ForwardTrueInt;
        double Purity = ForwardRecoTrueInt / ForwardRecoInt;
        double EffPurity = Efficiency * Purity;

        std::cout << "Efficiency: " << Efficiency << std::endl;
        std::cout << "Purity: " << Purity << std::endl;
        std::cout << "Efficiency * Purity: " << EffPurity << std::endl;
        std::cout << std::endl;

        EfficiencyHisto->Fill(EfficiencyHisto->GetBinCenter(i), Efficiency);
        PurityHisto->Fill(PurityHisto->GetBinCenter(i), Purity);
        EfficiencyPurityHisto->Fill(EfficiencyPurityHisto->GetBinCenter(i), EffPurity);
    }

    EfficiencyHisto->GetXaxis()->SetTitleFont(FontStyle);
    EfficiencyHisto->GetXaxis()->SetLabelFont(FontStyle);
    EfficiencyHisto->GetXaxis()->SetNdivisions(8);
    EfficiencyHisto->GetXaxis()->SetLabelSize(TextSize);
    EfficiencyHisto->GetXaxis()->SetTitleSize(TextSize);
    EfficiencyHisto->GetXaxis()->SetTitleOffset(1.1);
    EfficiencyHisto->GetXaxis()->CenterTitle();
    EfficiencyHisto->GetXaxis()->SetTitle("Cut position");

    EfficiencyHisto->GetYaxis()->SetTitleFont(FontStyle);
    EfficiencyHisto->GetYaxis()->SetLabelFont(FontStyle);
    EfficiencyHisto->GetYaxis()->SetNdivisions(6);
    EfficiencyHisto->GetYaxis()->SetLabelSize(TextSize);
    EfficiencyHisto->GetYaxis()->SetTitleSize(TextSize);
    EfficiencyHisto->GetYaxis()->SetTitleOffset(1.3);
    EfficiencyHisto->GetYaxis()->SetTickSize(0);
    EfficiencyHisto->GetYaxis()->CenterTitle();
    EfficiencyHisto->GetYaxis()->SetTitle("Efficiency");

    EfficiencyHisto->GetYaxis()->SetRangeUser(0.,EfficiencyHisto->GetMaximum()*1.3);
    
    PlotCanvas->cd();
    EfficiencyHisto->SetMarkerColor(kMagenta + 1);
    EfficiencyHisto->SetMarkerStyle(20);
    EfficiencyHisto->SetMarkerSize(1.);
    EfficiencyHisto->Draw("hist p");
    PlotCanvas->SaveAs(dir+"/Figs/CAFAna/CutScan/Efficiency.pdf");

    PurityHisto->GetXaxis()->SetTitleFont(FontStyle);
    PurityHisto->GetXaxis()->SetLabelFont(FontStyle);
    PurityHisto->GetXaxis()->SetNdivisions(8);
    PurityHisto->GetXaxis()->SetLabelSize(TextSize);
    PurityHisto->GetXaxis()->SetTitleSize(TextSize);
    PurityHisto->GetXaxis()->SetTitleOffset(1.1);
    PurityHisto->GetXaxis()->CenterTitle();
    PurityHisto->GetXaxis()->SetTitle("Cut position");

    PurityHisto->GetYaxis()->SetTitleFont(FontStyle);
    PurityHisto->GetYaxis()->SetLabelFont(FontStyle);
    PurityHisto->GetYaxis()->SetNdivisions(6);
    PurityHisto->GetYaxis()->SetLabelSize(TextSize);
    PurityHisto->GetYaxis()->SetTitleSize(TextSize);
    PurityHisto->GetYaxis()->SetTitleOffset(1.3);
    PurityHisto->GetYaxis()->SetTickSize(0);
    PurityHisto->GetYaxis()->CenterTitle();
    PurityHisto->GetYaxis()->SetTitle("Purity");

    PurityHisto->GetYaxis()->SetRangeUser(0.,PurityHisto->GetMaximum()*1.3);
    
    PlotCanvas->cd();
    PurityHisto->SetMarkerColor(kMagenta + 1);
    PurityHisto->SetMarkerStyle(20);
    PurityHisto->SetMarkerSize(1.);
    PurityHisto->Draw("hist p");
    PlotCanvas->SaveAs(dir+"/Figs/CAFAna/CutScan/Purity.pdf");

    EfficiencyPurityHisto->GetXaxis()->SetTitleFont(FontStyle);
    EfficiencyPurityHisto->GetXaxis()->SetLabelFont(FontStyle);
    EfficiencyPurityHisto->GetXaxis()->SetNdivisions(8);
    EfficiencyPurityHisto->GetXaxis()->SetLabelSize(TextSize);
    EfficiencyPurityHisto->GetXaxis()->SetTitleSize(TextSize);
    EfficiencyPurityHisto->GetXaxis()->SetTitleOffset(1.1);
    EfficiencyPurityHisto->GetXaxis()->CenterTitle();
    EfficiencyPurityHisto->GetXaxis()->SetTitle("Cut position");

    EfficiencyPurityHisto->GetYaxis()->SetTitleFont(FontStyle);
    EfficiencyPurityHisto->GetYaxis()->SetLabelFont(FontStyle);
    EfficiencyPurityHisto->GetYaxis()->SetNdivisions(6);
    EfficiencyPurityHisto->GetYaxis()->SetLabelSize(TextSize);
    EfficiencyPurityHisto->GetYaxis()->SetTitleSize(TextSize);
    EfficiencyPurityHisto->GetYaxis()->SetTitleOffset(1.3);
    EfficiencyPurityHisto->GetYaxis()->SetTickSize(0);
    EfficiencyPurityHisto->GetYaxis()->CenterTitle();
    EfficiencyPurityHisto->GetYaxis()->SetTitle("Efficiency times purity");

    EfficiencyPurityHisto->GetYaxis()->SetRangeUser(0.,EfficiencyPurityHisto->GetMaximum()*1.3);
    
    PlotCanvas->cd();
    EfficiencyPurityHisto->SetMarkerColor(kMagenta + 1);
    EfficiencyPurityHisto->SetMarkerStyle(20);
    EfficiencyPurityHisto->SetMarkerSize(1.);
    EfficiencyPurityHisto->Draw("hist p");;
    PlotCanvas->SaveAs(dir+"/Figs/CAFAna/CutScan/EfficiencyPurity.pdf");
}
