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
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"

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

void SelectionEfficiency() {
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	
    
    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // We will create efficiency plots using true variables and definining signal efficiency
    // as the number of reconstructed the events that pass our signal definition and are true
    // signal events over the total true signal events; these two histograms are plotted

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionEfficiency.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    // Construct all spectra
    std::vector<std::tuple<std::unique_ptr<Spectrum>, std::unique_ptr<Spectrum>>> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto TrueSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<2>(Vars.at(i)), kTruthIsSignal, kNoSpillCut);
        auto RecoTrueSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<2>(Vars.at(i)), kTruthIsSignal, kNoSpillCut, kRecoIsSignal);
        Spectra.push_back({std::move(TrueSignals), std::move(RecoTrueSignals)});
    }

    NuLoader.Go();

    for (std::size_t i = 0; i < Vars.size(); i++) {
        auto& [TrueSignals, RecoTrueSignals] = Spectra.at(i);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1D* TrueHisto = TrueSignals->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSignals->ToTH1(TargetPOT);

        // Manage under/overflow bins
        TrueHisto->SetBinContent(TrueHisto->GetNbinsX(), TrueHisto->GetBinContent(TrueHisto->GetNbinsX()) + TrueHisto->GetBinContent(TrueHisto->GetNbinsX() + 1));
        RecoTrueHisto->SetBinContent(RecoTrueHisto->GetNbinsX(), RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX()) + RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX() + 1));

        TrueHisto->SetBinContent(1, TrueHisto->GetBinContent(0) + TrueHisto->GetBinContent(1));
        RecoTrueHisto->SetBinContent(1, RecoTrueHisto->GetBinContent(0) + RecoTrueHisto->GetBinContent(1));

        // Change axes title so efficiency plot inherits them
        TrueHisto->GetYaxis()->SetTitle("Signal efficiency");
        RecoTrueHisto->GetYaxis()->SetTitle("Signal efficiency");

        // Canvas margins
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        if (PlotNames[i].Contains("Serial")) {
            // Flatten out double differential plots
            auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator["True"+PlotNames[i]+"Plot"];
            auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
            int StartIndex = 0;

            // Loop over slices
            for (int iSlice = 0; iSlice < NSlices; iSlice++) {
                // Slice name
                TString SlicePlotName = PlotNames[i] + "_" + TString(std::to_string(iSlice));

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
                TH1D* SlicedTrueHisto = tools.GetHistoBins(
                    TrueHisto,
                    SerialVectorLowBin.at(iSlice),
                    SerialVectorHighBin.at(iSlice),
                    SliceWidth,
                    SerialSliceBinning,
                    "True"
                );
                TH1D* SlicedRecoTrueHisto = tools.GetHistoBins(
                    RecoTrueHisto,
                    SerialVectorLowBin.at(iSlice),
                    SerialVectorHighBin.at(iSlice),
                    SliceWidth,
                    SerialSliceBinning,
                    "RecoTrue"
                );

                // Create efficiency plot
                TEfficiency* Eff = new TEfficiency(*SlicedRecoTrueHisto, *SlicedTrueHisto);
                std::string VarLabel = (std::string) VarLabels.at(i);
                VarLabel.erase(VarLabel.end() - 7, VarLabel.end()); // get rid of (bin #)
                Eff->SetTitle(";True " + (TString)VarLabel + SerialNameToUnit[PlotNames[i]] + ";Efficiency");

                PlotCanvas->cd();
                Eff->SetMarkerStyle(21);
                Eff->SetMarkerColor(kBlack);
                Eff->Draw("AP");
                
                gPad->Update();
                Eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(SerialSliceBinning.at(0),SerialSliceBinning.at(SerialSliceBinning.size() - 1));
                double max = Eff->GetPaintedGraph()->GetYaxis()->GetBinUpEdge(Eff->GetPaintedGraph()->GetYaxis()->GetNbins());
                Eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., max*1.2);

                Eff->GetPaintedGraph()->GetXaxis()->CenterTitle();
                Eff->GetPaintedGraph()->GetXaxis()->SetTitleFont(FontStyle);
                Eff->GetPaintedGraph()->GetXaxis()->SetLabelFont(FontStyle);
                Eff->GetPaintedGraph()->GetXaxis()->SetNdivisions(8);
                Eff->GetPaintedGraph()->GetXaxis()->SetLabelSize(TextSize);
                Eff->GetPaintedGraph()->GetXaxis()->SetTitleSize(TextSize);
                Eff->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.1);

                Eff->GetPaintedGraph()->GetYaxis()->CenterTitle();
                Eff->GetPaintedGraph()->GetYaxis()->SetTitleFont(FontStyle);
                Eff->GetPaintedGraph()->GetYaxis()->SetLabelFont(FontStyle);
                Eff->GetPaintedGraph()->GetYaxis()->SetNdivisions(6);
                Eff->GetPaintedGraph()->GetYaxis()->SetLabelSize(TextSize);
                Eff->GetPaintedGraph()->GetYaxis()->SetTitleSize(TextSize);
                Eff->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.3);
                Eff->GetPaintedGraph()->GetYaxis()->SetTickSize(0);

                // Slice label
                TLatex *textSlice = new TLatex();
                TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel["True"+PlotNames[i]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
                textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

                // Save as pdf
                PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Efficiency/"+SlicePlotName+".pdf");

                // Save to root file
                SaveFile->WriteObject(Eff, SlicePlotName+"_eff");

                StartIndex += (SliceNBins + 1);
            }
        } else {
            // Create efficiency plot
            TEfficiency* Eff = new TEfficiency(*RecoTrueHisto, *TrueHisto);
            Eff->SetTitle((";True " + VarLabels.at(i) + ";Efficiency").c_str());

            PlotCanvas->cd();
            Eff->SetMarkerStyle(21);
            Eff->SetMarkerColor(kBlack);
            Eff->Draw("AP");
            
            gPad->Update();
            Eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(VarBins.at(i).Min(),VarBins.at(i).Max());
            double max = Eff->GetPaintedGraph()->GetYaxis()->GetBinUpEdge(Eff->GetPaintedGraph()->GetYaxis()->GetNbins());
            Eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0., max*1.2);

            Eff->GetPaintedGraph()->GetXaxis()->CenterTitle();
            Eff->GetPaintedGraph()->GetXaxis()->SetTitleFont(FontStyle);
            Eff->GetPaintedGraph()->GetXaxis()->SetLabelFont(FontStyle);
            Eff->GetPaintedGraph()->GetXaxis()->SetNdivisions(8);
            Eff->GetPaintedGraph()->GetXaxis()->SetLabelSize(TextSize);
            Eff->GetPaintedGraph()->GetXaxis()->SetTitleSize(TextSize);
            Eff->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.1);

            Eff->GetPaintedGraph()->GetYaxis()->CenterTitle();
            Eff->GetPaintedGraph()->GetYaxis()->SetTitleFont(FontStyle);
            Eff->GetPaintedGraph()->GetYaxis()->SetLabelFont(FontStyle);
            Eff->GetPaintedGraph()->GetYaxis()->SetNdivisions(6);
            Eff->GetPaintedGraph()->GetYaxis()->SetLabelSize(TextSize);
            Eff->GetPaintedGraph()->GetYaxis()->SetTitleSize(TextSize);
            Eff->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.3);
            Eff->GetPaintedGraph()->GetYaxis()->SetTickSize(0);

            // Save as pdf
            PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Efficiency/"+PlotNames[i]+".pdf");

            // Save to root file
            SaveFile->WriteObject(Eff, PlotNames[i]+"_eff");
        }
        delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();
}
