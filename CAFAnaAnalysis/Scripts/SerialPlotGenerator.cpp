// ROOT includes.
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>

// std includes.
#include <vector>

// Utils includes.
#include "../../Utils/Constants.h"
#include "../../Utils/Tools.cxx"

using namespace std;
using namespace Constants;

void SerialPlotGenerator() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;	

    Tools tools;

    // Load root file with histograms
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Selection.root";
    std::unique_ptr<TFile> File(TFile::Open(RootFilePath));

    const int NPlots = PlotNames.size();

    // Samples for each plot (all reco, reco true, bkg)
    std::vector<TString> SampleNames; std::vector<TString> Labels; std::vector<int> Colors;
    SampleNames.push_back("_reco"); Labels.push_back("Reconstructed"); Colors.push_back(kBlue+2);
    SampleNames.push_back("_reco_true"); Labels.push_back("True"); Colors.push_back(kRed+1);
    SampleNames.push_back("_bkg"); Labels.push_back("Background"); Colors.push_back(kOrange+7);

    const int NSamples = SampleNames.size();

    for (int iPlot = 0; iPlot < NPlots; iPlot++) {
        
	if ( !(PlotNames[iPlot].Contains("Serial")) ) { continue; }
	// Load true plots
	
        std::vector<TH1D*> TruePlots; TruePlots.resize(NSamples);
        for (int iSample = 0; iSample < NSamples; iSample++) {
            TString PlotName = PlotNames[iPlot] + SampleNames[iSample];
            TruePlots[iSample] = (TH1D*)(File->Get<TH1D>(PlotName));
        }
        
        // Flatten out double differential plots
        auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator["True"+PlotNames[iPlot]+"Plot"];
        auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
        int StartIndex = 0;

        // Create vector to store deserialize plots
        std::vector<std::vector<TH1D*>> Histos;
        Histos.resize(NSlices);
        for (int iSlice = 0; iSlice < NSlices; iSlice++) {
            Histos[iSlice].resize(NSamples);
        }

        // Loop over slices
        for (int iSlice = 0; iSlice < NSlices; iSlice++) {
            TString SlicePlotName = PlotNames[iPlot] + "_" + TString(std::to_string(iSlice));
            double SliceWidth = SliceDiscriminators[iSlice + 1] - SliceDiscriminators[iSlice]; 

            // Get number of bins
            int SliceNBins = SerialVectorBins.at(iSlice);
            std::vector<double> SerialSliceBinning;

            for (int iBin = 0; iBin < SliceNBins + 1; iBin++) {
                double value = SerialVectorRanges.at(StartIndex + iBin);
                SerialSliceBinning.push_back(value);
            } // End of the number of bins and the bin ranges declaration

            // Declare canvas and legend
            TString CanvasName = "Canvas_" + SlicePlotName;
            TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1124,768);

            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);

            TLegend* leg = new TLegend(0.2,0.7,0.75,0.83);
            leg->SetBorderSize(0);
            leg->SetNColumns(2);
            leg->SetTextSize(TextSize*0.8);
            leg->SetTextFont(FontStyle);

            for (int iSample = 0; iSample < NSamples; iSample++) {
                Histos[iSlice][iSample]= tools.GetHistoBins(
                    TruePlots[iSample],
                    SerialVectorLowBin.at(iSlice),
                    SerialVectorHighBin.at(iSlice),
                    SliceWidth,
                    SerialSliceBinning,
                    Labels[iSample]
                );
                Histos[iSlice][iSample]->SetLineWidth(4);
                Histos[iSlice][iSample]->SetLineColor(Colors.at(iSample));

                Histos[iSlice][iSample]->GetXaxis()->SetTitleFont(FontStyle);
                Histos[iSlice][iSample]->GetXaxis()->SetLabelFont(FontStyle);
                Histos[iSlice][iSample]->GetXaxis()->SetNdivisions(8);
                Histos[iSlice][iSample]->GetXaxis()->SetLabelSize(TextSize);
                std::string VarLabel = (std::string) VarLabels.at(iPlot);
                VarLabel.erase(VarLabel.end() - 7, VarLabel.end()); // get rid of (bin #)
                Histos[iSlice][iSample]->GetXaxis()->SetTitle("Reco " + (TString)VarLabel + SerialNameToUnit[PlotNames[iPlot]]);
                Histos[iSlice][iSample]->GetXaxis()->SetTitleSize(TextSize);
                Histos[iSlice][iSample]->GetXaxis()->SetTitleOffset(1.1);
                Histos[iSlice][iSample]->GetXaxis()->CenterTitle();

                Histos[iSlice][iSample]->GetYaxis()->SetTitleFont(FontStyle);
                Histos[iSlice][iSample]->GetYaxis()->SetLabelFont(FontStyle);
                Histos[iSlice][iSample]->GetYaxis()->SetNdivisions(6);
                Histos[iSlice][iSample]->GetYaxis()->SetLabelSize(TextSize);
                Histos[iSlice][iSample]->GetYaxis()->SetTitle("Events");
                Histos[iSlice][iSample]->GetYaxis()->SetTitleSize(TextSize);
                Histos[iSlice][iSample]->GetYaxis()->SetTitleOffset(1.3);
                Histos[iSlice][iSample]->GetYaxis()->SetTickSize(0);
                Histos[iSlice][iSample]->GetYaxis()->CenterTitle();

                double imax = TMath::Max(Histos[iSlice][iSample]->GetMaximum(),Histos[iSlice][0]->GetMaximum());

                double YAxisRange = 1.15*imax;
                Histos[iSlice][iSample]->GetYaxis()->SetRangeUser(0.,YAxisRange);
                Histos[iSlice][0]->GetYaxis()->SetRangeUser(0.,YAxisRange);			

                PlotCanvas->cd();
                Histos[iSlice][iSample]->Draw("hist same");
                Histos[iSlice][0]->Draw("hist same");	

                leg->AddEntry(Histos[iSlice][iSample],Labels[iSample],"l");
            }
            PlotCanvas->cd();
            leg->Draw();

            TLatex *textSlice = new TLatex();
            textSlice->SetTextFont(FontStyle);
            textSlice->SetTextSize(TextSize);
            TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel["True"+PlotNames[iPlot]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
            textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

            TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";
            PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Serial/"+SlicePlotName+".png");
            delete PlotCanvas;

            // Update starting index to move to next slice
            StartIndex += (SliceNBins + 1);
        } // End loop over slices
    } // End loop over plots
}
