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
#include "TLatex.h"
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

void SelectionTopologyBreakdown() {
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
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionTopologyBreakdown.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    //////////////
    // Topologies
    //////////////

    std::vector<std::tuple<std::string, TruthCut>> Topologies = {
        {"CC2p0#pi", kTruthIsSignal},
        {"CC1p0#pi", kCC1p0pi},
        {"CC(N > 2)p0#pi", kCCNg2p0pi},
        {"CC(N #geq 0)p1#pi", kCCNgt0p1pi},
        {"CC0p0#pi", kCC0p0pi},
        {"Other", kOtherTopology}
        // {"CC(N #geq 0)p(N > 1)#pi", kCCNgt0pNg1pi},
    };
    std::vector<int> Colors{kBlack, kBlue, kRed+1, kOrange+7, kGreen+1, kMagenta+1, kCyan+2};

    std::vector<std::vector<std::unique_ptr<Spectrum>>> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        std::vector<std::unique_ptr<Spectrum>> InnerSpectra;

        // Without any topology discrimination
        auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsSignal); 
        InnerSpectra.push_back(std::move(RecoSignals));

        for (std::size_t j = 0; j < Topologies.size(); j++) {
            TruthCut TopologyCut = std::get<1>(Topologies[j]);

            const Cut kRecoSignalsCut([=](const caf::SRSliceProxy* slc) {
                return kRecoIsSignal(slc) && TopologyCut(&slc->truth);
            });
            auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoSignalsCut); 
            InnerSpectra.push_back(std::move(RecoSignals));
        }
        Spectra.push_back(std::move(InnerSpectra));
    }

    NuLoader.Go();

    for (std::size_t iVar = 0; iVar < Vars.size(); iVar++) {
        if (PlotNames[iVar].Contains("Serial")) {
            // Flatten out double differential plots
            auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator["True"+PlotNames[iVar]+"Plot"];
            auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
            int StartIndex = 0;

            // Create vector to store deserialize plots
            std::vector<std::vector<TH1D*>> Histos;
            Histos.resize(NSlices);
            for (int iSlice = 0; iSlice < NSlices; iSlice++) {
                Histos[iSlice].resize(Topologies.size());
            }

            for (int iSlice = 0; iSlice < NSlices; iSlice++) {
                TString SlicePlotName = PlotNames[iVar] + "_" + TString(std::to_string(iSlice));
                double SliceWidth = SliceDiscriminators[iSlice + 1] - SliceDiscriminators[iSlice]; 

                // Get number of bins
                int SliceNBins = SerialVectorBins.at(iSlice);
                std::vector<double> SerialSliceBinning;

                for (int iBin = 0; iBin < SliceNBins + 1; iBin++) {
                    double value = SerialVectorRanges.at(StartIndex + iBin);
                    SerialSliceBinning.push_back(value);
                } 

                // Declare canvas and legend
                TString CanvasName = "Canvas_" + SlicePlotName;
                TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1124,768);

                PlotCanvas->SetTopMargin(0.13);
                PlotCanvas->SetLeftMargin(0.17);
                PlotCanvas->SetRightMargin(0.05);
                PlotCanvas->SetBottomMargin(0.16);

                TLegend* leg = new TLegend(0.2,0.68,0.87,0.83);
                leg->SetBorderSize(0);
                leg->SetNColumns(2);
                leg->SetTextSize(TextSize*0.7);
                leg->SetTextFont(FontStyle);

                for (std::size_t iTop = 0; iTop < Topologies.size() + 1; iTop++) {
                    auto& TopRecoSignals = Spectra[iVar][iTop];
                    Histos[iSlice][iTop]= tools.GetHistoBins(
                        TopRecoSignals->ToTH1(TargetPOT),
                        SerialVectorLowBin.at(iSlice),
                        SerialVectorHighBin.at(iSlice),
                        SliceWidth,
                        SerialSliceBinning,
                        VarLabels[iTop]
                    );
                    Histos[iSlice][iTop]->SetLineWidth(4);
                    Histos[iSlice][iTop]->SetLineColor(Colors.at(iTop));

                    Histos[iSlice][iTop]->GetXaxis()->SetTitleFont(FontStyle);
                    Histos[iSlice][iTop]->GetXaxis()->SetLabelFont(FontStyle);
                    Histos[iSlice][iTop]->GetXaxis()->SetNdivisions(8);
                    Histos[iSlice][iTop]->GetXaxis()->SetLabelSize(TextSize);
                    std::string VarLabel = (std::string) VarLabels.at(iVar);
                    VarLabel.erase(VarLabel.end() - 7, VarLabel.end()); // get rid of (bin #)
                    Histos[iSlice][iTop]->GetXaxis()->SetTitle("Reco " + (TString)VarLabel + SerialNameToUnit[PlotNames[iVar]]);
                    Histos[iSlice][iTop]->GetXaxis()->SetTitleSize(TextSize);
                    Histos[iSlice][iTop]->GetXaxis()->SetTitleOffset(1.1);
                    Histos[iSlice][iTop]->GetXaxis()->CenterTitle();

                    Histos[iSlice][iTop]->GetYaxis()->SetTitleFont(FontStyle);
                    Histos[iSlice][iTop]->GetYaxis()->SetLabelFont(FontStyle);
                    Histos[iSlice][iTop]->GetYaxis()->SetNdivisions(6);
                    Histos[iSlice][iTop]->GetYaxis()->SetLabelSize(TextSize);
                    Histos[iSlice][iTop]->GetYaxis()->SetTitleSize(TextSize);
                    Histos[iSlice][iTop]->GetYaxis()->SetTitleOffset(1.3);
                    Histos[iSlice][iTop]->GetYaxis()->SetTickSize(0);
                    Histos[iSlice][iTop]->GetYaxis()->CenterTitle();

                    double imax = Histos[iSlice][0]->GetMaximum();
                    double YAxisRange = 1.4*imax;
                    Histos[iSlice][iTop]->GetYaxis()->SetRangeUser(0.,YAxisRange);

                    double frac = Histos[iSlice][iTop]->Integral() / Histos[iSlice][0]->Integral() * 100.;
                    std::string TopLabel = (iTop == 0) ? "All" : std::get<0>(Topologies[iTop - 1]);
                    TString LegLabel = (TString)TopLabel + " (" + tools.to_string_with_precision(frac,1) + "%)";
                    TLegendEntry* legReco = leg->AddEntry(Histos[iSlice][iTop],LegLabel,"l");

                    PlotCanvas->cd();
                    Histos[iSlice][iTop]->Draw("hist same");
                    Histos[iSlice][0]->Draw("hist same");

                    // Save to root file
                    SaveFile->WriteObject(Histos[iSlice][iTop], SlicePlotName+(TString)LegLabel+"_reco");
                }
                leg->Draw();

                TLatex *textSlice = new TLatex();
                textSlice->SetTextFont(FontStyle);
                textSlice->SetTextSize(TextSize);
                TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel["True"+PlotNames[iVar]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
                textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

                // Save as png
                PlotCanvas->SaveAs(dir+"/Figs/CAFAna/TopologyBreakdown/"+SlicePlotName+".png");

                delete PlotCanvas;

                StartIndex += (SliceNBins + 1);
            }
        } else {
            std::vector<TH1D*> Histos; Histos.resize(Topologies.size() + 1);
            TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);

            TLegend* leg = new TLegend(0.2,0.68,0.87,0.83);
            leg->SetBorderSize(0);
            leg->SetNColumns(2);
            leg->SetTextSize(TextSize*0.7);
            leg->SetTextFont(FontStyle);

            for (std::size_t iTop = 0; iTop < Topologies.size() + 1; iTop++) {
                auto& TopRecoSignals = Spectra[iVar][iTop];
                Histos[iTop] = TopRecoSignals->ToTH1(TargetPOT);

                // Manage under/overflow bins
                Histos[iTop]->SetBinContent(Histos[iTop]->GetNbinsX(), Histos[iTop]->GetBinContent(Histos[iTop]->GetNbinsX()) + Histos[iTop]->GetBinContent(Histos[iTop]->GetNbinsX() + 1));
                Histos[iTop]->SetBinContent(1, Histos[iTop]->GetBinContent(0) + Histos[iTop]->GetBinContent(1));

                double frac = Histos[iTop]->Integral() / Histos[0]->Integral() * 100.;
                std::string TopLabel = (iTop == 0) ? "All" : std::get<0>(Topologies[iTop - 1]);
                TString LegLabel = (TString)TopLabel + " (" + tools.to_string_with_precision(frac,1) + "%)";
                TLegendEntry* legReco = leg->AddEntry(Histos[iTop],LegLabel,"l");
                Histos[iTop]->SetLineColor(Colors.at(iTop));
                Histos[iTop]->SetLineWidth(4);

                // Style histograms
                if (iTop == 0) {
                    Histos[iTop]->GetXaxis()->SetTitleFont(FontStyle);
                    Histos[iTop]->GetXaxis()->SetLabelFont(FontStyle);
                    Histos[iTop]->GetXaxis()->SetNdivisions(8);
                    Histos[iTop]->GetXaxis()->SetLabelSize(TextSize);
                    Histos[iTop]->GetXaxis()->SetTitleSize(TextSize);
                    Histos[iTop]->GetXaxis()->SetTitleOffset(1.1);
                    Histos[iTop]->GetXaxis()->CenterTitle();
                    Histos[iTop]->GetXaxis()->SetTitle(("Reco " + VarLabels.at(iVar)).c_str());

                    Histos[iTop]->GetYaxis()->SetTitleFont(FontStyle);
                    Histos[iTop]->GetYaxis()->SetLabelFont(FontStyle);
                    Histos[iTop]->GetYaxis()->SetNdivisions(6);
                    Histos[iTop]->GetYaxis()->SetLabelSize(TextSize);
                    Histos[iTop]->GetYaxis()->SetTitleSize(TextSize);
                    Histos[iTop]->GetYaxis()->SetTitleOffset(1.3);
                    Histos[iTop]->GetYaxis()->SetTickSize(0);
                    Histos[iTop]->GetYaxis()->CenterTitle();
                }
                double imax = Histos[0]->GetMaximum();
                double YAxisRange = 1.4*imax;
                Histos[iTop]->GetYaxis()->SetRangeUser(0.,YAxisRange);

                PlotCanvas->cd();
                Histos[iTop]->Draw("hist same");

                // Save to root file
                SaveFile->WriteObject(Histos[iTop], PlotNames[iVar]+(TString)LegLabel+"_reco");
            }
            leg->Draw();

            // Save as png
            PlotCanvas->SaveAs(dir+"/Figs/CAFAna/TopologyBreakdown/"+PlotNames[iVar]+".png");

            delete PlotCanvas;
        }
    }
}
