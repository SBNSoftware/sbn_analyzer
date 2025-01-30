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

void SelectionInteBreakdown() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;	

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionInteBreakdown.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    ////////////////////
    // Interaction modes
    ////////////////////

    std::vector<std::tuple<std::string, int>> IntModes = {
        {"QE", 0},
        {"MEC", 10},
        {"RES", 1},
        {"DIS", 2}
    };
    std::vector<int> Colors{kBlack, kBlue, kRed+1, kOrange+7, kGreen+1};

    std::vector<std::vector<std::unique_ptr<Spectrum>>> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        std::vector<std::unique_ptr<Spectrum>> InnerSpectra;

        // Without any interaction discrimination
        auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsSignal); 
        InnerSpectra.push_back(std::move(RecoSignals));

        for (std::size_t j = 0; j < IntModes.size(); j++) {
            int IntCode = std::get<1>(IntModes[j]);

            const Cut kRecoSignalsCut([=](const caf::SRSliceProxy* slc) {
                return (
                    kRecoIsSignal(slc) &&
                    slc->truth.genie_mode == IntCode
                );
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
                Histos[iSlice].resize(IntModes.size());
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

                TLegend* leg = new TLegend(0.2,0.73,0.85,0.83);
                leg->SetBorderSize(0);
                leg->SetNColumns(3);
                leg->SetTextSize(TextSize*0.8);
                leg->SetTextFont(FontStyle);

                for (std::size_t iInt = 0; iInt < IntModes.size() + 1; iInt++) {
                    auto& IntRecoSignals = Spectra[iVar][iInt];
                    Histos[iSlice][iInt]= tools.GetHistoBins(
                        IntRecoSignals->ToTH1(TargetPOT),
                        SerialVectorLowBin.at(iSlice),
                        SerialVectorHighBin.at(iSlice),
                        SliceWidth,
                        SerialSliceBinning,
                        VarLabels[iInt]
                    );
                    Histos[iSlice][iInt]->SetLineWidth(4);
                    Histos[iSlice][iInt]->SetLineColor(Colors.at(iInt));

                    Histos[iSlice][iInt]->GetXaxis()->SetTitleFont(FontStyle);
                    Histos[iSlice][iInt]->GetXaxis()->SetLabelFont(FontStyle);
                    Histos[iSlice][iInt]->GetXaxis()->SetNdivisions(8);
                    Histos[iSlice][iInt]->GetXaxis()->SetLabelSize(TextSize);
                    std::string VarLabel = (std::string) VarLabels.at(iVar);
                    VarLabel.erase(VarLabel.end() - 7, VarLabel.end()); // get rid of (bin #)
                    Histos[iSlice][iInt]->GetXaxis()->SetTitle("Reco " + (TString)VarLabel + SerialNameToUnit[PlotNames[iVar]]);
                    Histos[iSlice][iInt]->GetXaxis()->SetTitleSize(TextSize);
                    Histos[iSlice][iInt]->GetXaxis()->SetTitleOffset(1.1);
                    Histos[iSlice][iInt]->GetXaxis()->CenterTitle();

                    Histos[iSlice][iInt]->GetYaxis()->SetTitleFont(FontStyle);
                    Histos[iSlice][iInt]->GetYaxis()->SetLabelFont(FontStyle);
                    Histos[iSlice][iInt]->GetYaxis()->SetNdivisions(6);
                    Histos[iSlice][iInt]->GetYaxis()->SetLabelSize(TextSize);
                    Histos[iSlice][iInt]->GetYaxis()->SetTitleSize(TextSize);
                    Histos[iSlice][iInt]->GetYaxis()->SetTitleOffset(1.3);
                    Histos[iSlice][iInt]->GetYaxis()->SetTickSize(0);
                    Histos[iSlice][iInt]->GetYaxis()->CenterTitle();

                    double imax = TMath::Max(Histos[iSlice][iInt]->GetMaximum(),Histos[iSlice][0]->GetMaximum());
                    double YAxisRange = 1.15*imax;
                    Histos[iSlice][iInt]->GetYaxis()->SetRangeUser(0.,YAxisRange);
                    Histos[iSlice][0]->GetYaxis()->SetRangeUser(0.,YAxisRange);			

                    double frac = Histos[iSlice][iInt]->Integral() / Histos[iSlice][0]->Integral() * 100.;
                    std::string IntLabel = (iInt == 0) ? "All" : std::get<0>(IntModes[iInt - 1]);
                    TString LegLabel = (TString)IntLabel + " (" + tools.to_string_with_precision(frac,1) + "%)";
                    TLegendEntry* legReco = leg->AddEntry(Histos[iSlice][iInt],LegLabel,"l");

                    PlotCanvas->cd();
                    Histos[iSlice][iInt]->Draw("hist same");
                    Histos[iSlice][0]->Draw("hist same");

                    // Save to root file
                    SaveFile->WriteObject(Histos[iSlice][iInt], SlicePlotName+(TString)LegLabel+"_reco");
                }
                leg->Draw();

                TLatex *textSlice = new TLatex();
                textSlice->SetTextFont(FontStyle);
                textSlice->SetTextSize(TextSize);
                TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel["True"+PlotNames[iVar]+"Plot"] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
                textSlice->DrawLatexNDC(0.4,0.92,SliceLabel);

                // Save as pdf
                PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/InteBreakdown/"+SlicePlotName+".pdf");

                delete PlotCanvas;

                StartIndex += (SliceNBins + 1);
            }
        } else {
            std::vector<TH1D*> Histos; Histos.resize(IntModes.size() + 1);
            TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);

            TLegend* leg = new TLegend(0.2,0.73,0.85,0.83);
            leg->SetBorderSize(0);
            leg->SetNColumns(3);
            leg->SetTextSize(TextSize*0.8);
            leg->SetTextFont(FontStyle);

            for (std::size_t iInt = 0; iInt < IntModes.size() + 1; iInt++) {
                auto& IntRecoSignals = Spectra[iVar][iInt];
                Histos[iInt] = IntRecoSignals->ToTH1(TargetPOT);

                // Manage under/overflow bins
                Histos[iInt]->SetBinContent(Histos[iInt]->GetNbinsX(), Histos[iInt]->GetBinContent(Histos[iInt]->GetNbinsX()) + Histos[iInt]->GetBinContent(Histos[iInt]->GetNbinsX() + 1));
                Histos[iInt]->SetBinContent(1, Histos[iInt]->GetBinContent(0) + Histos[iInt]->GetBinContent(1));

                double frac = Histos[iInt]->Integral() / Histos[0]->Integral() * 100.;
                std::string IntLabel = (iInt == 0) ? "All" : std::get<0>(IntModes[iInt - 1]);
                TString LegLabel = (TString)IntLabel + " (" + tools.to_string_with_precision(frac,1) + "%)";
                TLegendEntry* legReco = leg->AddEntry(Histos[iInt],LegLabel,"l");
                Histos[iInt]->SetLineColor(Colors.at(iInt));
                Histos[iInt]->SetLineWidth(4);

                // Style histograms
                if (iInt == 0) {
                    Histos[iInt]->GetXaxis()->SetTitleFont(FontStyle);
                    Histos[iInt]->GetXaxis()->SetLabelFont(FontStyle);
                    Histos[iInt]->GetXaxis()->SetNdivisions(8);
                    Histos[iInt]->GetXaxis()->SetLabelSize(TextSize);
                    Histos[iInt]->GetXaxis()->SetTitleSize(TextSize);
                    Histos[iInt]->GetXaxis()->SetTitleOffset(1.1);
                    Histos[iInt]->GetXaxis()->CenterTitle();
                    Histos[iInt]->GetXaxis()->SetTitle(("Reco " + VarLabels.at(iVar)).c_str());

                    Histos[iInt]->GetYaxis()->SetTitleFont(FontStyle);
                    Histos[iInt]->GetYaxis()->SetLabelFont(FontStyle);
                    Histos[iInt]->GetYaxis()->SetNdivisions(6);
                    Histos[iInt]->GetYaxis()->SetLabelSize(TextSize);
                    Histos[iInt]->GetYaxis()->SetTitleSize(TextSize);
                    Histos[iInt]->GetYaxis()->SetTitleOffset(1.3);
                    Histos[iInt]->GetYaxis()->SetTickSize(0);
                    Histos[iInt]->GetYaxis()->CenterTitle();
                }
                double imax = Histos[0]->GetMaximum();
                double YAxisRange = 1.3*imax;
                Histos[iInt]->GetYaxis()->SetRangeUser(0.,YAxisRange);

                PlotCanvas->cd();
                Histos[iInt]->Draw("hist same");

                // Save to root file
                SaveFile->WriteObject(Histos[iInt], PlotNames[iVar]+(TString)LegLabel+"_reco");
            }
            leg->Draw();

            // Save as pdf
            PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/InteBreakdown/"+PlotNames[iVar]+".pdf");

            delete PlotCanvas;
        }

        
    }
}
