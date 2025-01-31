// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"

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
#include "Utils.h"
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionMigrationMatrix() {
    int FontStyle = 132;
    double TextSize = 0.06;	

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // We now create migration matrices for several variables, where these represent the ratio
    // between the events with a reconstructed value in bin i over the total true signal events
    // with true value in bin i

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Matrix.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    //Define set of analysis specific variables
    static const std::vector<std::tuple<Var, Var, TruthVar>> Vars=GetSISVars();
    static const std::vector<Binning> VarBins = GetSISBins();
    static const std::vector<std::string> VarLabels= GetSISVarLabels();
    static const std::vector<std::string> varnamesshort=GetSISVarNames();
    
    //Part thats really for the selection efficency 
    std::cout<<Vars.size()<<", "<<VarBins.size()<< ", "<<VarLabels.size()<<endl;

    // Construct spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>, 
        std::unique_ptr<Spectrum>,
        std::vector<std::unique_ptr<Spectrum>>
    >> Spectra;
    for (std::size_t i = 0; i < Vars.size(); i++) {
        const std::vector<double>& BinEdges = VarBins.at(i).Edges();
        std::vector<std::unique_ptr<Spectrum>> InnerSpectra;
        Var kCurrentVar = std::get<1>(Vars.at(i));

        for (int j = 0; j < VarBins.at(i).NBins(); j++) {
            double BinMin = BinEdges.at(j);
            double BinMax = (j == VarBins.at(i).NBins() - 1) ?  VarBins.at(i).Max() : BinEdges.at(j + 1);

            const Cut TempCut([=](const caf::SRSliceProxy* slc) {
                return (
                    kRecoIsTrueReco(slc) && 
                    kCurrentVar(slc) >= BinMin &&
                    kCurrentVar(slc) < BinMax
                );
            });
            auto RecoBinValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, TempCut);
            InnerSpectra.push_back(std::move(RecoBinValues));
        }
        auto RecoTruthValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<1>(Vars.at(i)), kNoSpillCut, kRecoIsTrueReco);
        auto TruthValues = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<2>(Vars.at(i)), kTruthIsSignal, kNoSpillCut);
        Spectra.push_back({std::move(TruthValues), std::move(RecoTruthValues), std::move(InnerSpectra)});
    }

    NuLoader.Go();

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {
        double RecoTotalEvents = 0.;
        auto& [TruthValues, RecoTruthValues, InnerSpectra] = Spectra.at(i);

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        TH1* TruthValuesHist = TruthValues->ToTH1(TargetPOT);
        TH1* RecoTruthValuesHist = RecoTruthValues->ToTH1(TargetPOT);

        // Manage under/overflow bins
        TruthValuesHist->SetBinContent(TruthValuesHist->GetNbinsX(), TruthValuesHist->GetBinContent(TruthValuesHist->GetNbinsX()) + TruthValuesHist->GetBinContent(TruthValuesHist->GetNbinsX() + 1));
        RecoTruthValuesHist->SetBinContent(RecoTruthValuesHist->GetNbinsX(), RecoTruthValuesHist->GetBinContent(RecoTruthValuesHist->GetNbinsX()) + RecoTruthValuesHist->GetBinContent(RecoTruthValuesHist->GetNbinsX() + 1));

        TruthValuesHist->SetBinContent(1, TruthValuesHist->GetBinContent(0) + TruthValuesHist->GetBinContent(1));
        RecoTruthValuesHist->SetBinContent(1, RecoTruthValuesHist->GetBinContent(0) + RecoTruthValuesHist->GetBinContent(1));

        // Get bins for matrices
        const int NBins = VarBins.at(i).NBins();
        const std::vector<double>& BinEdges = VarBins.at(i).Edges();

        TH2* MigrationMatrix = new TH2D(
            "Migration",
            "Migration",
            NBins, BinEdges.data(),
            NBins, BinEdges.data()
        );
        TH2* ResponseMatrix = new TH2D(
            "Response",
            "Response",
            NBins, BinEdges.data(),
            NBins, BinEdges.data()
        );

        std::vector<TH1*> RecoValuesHistos;
        for (int j = 0; j < VarBins.at(i).NBins(); j++) {
            RecoValuesHistos.push_back(InnerSpectra.at(j)->ToTH1(TargetPOT));
        }

        // Debugging
        std::cout << PlotNames.at(i) << std::endl;

        for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
            double RecoTruthCounts = RecoTruthValuesHist->GetBinContent(x);
            double TruthCounts = TruthValuesHist->GetBinContent(x);
            TH1* RecoValuesHist = RecoValuesHistos.at(x - 1); // -1 because ROOT labels bins starting from 1

            // Manage under/overflow bins
            RecoValuesHist->SetBinContent(RecoValuesHist->GetNbinsX(), RecoValuesHist->GetBinContent(RecoValuesHist->GetNbinsX()) + RecoValuesHist->GetBinContent(RecoValuesHist->GetNbinsX() + 1));
            RecoValuesHist->SetBinContent(1, RecoValuesHist->GetBinContent(0) + RecoValuesHist->GetBinContent(1));

            for (int y = 1; y < VarBins.at(i).NBins() + 1; y++) {
                double RecoCounts = RecoValuesHist->GetBinContent(y);

                // Migration matrix
                double MigrationRatio = RecoCounts / RecoTruthCounts;
                if (RecoTruthCounts == 0.) MigrationRatio = 0.0;

                // Response matrix
                double ResponseRatio = RecoCounts / TruthCounts;
                if (TruthCounts == 0.) ResponseRatio = 0.0;

                // Debugging
                std::cout << "Reco low bin: " << y << ". Counts: " << RecoCounts << std::endl;
                std::cout << "True low bin: " << x << ". Counts: " << TruthCounts << std::endl;
                std::cout << "Reco truth counts: " << RecoTruthCounts << std::endl;
                std::cout << "Migration ratio: " << MigrationRatio << std::endl;
                std::cout << "Response ratio: " << ResponseRatio << std::endl;
                std::cout << std::endl;

                MigrationMatrix->Fill(
                    TruthValuesHist->GetXaxis()->GetBinCenter(x),
                    RecoValuesHist->GetXaxis()->GetBinCenter(y),
                    MigrationRatio
                );
                ResponseMatrix->Fill(
                    TruthValuesHist->GetXaxis()->GetBinCenter(x),
                    RecoValuesHist->GetXaxis()->GetBinCenter(y),
                    ResponseRatio
                );
            }
            RecoTotalEvents += RecoValuesHist->Integral();
        }
        // Sanity check, these should be the same
        std::cout << "Total true var events: " << RecoTruthValuesHist->Integral() << std::endl;
        std::cout << "Total reco var events: " << RecoTotalEvents << std::endl;

        MigrationMatrix->GetXaxis()->SetTitle(("True " + VarLabels.at(i)).c_str());
        MigrationMatrix->GetXaxis()->CenterTitle();
        MigrationMatrix->GetXaxis()->SetTitleOffset(1.1);
        MigrationMatrix->GetXaxis()->SetTitleFont(FontStyle);
        MigrationMatrix->GetXaxis()->SetTitleSize(TextSize);
        MigrationMatrix->GetXaxis()->SetLabelFont(FontStyle);
        MigrationMatrix->GetXaxis()->SetLabelSize(TextSize);
        MigrationMatrix->GetXaxis()->SetNdivisions(5);

        MigrationMatrix->GetYaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());
        MigrationMatrix->GetYaxis()->CenterTitle();
        MigrationMatrix->GetYaxis()->SetTitleOffset(1.1);
        MigrationMatrix->GetYaxis()->SetTitleFont(FontStyle);
        MigrationMatrix->GetYaxis()->SetTitleSize(TextSize);
        MigrationMatrix->GetYaxis()->SetLabelFont(FontStyle);
        MigrationMatrix->GetYaxis()->SetLabelSize(TextSize);
        MigrationMatrix->GetYaxis()->SetNdivisions(5);

        double MigrationMin = MigrationMatrix->GetMinimum();
        double MigrationMax = MigrationMatrix->GetMaximum();
        MigrationMatrix->GetZaxis()->SetRangeUser(MigrationMin,MigrationMax);
        MigrationMatrix->GetZaxis()->CenterTitle();
        MigrationMatrix->GetZaxis()->SetTitleFont(FontStyle);
        MigrationMatrix->GetZaxis()->SetTitleSize(TextSize);
        MigrationMatrix->GetZaxis()->SetLabelFont(FontStyle);
        MigrationMatrix->GetZaxis()->SetLabelSize(TextSize);
        MigrationMatrix->GetZaxis()->SetNdivisions(5);

        ResponseMatrix->GetXaxis()->SetTitle(("True " + VarLabels.at(i)).c_str());
        ResponseMatrix->GetXaxis()->CenterTitle();
        ResponseMatrix->GetXaxis()->SetTitleOffset(1.1);
        ResponseMatrix->GetXaxis()->SetTitleFont(FontStyle);
        ResponseMatrix->GetXaxis()->SetTitleSize(TextSize);
        ResponseMatrix->GetXaxis()->SetLabelFont(FontStyle);
        ResponseMatrix->GetXaxis()->SetLabelSize(TextSize);
        ResponseMatrix->GetXaxis()->SetNdivisions(5);

        ResponseMatrix->GetYaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());
        ResponseMatrix->GetYaxis()->CenterTitle();
        ResponseMatrix->GetYaxis()->SetTitleOffset(1.1);
        ResponseMatrix->GetYaxis()->SetTitleFont(FontStyle);
        ResponseMatrix->GetYaxis()->SetTitleSize(TextSize);
        ResponseMatrix->GetYaxis()->SetLabelFont(FontStyle);
        ResponseMatrix->GetYaxis()->SetLabelSize(TextSize);
        ResponseMatrix->GetYaxis()->SetNdivisions(5);

        double ResponseMin = ResponseMatrix->GetMinimum();
        double ResponseMax = ResponseMatrix->GetMaximum();
        ResponseMatrix->GetZaxis()->SetRangeUser(ResponseMin,ResponseMax);
        ResponseMatrix->GetZaxis()->CenterTitle();
        ResponseMatrix->GetZaxis()->SetTitleFont(FontStyle);
        ResponseMatrix->GetZaxis()->SetTitleSize(TextSize);
        ResponseMatrix->GetZaxis()->SetLabelFont(FontStyle);
        ResponseMatrix->GetZaxis()->SetLabelSize(TextSize);
        ResponseMatrix->GetZaxis()->SetNdivisions(5);

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        PlotCanvas->cd();
        MigrationMatrix->Draw("colz");

        // Save migration matrix as pdf
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Matrices/Migration"+PlotNames[i]+".pdf");

        PlotCanvas->cd();
        ResponseMatrix->Draw("colz");

        // Save response matrix as pdf
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Matrices/Response"+PlotNames[i]+".pdf");
        
        // Save both to root file
        SaveFile->WriteObject(MigrationMatrix, PlotNames[i]+"_migration");
        SaveFile->WriteObject(ResponseMatrix, PlotNames[i]+"_response");
        SaveFile->WriteObject(TruthValuesHist, PlotNames[i]+"_true");

        delete PlotCanvas;
        delete MigrationMatrix;
        delete ResponseMatrix;
    }
    // Close file
    SaveFile->Close();
}
