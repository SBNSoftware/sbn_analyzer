#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

#include "../../Utils/Constants.h"
#include "../../Utils/Tools.cxx"

using namespace std;
using namespace Constants;

void GeneratorInteBreakDown() {

    //------------------------------//

    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;			

    TString OutFilePath = "/pnfs/sbnd/persistent/users/" + (TString)UserName + "/HighSamples/FlatTree/";

    Tools tools;

    //------------------------------//

    // Event generators

    std::vector<TString> Names; std::vector<TString> Labels; 
    std::vector<TString> Process{"","QE","MEC","RES","DIS"};
    std::vector<int> Colors{kBlack,kBlue,kRed+1,kOrange+7,kGreen+1,kMagenta+1};

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_AR23.root"); 
    Labels.push_back("GENIE_AR23");

    // Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_AR23_Emp2015.root"); 
    // Labels.push_back("GENIE_AR23_Emp2015");
    // Colors.push_back(kBlue+16);

    // Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_AR23_Nieves2016.root"); 
    // Labels.push_back("GENIE_AR23_Nieves2016");
    // Colors.push_back(kBlue+24);

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_G18.root"); 
    Labels.push_back("GENIE_G18");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NuWro.root"); 
    Labels.push_back("NuWro");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT.root"); 
    Labels.push_back("NEUT");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU.root"); 
    Labels.push_back("GiBUU");

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_NoFSI.root"); 
    Labels.push_back("GiBUU_NoFSI");

    const int NSamples = Names.size();
    const int NColors = Colors.size();
    const int NProcesses = Process.size();

    // Sanity check
    if (NColors < NProcesses) { cout << "Give me some more colors for all the processes!" << endl; return; }

    std::vector<TFile*> Files; Files.resize(NSamples);

    //------------------------------//

    // Plots to overlay

    std::vector<TString> PlotNames;
    std::vector<TString> YAxisLabel;
    std::vector<TString> XAxisLabel;

    //------------------------------//

    // Pre FSI
    PlotNames.push_back("TrueNoFSILeadingProtonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L}}");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIRecoilProtonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{R}}");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSILeadingProtonMomentumPlot");
    XAxisLabel.push_back("|#vec{p}_{L}| [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{L}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIRecoilProtonMomentumPlot");
    XAxisLabel.push_back("|#vec{p}_{R}| [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{R}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIMuonMomentumPlot");
    XAxisLabel.push_back("|#vec{p}_{#mu}| [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{#mu}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSICosOpeningAngleProtonsPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSICosOpeningAngleMuonTotalProtonPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSITransverseMomentumPlot");
    XAxisLabel.push_back("#delta P_{T} [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIDeltaAlphaTPlot");
    XAxisLabel.push_back("#delta #alpha_{T} [deg]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta #alpha_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIInvariantMassPlot");
    XAxisLabel.push_back("W [GeV]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta W} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSICosOpeningAngleLProtonMuonPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{#mu}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{#mu}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSICosOpeningAngleRProtonMuonPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{R},#vec{p}_{#mu}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{R},#vec{p}_{#mu}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    // GKI
    PlotNames.push_back("TrueNoFSICosOpeningAngleMomentumTransferTotalProtonPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{q},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{q},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIMissingMomentumPlot");
    XAxisLabel.push_back("p_{n} [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueNoFSIAlphaThreeDPlot");
    XAxisLabel.push_back("#alpha_{3D} [deg]");
    YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    // Post FSI
    PlotNames.push_back("TrueVertexXPlot");
    XAxisLabel.push_back("#vec{v}_{x} [cm]");
    YAxisLabel.push_back("#frac{d#sigma}{d #vec{v}_{x}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueVertexYPlot");
    XAxisLabel.push_back("#vec{v}_{y} [cm]");
    YAxisLabel.push_back("#frac{d#sigma}{d #vec{v}_{y}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueVertexZPlot");
    XAxisLabel.push_back("#vec{v}_{z} [cm]");
    YAxisLabel.push_back("#frac{d#sigma}{d #vec{v}_{z}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueMuonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{#mu}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueLeadingProtonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L}}");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueRecoilProtonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{R}}");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueLeadingProtonMomentumPlot");
    XAxisLabel.push_back("|#vec{p}_{L}| [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{L}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueRecoilProtonMomentumPlot");
    XAxisLabel.push_back("|#vec{p}_{R}| [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{R}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueMuonMomentumPlot");
    XAxisLabel.push_back("|#vec{p}_{#mu}| [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d|#vec{p}_{#mu}|} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueCosOpeningAngleProtonsPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueCosOpeningAngleMuonTotalProtonPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueTransverseMomentumPlot");
    XAxisLabel.push_back("#delta P_{T} [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueDeltaAlphaTPlot");
    XAxisLabel.push_back("#delta #alpha_{T} [deg]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta #alpha_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
	
    PlotNames.push_back("TrueInvariantMassPlot");
    XAxisLabel.push_back("W [GeV]");
    YAxisLabel.push_back("#frac{d#sigma}{dW} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueCosOpeningAngleLProtonMuonPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{#mu}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{#mu}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueCosOpeningAngleRProtonMuonPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{R},#vec{p}_{#mu}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{R},#vec{p}_{#mu}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");
   
    // GKI
    PlotNames.push_back("TrueCosOpeningAngleMomentumTransferTotalProtonPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{q},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{q},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueMissingMomentumPlot");
    XAxisLabel.push_back("p_{n} [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueAlphaThreeDPlot");
    XAxisLabel.push_back("#alpha_{3D} [deg]");
    YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    // Double differential final state
    PlotNames.push_back("TrueSerialTransverseMomentum_InMuonCosThetaPlot");
    XAxisLabel.push_back("#delta P_{T} [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialDeltaAlphaT_InMuonCosThetaPlot");
    XAxisLabel.push_back("#delta #alpha_{T} [deg]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta #alpha_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    // GKI
    PlotNames.push_back("TrueSerialMissingMomentum_InMuonCosThetaPlot");
    XAxisLabel.push_back("p_{n} [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialAlphaThreeD_InMuonCosThetaPlot");
    XAxisLabel.push_back("#alpha_{3D} [deg]");
    YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialCosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{q},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{q},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    // Double differential pre FSI
    PlotNames.push_back("TrueSerialNoFSITransverseMomentum_InMuonCosThetaPlot");
    XAxisLabel.push_back("#delta P_{T} [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta P_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialNoFSIDeltaAlphaT_InMuonCosThetaPlot");
    XAxisLabel.push_back("#delta #alpha_{T} [deg]");
    YAxisLabel.push_back("#frac{d#sigma}{d#delta #alpha_{T}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialNoFSICosOpeningAngleProtons_InMuonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{L},#vec{p}_{R}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{L},#vec{p}_{R}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialNoFSICosOpeningAngleMuonTotalProton_InMuonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    // GKI
    PlotNames.push_back("TrueSerialNoFSIMissingMomentum_InMuonCosThetaPlot");
    XAxisLabel.push_back("p_{n} [GeV/c]");
    YAxisLabel.push_back("#frac{d#sigma}{dp_{n}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialNoFSIAlphaThreeD_InMuonCosThetaPlot");
    XAxisLabel.push_back("#alpha_{3D} [deg]");
    YAxisLabel.push_back("#frac{d#sigma}{d#alpha_{3D}} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");

    PlotNames.push_back("TrueSerialNoFSICosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot");
    XAxisLabel.push_back("cos(#theta_{#vec{q},#vec{p}_{sum}})");
    YAxisLabel.push_back("#frac{d#sigma}{dcos(#theta_{#vec{q},#vec{p}_{sum}})} #left[10^{-38} #frac{cm^{2}}{Ar}#right]");


    //------------------------------//

    const int NPlots = PlotNames.size();
    const int NXLabels = XAxisLabel.size();
    const int NYLabels = YAxisLabel.size();

    // sanity check
    if ((NPlots != NXLabels) || (NPlots != NYLabels)) { cout << "Inconsistent number of plots and labels! Aborting !" << endl; return; }

    //------------------------------//	

    // Loop over the samples to open the files and the TTree

    for (int iSample = 0; iSample < NSamples; iSample++) {
        Files[iSample] = new TFile(Names[iSample],"readonly");
    } // End of the loop over the samples

    //------------------------------//

    // Loop over the plots to be compared

    for (int iPlot = 0; iPlot < NPlots; iPlot++) {

        TString PlotNameDuplicate = PlotNames[iPlot];
        TString GeneralPlotName = PlotNameDuplicate.ReplaceAll("NoFSI","");

        for (int iSample = 0; iSample < NSamples; iSample++) {	
            if (PlotNames[iPlot].Contains("NoFSI") && Labels[iSample] == "GiBUU_NoFSI") {
                continue;
            }

            TString LabelCopy = Labels[iSample];
            TString CanvasName = "ThreeDKI_"+LabelCopy.ReplaceAll(" ","_")+"_InteBreakDown_" + PlotNames[iPlot];
            TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
            PlotCanvas->cd();
            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);		
            PlotCanvas->Draw();	

            TLegend* leg = new TLegend(0.13,0.88,0.95,0.99);
            leg->SetBorderSize(0);
            leg->SetNColumns(3);
            leg->SetTextSize(TextSize);	
            leg->SetTextFont(FontStyle);						
            leg->SetMargin(0.2);				
            leg->SetFillColor(0);				

            // Loop over the interaction processes

            std::vector<TH1D*> Histos; Histos.resize(NProcesses);

            for (int iProcess = 0; iProcess < NProcesses; iProcess++) {	
                Histos[iProcess] = (TH1D*)(Files[iSample]->Get(Process[iProcess]+PlotNames[iPlot]));

                Histos[iProcess]->SetLineWidth(4);
                Histos[iProcess]->SetLineColor( Colors.at(iProcess) );	

                Histos[iProcess]->GetXaxis()->SetTitleFont(FontStyle);
                Histos[iProcess]->GetXaxis()->SetLabelFont(FontStyle);
                Histos[iProcess]->GetXaxis()->SetNdivisions(8);
                Histos[iProcess]->GetXaxis()->SetLabelSize(TextSize);
                Histos[iProcess]->GetXaxis()->SetTitleSize(TextSize);	
                Histos[iProcess]->GetXaxis()->SetTitleOffset(1.1);					
                Histos[iProcess]->GetXaxis()->CenterTitle();						

                Histos[iProcess]->GetYaxis()->SetTitleFont(FontStyle);
                Histos[iProcess]->GetYaxis()->SetLabelFont(FontStyle);
                Histos[iProcess]->GetYaxis()->SetNdivisions(6);
                Histos[iProcess]->GetYaxis()->SetLabelSize(TextSize);
                Histos[iProcess]->GetYaxis()->SetTitle(YAxisLabel.at(iPlot));
                Histos[iProcess]->GetYaxis()->SetTitleSize(TextSize);
                Histos[iProcess]->GetYaxis()->SetTitleOffset(1.25);
                Histos[iProcess]->GetYaxis()->SetTickSize(0);
                Histos[iProcess]->GetYaxis()->CenterTitle();	
                Histos[iProcess]->GetYaxis()->SetRangeUser(0.,1.15*Histos[0]->GetMaximum());

                Histos[iProcess]->Draw("hist same");
                Histos[0]->Draw("hist same");

                double frac = Histos[iProcess]->Integral("width")/Histos[0]->Integral("width") * 100.;
                TString LegLabel = Process[iProcess] + " (" + tools.to_string_with_precision(frac,1) + "%)";
                if (iProcess == 0) { LegLabel = "Total (" + tools.to_string_with_precision(frac,1) + "%)"; }
                TLegendEntry* legColor = leg->AddEntry(Histos[iProcess],LegLabel,"l");
                legColor->SetTextColor( Colors.at(iProcess) ); 
            } // End of the loop over the processes
            PlotCanvas->cd();
            leg->Draw();

            TLatex *textSlice = new TLatex();
            textSlice->SetTextFont(FontStyle);
            textSlice->SetTextSize(TextSize);
            TString PlotNameDuplicate = PlotNames[iPlot];
            TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("True","") ;
            textSlice->DrawLatexNDC(0.2, 0.81, Labels[iSample] + "      " + LatexLabel[ReducedPlotName].ReplaceAll("All events",""));

            gPad->RedrawAxis();
            TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";
            TString SaveDirectory = (PlotNames[iPlot].Contains("NoFSI")) ? "PreFSI" : "PostFSI";
            PlotCanvas->SaveAs(dir+"/Figs/InteBreakDown/"+SaveDirectory+"/InteBreakDown_"+Labels[iSample]+"_"+PlotNames[iPlot]+".png");
            delete PlotCanvas;
        } // End of the loop over the samples grabing the plots	

        // If variable is double differential, also generate unserialized plots
        if (PlotNames[iPlot].Contains("Serial")) {
            std::vector<std::vector<std::vector<TH1D*>>> Histos;
            std::vector<std::vector<TH1D*>> TruePlots; TruePlots.resize(NSamples);
            for (int iSample = 0; iSample < NSamples; iSample++) TruePlots[iSample].resize(NProcesses);

            //------------------------------//
        
            // Flatten out double differential bins
            TString PlotNameDuplicate = PlotNames[iPlot]; PlotNameDuplicate.ReplaceAll("NoFSI", "");
            auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator[PlotNameDuplicate];
            auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
            
            int StartIndex = 0;

            //------------------------------//
            // Load true plots
            for (int iSample = 0; iSample < NSamples; iSample++) {
                for (int iProcess = 0; iProcess < NProcesses; iProcess++) {	
                    TString PlotName = Process[iProcess]+PlotNames[iPlot];
                    TruePlots[iSample][iProcess] = (TH1D*)(Files[iSample]->Get(PlotName));
                }
            }

            // Resize histos
            Histos.resize(NSlices);
            for (int iSlice = 0; iSlice < NSlices; iSlice++) {
                Histos[iSlice].resize(NSamples);
                for (int iSample = 0; iSample < NSamples; iSample++) {
                    Histos[iSlice][iSample].resize(NProcesses);
                }
            }

            //------------------------------//

            // Loop over the slices
            for (int iSlice = 0; iSlice < NSlices; iSlice++) {
                int SliceNBins = SerialVectorBins.at(iSlice);
                for (int iSample = 0; iSample < NSamples; iSample++) {
                    if (PlotNames[iPlot].Contains("NoFSI") && Labels[iSample] == "GiBUU_NoFSI") {
                        continue;
                    }

                    TString SlicePlotName = Labels[iSample] + "_" + PlotNames[iPlot] + "_" + TString(std::to_string(iSlice));
                    double SliceWidth = SliceDiscriminators[iSlice + 1] - SliceDiscriminators[iSlice]; 
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

                    TLegend* leg = new TLegend(0.13,0.88,0.95,0.99);
                    leg->SetBorderSize(0);
                    leg->SetNColumns(3);
                    leg->SetTextSize(TextSize);	
                    leg->SetTextFont(FontStyle);						
                    leg->SetMargin(0.2);				
                    leg->SetFillColor(0);	

                    for (int iProcess = 0; iProcess < NProcesses; iProcess++) {
                        Histos[iSlice][iSample][iProcess] = tools.GetHistoBins(
                            TruePlots[iSample][iProcess],
                            SerialVectorLowBin.at(iSlice),
                            SerialVectorHighBin.at(iSlice),
                            SliceWidth,
                            SerialSliceBinning,
                            Labels[iSample]
                        );
                        Histos[iSlice][iSample][iProcess]->SetLineWidth(4);
                        Histos[iSlice][iSample][iProcess]->SetLineColor(Colors.at(iProcess));

                        Histos[iSlice][iSample][iProcess]->GetXaxis()->SetTitleFont(FontStyle);
                        Histos[iSlice][iSample][iProcess]->GetXaxis()->SetLabelFont(FontStyle);
                        Histos[iSlice][iSample][iProcess]->GetXaxis()->SetNdivisions(8);
                        Histos[iSlice][iSample][iProcess]->GetXaxis()->SetLabelSize(TextSize);
                        Histos[iSlice][iSample][iProcess]->GetXaxis()->SetTitle(XAxisLabel.at(iPlot));
                        Histos[iSlice][iSample][iProcess]->GetXaxis()->SetTitleSize(TextSize);
                        Histos[iSlice][iSample][iProcess]->GetXaxis()->SetTitleOffset(1.1);
                        Histos[iSlice][iSample][iProcess]->GetXaxis()->CenterTitle();

                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetTitleFont(FontStyle);
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetLabelFont(FontStyle);
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetNdivisions(6);
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetLabelSize(TextSize);
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetTitle(YAxisLabel.at(iPlot));
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetTitleSize(TextSize);
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetTitleOffset(1.3);
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetTickSize(0);
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->CenterTitle();
                        Histos[iSlice][iSample][iProcess]->GetYaxis()->SetRangeUser(0.,1.15*Histos[iSlice][iSample][0]->GetMaximum());

                        Histos[iSlice][iSample][iProcess]->Draw("hist same");
                        Histos[iSlice][iSample][0]->Draw("hist same");

                        double frac = Histos[iSlice][iSample][iProcess]->Integral("width")/Histos[iSlice][iSample][0]->Integral("width") * 100.;
                        TString LegLabel = Process[iProcess] + " (" + tools.to_string_with_precision(frac,1) + "%)";
                        if (iProcess == 0) { LegLabel = "Total (" + tools.to_string_with_precision(frac,1) + "%)"; }
                        TLegendEntry* legColor = leg->AddEntry(Histos[iSlice][iSample][0],LegLabel,"l");
                        legColor->SetTextColor(Colors.at(iProcess)); 
                    }
                    PlotCanvas->cd();
                    leg->Draw();

                    TLatex *textSlice = new TLatex();
                    textSlice->SetTextFont(FontStyle);
                    textSlice->SetTextSize(TextSize);
                    TString ReducedPlotName = PlotNameDuplicate; ReducedPlotName.ReplaceAll("True","");
                    TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[iSlice], 1) + " < " + PlotNameToSliceLabel[PlotNameDuplicate] + " < " + tools.to_string_with_precision(SliceDiscriminators[iSlice + 1], 1);
                    textSlice->DrawLatexNDC(0.2, 0.81, Labels[iSample] + "      " + LatexLabel[ReducedPlotName].ReplaceAll("All events","") + SliceLabel);

                    gPad->RedrawAxis();
                    TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";
                    TString SaveDirectory = (PlotNames[iPlot].Contains("NoFSI")) ? "PreFSI" : "PostFSI";
                    PlotCanvas->SaveAs(dir+"/Figs/InteBreakDown/"+SaveDirectory+"/InteBreakDown_"+SlicePlotName+".png");
                    delete PlotCanvas;

                }
                // Update start index for next slice
                StartIndex += (SliceNBins + 1);
            }
        }

    } // End of the loop over the plots
} // End of the program
