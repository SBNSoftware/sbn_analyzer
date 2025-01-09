#include <TFile.h>
#include <TTree.h>
#include <TString.h>

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

void GeneratorOverlay() {

    //------------------------------//

    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    int FontStyle = 132;
    double TextSize = 0.06;			

    TString OutFilePath = "/pnfs/sbnd/persistent/users/" + (TString)UserName + "/HighSamples/FlatTree/";

    //------------------------------//

    Tools tools;

    //------------------------------//

    // Event generators

    std::vector<TString> Names; std::vector<TString> Labels; std::vector<int> Colors;
    
    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_AR23.root"); 
    Labels.push_back("GENIE AR23");
    Colors.push_back(kBlue+8);

    // Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_AR23_Emp2015.root"); 
    // Labels.push_back("GENIE AR23 Emp2015");
    // Colors.push_back(kBlue+16);

    // Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_AR23_Nieves2016.root"); 
    // Labels.push_back("GENIE AR23 Nieves2016");
    // Colors.push_back(kBlue+24);

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GENIE_G18.root"); 
    Labels.push_back("GENIE G18");
    Colors.push_back(kBlue+2);

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NuWro.root"); 
    Labels.push_back("NuWro");
    Colors.push_back(kRed+1);

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_NEUT.root"); 
    Labels.push_back("NEUT");
    Colors.push_back(kOrange+7);

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU.root"); 
    Labels.push_back("GiBUU");
    Colors.push_back(kGreen+1);

    Names.push_back(OutFilePath+"FlatTreeAnalyzerOutput_GiBUU_NoFSI.root"); 
    Labels.push_back("GiBUU NoFSI");
    Colors.push_back(kGreen+1);

    const int NSamples = Names.size();
    std::vector<TFile*> Files; Files.resize(NSamples);

    //------------------------------//

    // Plots to overlay

    std::vector<TString> PlotNames;

    // Pre FSI
    PlotNames.push_back("TrueNoFSILeadingProtonCosThetaPlot");
    PlotNames.push_back("TrueNoFSIRecoilProtonCosThetaPlot");
    PlotNames.push_back("TrueNoFSILeadingProtonMomentumPlot");
    PlotNames.push_back("TrueNoFSIRecoilProtonMomentumPlot");
    PlotNames.push_back("TrueNoFSIMuonMomentumPlot");
    PlotNames.push_back("TrueNoFSICosOpeningAngleProtonsPlot");
    PlotNames.push_back("TrueNoFSICosOpeningAngleMuonTotalProtonPlot");
    PlotNames.push_back("TrueNoFSITransverseMomentumPlot");
    PlotNames.push_back("TrueNoFSIDeltaAlphaTPlot");
    PlotNames.push_back("TrueNoFSICosOpeningAngleMomentumTransferTotalProtonPlot");
    PlotNames.push_back("TrueNoFSIMissingMomentumPlot");
    PlotNames.push_back("TrueNoFSIAlphaThreeDPlot");
    PlotNames.push_back("TrueNoFSIInvariantMassPlot");
    PlotNames.push_back("TrueNoFSICosOpeningAngleLProtonMuonPlot"); 
    PlotNames.push_back("TrueNoFSICosOpeningAngleRProtonMuonPlot");      
    
    // Post FSI
    PlotNames.push_back("TrueVertexXPlot");
    PlotNames.push_back("TrueVertexYPlot");
    PlotNames.push_back("TrueVertexZPlot");
    PlotNames.push_back("TrueMuonCosThetaPlot");
    PlotNames.push_back("TrueLeadingProtonCosThetaPlot");
    PlotNames.push_back("TrueRecoilProtonCosThetaPlot");
    PlotNames.push_back("TrueLeadingProtonMomentumPlot");
    PlotNames.push_back("TrueRecoilProtonMomentumPlot");
    PlotNames.push_back("TrueMuonMomentumPlot");
    PlotNames.push_back("TrueCosOpeningAngleProtonsPlot");
    PlotNames.push_back("TrueCosOpeningAngleMuonTotalProtonPlot");
    PlotNames.push_back("TrueTransverseMomentumPlot");
    PlotNames.push_back("TrueDeltaAlphaTPlot");
    PlotNames.push_back("TrueCosOpeningAngleMomentumTransferTotalProtonPlot");
    PlotNames.push_back("TrueMissingMomentumPlot");
    PlotNames.push_back("TrueAlphaThreeDPlot");
    PlotNames.push_back("TrueInvariantMassPlot");
    PlotNames.push_back("TrueCosOpeningAngleLProtonMuonPlot");   
    PlotNames.push_back("TrueCosOpeningAngleRProtonMuonPlot");

    // Double differential final state
    PlotNames.push_back("TrueSerialTransverseMomentum_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialDeltaAlphaT_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialMissingMomentum_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialAlphaThreeD_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialCosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot");

    // Double differential pre FSI
    PlotNames.push_back("TrueSerialNoFSITransverseMomentum_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialNoFSIDeltaAlphaT_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialNoFSICosOpeningAngleProtons_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialNoFSICosOpeningAngleMuonTotalProton_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialNoFSIMissingMomentum_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialNoFSIAlphaThreeD_InMuonCosThetaPlot");
    PlotNames.push_back("TrueSerialNoFSICosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot");

    const int NPlots = PlotNames.size();

    //------------------------------//	

    // Loop over the samples to open the files and the TTree

    for (int iSample = 0; iSample < NSamples; iSample++) {

        Files[iSample] = new TFile(Names[iSample],"readonly");

    } // End of the loop over the samples

    //------------------------------//

    // Loop over the plots to be compared

    for (int iPlot = 0; iPlot < NPlots; iPlot++) {

        TString CanvasName = "Canvas_" + PlotNames[iPlot];
        TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
        PlotCanvas->cd();
        PlotCanvas->SetTopMargin(0.12);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetBottomMargin(0.15);		
        PlotCanvas->Draw();	

        TLegend* leg = new TLegend(0.2,0.7,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(2);
        leg->SetTextSize(TextSize*0.8);	
        leg->SetTextFont(FontStyle);						

        // Loop over the samples to open the files and to get the corresponding plot

        std::vector<TH1D*> Histos; Histos.resize(NSamples);

        for (int iSample = 0; iSample < NSamples; iSample++) {	
            // Exclude normal GiBUU in pre-FSI plots and exclude pre-FSI plots in final state plots
            if (PlotNames[iPlot].Contains("NoFSI") && Labels[iSample]=="GiBUU") { continue; } 
            if (!PlotNames[iPlot].Contains("NoFSI") && Labels[iSample]=="GiBUU NoFSI") { continue; }
            
            // For pre-FSI plots, use GiBUU NoFSI's final state variables
            TString PlotName = PlotNames[iPlot];
            if (PlotNames[iPlot].Contains("NoFSI") && Labels[iSample]=="GiBUU NoFSI") { PlotName.ReplaceAll("NoFSI",""); }
            Histos[iSample] = (TH1D*)(Files[iSample]->Get(PlotName));

            Histos[iSample]->SetLineWidth(4);
            Histos[iSample]->SetLineColor( Colors.at(iSample) );	

            Histos[iSample]->GetXaxis()->SetTitleFont(FontStyle);
            Histos[iSample]->GetXaxis()->SetLabelFont(FontStyle);
            Histos[iSample]->GetXaxis()->SetNdivisions(8);
            Histos[iSample]->GetXaxis()->SetLabelSize(TextSize);
            Histos[iSample]->GetXaxis()->SetTitleSize(TextSize);	
            Histos[iSample]->GetXaxis()->SetTitleOffset(1.1);					
            Histos[iSample]->GetXaxis()->CenterTitle();						

            Histos[iSample]->GetYaxis()->SetTitleFont(FontStyle);
            Histos[iSample]->GetYaxis()->SetLabelFont(FontStyle);
            Histos[iSample]->GetYaxis()->SetNdivisions(6);
            Histos[iSample]->GetYaxis()->SetLabelSize(TextSize);
            Histos[iSample]->GetYaxis()->SetTitle("Cross Section [10^{-38} cm^{2}/Ar]");
            Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
            Histos[iSample]->GetYaxis()->SetTitleOffset(1.3);
            Histos[iSample]->GetYaxis()->SetTickSize(0);
            Histos[iSample]->GetYaxis()->CenterTitle();	

            double imax = TMath::Max(Histos[iSample]->GetMaximum(),Histos[0]->GetMaximum());

            double YAxisRange = 1.15*imax;
            if (PlotNames[iPlot] == "TrueLeadingProtonMomentumPlot") { YAxisRange *= 1.05; };
            Histos[iSample]->GetYaxis()->SetRangeUser(0.,YAxisRange);
            Histos[0]->GetYaxis()->SetRangeUser(0.,YAxisRange);			

            PlotCanvas->cd();
            Histos[iSample]->Draw("hist same");
            Histos[0]->Draw("hist same");	

            leg->AddEntry(Histos[iSample],Labels[iSample],"l");
            
            //----------------------------------------//					

        } // End of the loop over the samples grabing the plots	

        PlotCanvas->cd();
        leg->Draw();
        
        TLatex *textSlice = new TLatex();
        textSlice->SetTextFont(FontStyle);
        textSlice->SetTextSize(TextSize);
        TString PlotNameDuplicate = PlotNames[iPlot];
        TString GeneralPlotName = PlotNameDuplicate.ReplaceAll("NoFSI","");
        TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("True","");

        if (GeneralPlotName.Contains("TrueSerial")) {
            auto [SliceDiscriminators, SliceBinning] = PlotNameToDiscriminator[GeneralPlotName];
            auto [NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin] = tools.FlattenNDBins(SliceDiscriminators, SliceBinning);
            TString SliceLabel = tools.to_string_with_precision(SliceDiscriminators[0], 1) + " < " + PlotNameToSliceLabel[GeneralPlotName] + " < " + tools.to_string_with_precision(SliceDiscriminators[NSlices], 1);
            textSlice->DrawLatexNDC(0.16, 0.93, SliceLabel);
        }

        TString dir = "/exp/sbnd/app/users/" + (TString)UserName + "/CC1muAnalysis";
        TString SaveDirectory = (PlotNames[iPlot].Contains("NoFSI")) ? "PreFSI" : "PostFSI";
        PlotCanvas->SaveAs(dir+"/Figs/Overlay/"+SaveDirectory+"/Overlay_"+PlotNames[iPlot]+".png");
        delete PlotCanvas;

    } // End of the loop over the plots

    //------------------------------//

} // End of the program
