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
#include "THStack.h"
#include "TPad.h"

// std includes.
#include <vector>
#include <memory>

// Definitions for CutVars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"
#include "../../Utils/Tools.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionCutPlots() {

    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    //SpectrumLoader NuLoader(mc_wc);
    //SpectrumLoader DataNuLoader(data_wc);
    //SpectrumLoader CosmicNuLoader(cosmic_wc);

//    SpectrumLoader NuLoader(InputFiles);
//    SpectrumLoader DataNuLoader(DataInputFiles); 
//    SpectrumLoader CosmicNuLoader(CosmicInputFiles);    

    SpectrumLoader NuLoader("/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025Av3/v10_04_06_01/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/00/02/caf-46414215-0544-42a7-b617-88835bbd7c30.root");
    SpectrumLoader DataNuLoader("/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025Av3/v10_04_06_01/MCP2025Av3_DevSample/flatcaf/bnblight/1d/reco2_reco1_filtered_decoded-raw_filtered_data_EventBuilder2_art2_run18255_64_strmBNBLight_20250218T121038-1da89634-988e-9dce-f4bb-b0336e22b333.flat.caf.root");
    SpectrumLoader CosmicNuLoader("/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025Av3/v10_04_06_01/prodcorsika_proton_intime_sbnd/CV/caf/da/68/caf.flat.caf-7c765d51-47f1-485c-aded-244f5d04a2c1.root"); 

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionCutPlots.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    // Colors to use
    std::vector<int> Colors;
    Colors.push_back(kAzure+7);
    Colors.push_back(kOrange-3);
    Colors.push_back(kGreen+1);
    Colors.push_back(kBlue-4);

    // Cuts and CutVars to plot
    std::vector<std::vector<Var>> CutVars; std::vector<Binning> CutVarBins; std::vector<Cut> SignalCuts;
    std::vector<TString> PlotNames; std::vector<std::string> CutVarLabels;
    std::vector<std::vector<Cut>> HistosCuts; std::vector<std::vector<TString>> HistosLabels;

    // Nu score before any cut
    CutVars.push_back({kNuScore}); CutVarBins.push_back(bNuScore); SignalCuts.push_back(kNoCut);
    PlotNames.push_back("NoCutNuScore"); CutVarLabels.push_back("#nu score");
    HistosCuts.push_back({kIsNotCosmic, kIsCosmic}); HistosLabels.push_back({"MC non-cosmic", "MC cosmic"});

    // Fmatch score before any cut
 /*   CutVars.push_back({kFMatchScore}); CutVarBins.push_back(bFMatchScore); SignalCuts.push_back(kNoCut);
    PlotNames.push_back("NoCutFMatchScore"); CutVarLabels.push_back("flash matching score");
    HistosCuts.push_back({kIsNotCosmic, kIsCosmic}); HistosLabels.push_back({"Not cosmic", "Cosmic"});

    // Fmatch time before any cut
    CutVars.push_back({kFMatchTime}); CutVarBins.push_back(bFMatchTime); SignalCuts.push_back(kNoCut);
    PlotNames.push_back("NoCutFMatchTime"); CutVarLabels.push_back("flash matching time");
    HistosCuts.push_back({kIsNotCosmic, kIsCosmic}); HistosLabels.push_back({"Not cosmic", "Cosmic"});*/

    // Muon chi square for muon cut
    CutVars.push_back({kMuMuChi2, kProtonMuChi2, kSecondProtonMuChi2, kPionMuChi2}); CutVarBins.push_back(bMuChi2); SignalCuts.push_back(kCosmicCut);
    PlotNames.push_back("ParticleCutChi2Muon"); CutVarLabels.push_back("#chi^{2}_{#mu}"); 
    HistosCuts.push_back({kHasMuon, kHasProton, kHasSecondProton, kHasPion}); HistosLabels.push_back({"#mu", "p1", "p2", "#pi"});

    // Proton chi square for muon cut
    CutVars.push_back({kMuProtonChi2, kProtonProtonChi2, kSecondProtonProtonChi2, kPionProtonChi2}); CutVarBins.push_back(bProtonChi2); SignalCuts.push_back(kCosmicCut);
    PlotNames.push_back("ParticleCutChi2Proton"); CutVarLabels.push_back("#chi^{2}_{p}"); 
    HistosCuts.push_back({kHasMuon, kHasProton, kHasSecondProton, kHasPion}); HistosLabels.push_back({"#mu", "p1", "p2", "#pi"});

    //------------------------------------//  

    // Construct all spectra
    // For MC, we need the breakdown of the contributions
    // Aka we need a vector of vectors
    std::vector< std::vector<std::unique_ptr<Spectrum> > > Spectra;
    // That is not the case for data
    // where we don't have the contributions
    std::vector< std::unique_ptr<Spectrum> > DataSpectra;  
    std::vector< std::unique_ptr<Spectrum> > CosmicSpectra;        

    //------------------------------------//  

    for (std::size_t i = 0; i < CutVars.size(); i++) {

        // MC

        std::vector<std::unique_ptr<Spectrum>> InnerSpectra;

        for (std::size_t j = 0; j < HistosCuts[i].size(); j++) {

            std::unique_ptr<Spectrum> Signals;            

            if (CutVars[i].size() == 1) {
 
                Signals = std::make_unique<Spectrum>(CutVarLabels[i], CutVarBins[i], NuLoader, CutVars[i][0], kNoSpillCut, SignalCuts[i] && HistosCuts[i][j]);
 
            } else {
 
                Signals = std::make_unique<Spectrum>(CutVarLabels[i], CutVarBins[i], NuLoader, CutVars[i][j], kNoSpillCut, SignalCuts[i] && HistosCuts[i][j]);

            }
 
            InnerSpectra.push_back(std::move(Signals));
 
        }
 
        Spectra.push_back(std::move(InnerSpectra));

        //------------------------------------//

        // Data

        std::unique_ptr<Spectrum> DataSignals;
        DataSignals = std::make_unique<Spectrum>(CutVarLabels[i], CutVarBins[i], DataNuLoader, CutVars[i][0], kNoSpillCut, SignalCuts[i]);
        DataSpectra.push_back(std::move(DataSignals)); 
        
        //------------------------------------//

        // Cosmics

        std::unique_ptr<Spectrum> CosmicSignals;
        CosmicSignals = std::make_unique<Spectrum>(CutVarLabels[i], CutVarBins[i], CosmicNuLoader, CutVars[i][0], kNoSpillCut, SignalCuts[i]);
        CosmicSpectra.push_back(std::move(CosmicSignals)); 
        
        //------------------------------------//        
 
    }

    //------------------------------------//      

    NuLoader.Go();
    DataNuLoader.Go();  
    CosmicNuLoader.Go();        

    //------------------------------------//   

    // Create all the vectors for all the objects needed 
      
    vector<TCanvas*> PlotCanvas; PlotCanvas.clear(); PlotCanvas.resize( (int)CutVars.size() );
    vector<TLegend*> leg; leg.clear(); leg.resize( (int)CutVars.size() ); 
    vector<TH1D*> DataHistos; DataHistos.clear(); DataHistos.resize( (int)CutVars.size() );  
    vector<TH1D*> CosmicHistos; CosmicHistos.clear(); CosmicHistos.resize( (int)CutVars.size() );  
    vector< vector<TH1D*> > Histos; Histos.clear(); Histos.resize( (int)CutVars.size() ); 
    vector<THStack*> stacked_histos; stacked_histos.clear(); stacked_histos.resize( (int)CutVars.size() );

    //------------------------------------// 

    // loop over the selection cuts

    for (std::size_t i = 0; i < CutVars.size(); i++) {

        //------------------------------------//         

        // Canvas 

        PlotCanvas[i] = new TCanvas("Selection_"+PlotNames[i],"Selection_"+PlotNames[i],205,34,1124,768);

        PlotCanvas[i]->SetTopMargin(0.13);
        PlotCanvas[i]->SetLeftMargin(0.15);
        PlotCanvas[i]->SetRightMargin(0.05);
        PlotCanvas[i]->SetBottomMargin(0.14);

        //------------------------------------//  
        
        // Legend

        leg[i] = new TLegend(0.22,0.88,0.98,0.98);
        leg[i]->SetBorderSize(0);
        leg[i]->SetNColumns(2);
        if (PlotNames[i] != "NoCutNuScore") { leg[i]->SetNColumns(3); }
        leg[i]->SetMargin(0.05);
        leg[i]->SetTextSize(TextSize);
        leg[i]->SetTextFont(FontStyle);

        //------------------------------------//      
        
        // Data plot
        
        auto& DataSignals = DataSpectra[i];
        DataHistos[i] = DataSignals->ToTH1(TargetPOT);  
        
        DataHistos[i]->SetLineColor(kBlack);
        DataHistos[i]->SetMarkerColor(kBlack);
        DataHistos[i]->SetMarkerStyle(20);            
        DataHistos[i]->SetMarkerSize(2.);    

        SaveFile->WriteObject(DataHistos[i], PlotNames[i]+"_data");

        //------------------------------------//      
        
        // Cosmic plot
        
        auto& CosmicSignals = CosmicSpectra[i];
        CosmicHistos[i] = CosmicSignals->ToTH1(TargetPOT);  
        CosmicHistos[i]->Scale(TargetPOT/cosmic_pot);
        
        CosmicHistos[i]->SetLineColor(kGray);
        CosmicHistos[i]->SetMarkerColor(kGray);
        CosmicHistos[i]->SetFillColor(kGray);  
        CosmicHistos[i]->SetFillStyle(3244);   
        
        SaveFile->WriteObject(CosmicHistos[i], PlotNames[i]+"_cosmic");        
        
        //------------------------------------//          
        
        Histos[i].resize(HistosCuts[i].size());
        stacked_histos[i] = new THStack(PlotNames[i],"");
        stacked_histos[i]->Add(CosmicHistos[i],"hist");

        //------------------------------------//  
        
        // loop over the different breakdown contributions given a cut

        for (std::size_t j = 0; j < HistosCuts[i].size(); j++) {

            auto& Signals = Spectra[i][j];
            Histos[i][j] = Signals->ToTH1(TargetPOT);    
            
            Histos[i][j]->SetLineColor(Colors[j]);
            Histos[i][j]->SetFillColor(Colors[j]); 
            Histos[i][j]->SetLineWidth(4);

            Histos[i][j]->GetXaxis()->SetTitleFont(FontStyle);
            Histos[i][j]->GetXaxis()->SetLabelFont(FontStyle);
            Histos[i][j]->GetXaxis()->SetNdivisions(10);
            Histos[i][j]->GetXaxis()->SetLabelSize(TextSize);
            Histos[i][j]->GetXaxis()->SetTitleSize(TextSize);
            Histos[i][j]->GetXaxis()->SetTitleOffset(1.1);
            Histos[i][j]->GetXaxis()->CenterTitle();
            Histos[i][j]->GetXaxis()->SetTitle((CutVarLabels[i]).c_str());

            Histos[i][j]->GetYaxis()->SetTitleFont(FontStyle);
            Histos[i][j]->GetYaxis()->SetLabelFont(FontStyle);
            Histos[i][j]->GetYaxis()->SetNdivisions(6);
            Histos[i][j]->GetYaxis()->SetLabelSize(TextSize);
            Histos[i][j]->GetYaxis()->SetTitleSize(TextSize);
            Histos[i][j]->GetYaxis()->SetTitleOffset(1.3);
            Histos[i][j]->GetYaxis()->CenterTitle();

            stacked_histos[i]->Add(Histos[i][j],"hist");
            stacked_histos[i]->Draw("hist");

            // Save to root file

            SaveFile->WriteObject(Histos[i][j], PlotNames[i]+HistosLabels[i][j]+"_reco");

        } // end of the loop over the different breakdown contributions given a cut

        //------------------------------------//

        Tools tools; 

        TString data_events = tools.to_string_with_precision(DataHistos[i]->Integral(),0);
        TLegendEntry* legData = leg[i]->AddEntry(DataHistos[i],"Data [" + data_events+"]","e1x0p");

        TH1D* total_mc_cosmic = (TH1D*) (stacked_histos[i]->GetStack()->Last());

	    TString cosmic_events = tools.to_string_with_precision(CosmicHistos[i]->Integral(),1);
        TString cosmic_frac = tools.to_string_with_precision(CosmicHistos[i]->Integral()/total_mc_cosmic->Integral() * 100.,1);
        leg[i]->AddEntry(CosmicHistos[i],"cosmics (" + cosmic_frac + "%) [" + cosmic_events +"]","f");
        
	    //------------------------------------//

        double imax = TMath::Max(DataHistos[i]->GetMaximum(), total_mc_cosmic->GetMaximum());
        double YAxisRange = 1.15*imax;   

        total_mc_cosmic->GetYaxis()->SetRangeUser(0,YAxisRange);
        total_mc_cosmic->SetLineColor(kAzure + 7);
        total_mc_cosmic->SetLineWidth(1);        
        PlotCanvas[i]->cd();
        total_mc_cosmic->Draw("hist");

       for (int j = (int)(HistosCuts[i].size())-1; j >= 0; j--) {  

	        TString frac = tools.to_string_with_precision(Histos[i][j]->Integral()/total_mc_cosmic->Integral() * 100.,1); 

            leg[i]->AddEntry(Histos[i][j],HistosLabels[i][j] + " (" + frac + "%)","f");

        }

	    stacked_histos[i]->Draw("same hist");
        DataHistos[i]->Draw("e1x0 same"); 

        leg[i]->Draw();

	    gPad->RedrawAxis();

        // Save as pdf
        PlotCanvas[i]->SaveAs(dir_figs+"/Figs/CAFAna/Cuts/"+PlotNames[i]+".pdf");

        delete PlotCanvas[i];

    } // end of the loop over the selection cuts

    // Close file
    SaveFile->Close();

}
