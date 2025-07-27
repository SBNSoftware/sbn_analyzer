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
#include "Utils.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void Selection_arg(bool doSIS=true){//, bool do2025B=true) {
  signalIsSIS=doSIS; //true; //defines what TrueSignal is in Defnitions.h

  bool debug=false;//true;
  bool runData=true;
  //signalIsSis (set it defs overridable here)
  // Set defaults and load tools
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  if(debug) cout<<"l34"<<endl;
  int FontStyle = 132;
  double TextSize = 0.06;	
  if(debug)cout<<"l37"<<endl;
  std::filesystem::create_directory(plotdir.Data() );
  // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
  //SpectrumLoader NuLoader("/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025B_5e18_02/v10_06_00_02/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/*/*/caf.flat*.root")
    //spring25 2025B 5e18POT
    //"/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025Av3/v10_04_06_01/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/*/*/*caf.flat.caf*.root");//2025A in a way thats big hopefully??? 
  //InputFiles);
  //SpectrumLoader NuLoader("/pnfs/sbn/data_add/sbn_nd/poms_production/official/MCP2024B/v09_91_02_02/prodoverlay_corsika_cosmics_proton_genie_rockbox_sce/caf/*/*/caf.flat.caf*.root");//2024B MC //InputFiles);
  //SpectrumLoader DataLoader("/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025B/v10_06_00/DevSample/flatcaf/bnblight/*/*/*flat.caf.root");
  //InputFilesData);
  //hack needs to get fixed just checking if these files are in fact the issue here. 
  // SpectrumLoader DataLoader("/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025Av3/v10_04_06_01/MCP2025Av3_DevSample/flatcaf/bnblight/00/*.flat.caf.root");//spring25 validations sample 
  string infile;
  vector<string> dnames, mcnames;

  std::ifstream mcinfiles("list_mc_spring25_5e18POT_2025B.txt"); //list_Data_MCP2025v3_DevSample.txt");

  while (std::getline(mcinfiles, infile)){
    mcnames.push_back(infile);
  }
  SpectrumLoader NuLoader(mcnames);
  

  std::ifstream datainfiles("list_data_spring25FixedDev_5e18POT_2025B.txt"); //list_Data_MCP2025v3_DevSample.txt");

  while (std::getline(datainfiles, infile)){
    dnames.push_back(infile);
  }
  SpectrumLoader DataLoader(dnames);


  //"/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025Av3/v10_04_06_01/MCP2025Av3_DevSample/flatcaf/bnblight/cd/reco2_reco1_filtered_decoded-raw_filtered_data_EventBuilder3_art2_run18255_54_strmBNBLight_20250218T062324-cd4bf471-e5e4-c724-d453-1c3f25c3fe35.flat.caf.root");//signle file in sprint25 validations sample 
  //"/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025Av3/v10_04_06_01/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/*/*/*caf.flat.caf*.root");//2025A in a way thats big hopefully??? 
  //"/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025Av3/v10_04_06_01/prodgenie_corsika_proton_rockbox_sbnd/CV/caf/*/*/*caf.flat.caf*.root");//MC trying a hack
    
  if(debug)cout<<"l40"<<endl;
  // We now create overlaid plots for several reconstructed variables and three lines:
  //     1. all selected reconstructed events
  //     2. reco signal events
  //     3. reco background events
    
  // Root file to store objects in
  TString dirPath = "/exp/sbnd/data/users/" + (TString)UserName + Form("/CAFAnaOutput/%s", tag.c_str());
  std::filesystem::create_directory(dirPath.Data() );
  TString RootFilePath=dirPath;
  RootFilePath+= signalIsSIS? "/SelectionSIS.root":"/SelectionCC.root";
  TFile* SaveFile = new TFile(RootFilePath, "RECREATE");
  if(debug)  cout<<"l49"<<endl;
  // Construct all spectra
  std::vector<std::unique_ptr<Spectrum>  > dataSpectra;
  std::vector<std::tuple<
                std::unique_ptr<Spectrum>,
                std::unique_ptr<Spectrum>,
                std::unique_ptr<Spectrum>
                >> Spectra;

  std::vector<std::vector< std::unique_ptr<Spectrum> > > channelSpectra, dataCutsSpectra, mcCutsSpectra;
  
  static const std::vector<std::tuple<Var, Var, TruthVar>> Vars=GetSISVars();
  static const std::vector<Binning> VarBins = GetSISBins();
  static const std::vector<std::string>VarLabels= GetSISVarLabels();
  std::vector<std::string> VarNames= GetSISVarNames();
  

  int ncuts=kRecoCuts.size();
  vector<string> cutnames=getRecoCutsShortNames();

  cout<<"vars:"<<VarsOld.size()<<"\t VarLabels"<<VarLabels.size()<<"\t VarBins"<<VarBins.size()<<endl;
    
  std::vector<Cut> channelsRecoSignal=getRecoTrueChannelCuts();
  std::vector<string> chanNames=getChannelShortNames();

  for (std::size_t i = 0; i < Vars.size(); i++) {
    //======MC spectra=============
    cout<<"l57 i:"<<i<<"  VarBins.NBins(): "<<VarBins[i].NBins()<<" ["<<VarBins[i].Min()<<","<<VarBins[i].Max()<<"]"<<endl;
	auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsSignal); 
	auto RecoTrueSignals = std::make_unique<Spectrum> (VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsTrueReco); 
	auto RecoBkgSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsBackground); 
	Spectra.push_back({std::move(RecoSignals), std::move(RecoTrueSignals), std::move(RecoBkgSignals)});

    //========Data spectra============
    if(runData){
      cout<<"Going to make Data spectra"<<endl;
      auto RecoData=std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), DataLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsSignal); 
      dataSpectra.push_back(std::move(RecoData));
      dataCutsSpectra.push_back({});//create empty vector for this variable.
      for(int c=0; c<ncuts; c++){
        auto DataCut=std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), DataLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoCuts[c]);
        dataCutsSpectra[i].push_back(std::move(DataCut)); 
      }

    }

    //=========cuts spectra for MC=========
     if(runData){
      mcCutsSpectra.push_back({});//create empty vector for this variable.
      for(int c=0; c<ncuts; c++){
        auto mcCut=std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoCuts[c]);
        mcCutsSpectra[i].push_back(std::move(mcCut)); 
      }

    }
    // loop over channels 
    channelSpectra.push_back({});//create an empty vector for this variable
    for(std::size_t c=0; c<channelsRecoSignal.size(); c++){
      auto RecoTrueChannelSignal = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, channelsRecoSignal[c]); 
      channelSpectra[i].push_back(std::move(RecoTrueChannelSignal)); 
      //channelSpectra[iVar][jChannel]
    }
	if(debug)cout<<"l110 i:"<<i<<endl;
  }
  if(debug)cout<<"l112"<<endl;
  // We now create spectra that will help us get the efficiency and purity data for each of the cuts

  // Spectrum with all events
  /*    Spectrum sAllEvents("AllEvents", bEventCount, NuLoader, kTrueEventCount, kNoTruthCut, kNoSpillCut);
  // Spectrum with all events that were reconstructed
  Spectrum sAllRecoEvents("AllRecoEvents", bEventCount, NuLoader, kTrueEventCount, kNoTruthCut, kNoSpillCut, kNoCut);
  // Spectrum with all true signal events
  Spectrum sAllTrueEvents("AllTrueEvents", bEventCount, NuLoader, kTrueEventCount, kTruthIsSignal, kNoSpillCut);
  // Spectrum with all true signal events that were reconstructed
  Spectrum sAllTrueRecoEvents("AllTrueRecoEvents", bEventCount, NuLoader, kTrueEventCount, kTruthIsSignal, kNoSpillCut, kNoCut);
  // Spectrum with first cut (cosmic)
  Spectrum sFirstCut("FirstCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFirstCut);
  Spectrum sFirstCutTrue("FirstCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFirstCutTrue);
  // Spectrum with second cut (cosmic and vertex FV)
  Spectrum sSecondCut("SecondCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kSecondCut);
  Spectrum sSecondCutTrue("SecondCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kSecondCutTrue);
  // Spectrum with second cut (cosmic, vertex FV, and one muon)
  Spectrum sThirdCut("ThirdCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kThirdCut);
  Spectrum sThirdCutTrue("ThirdCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kThirdCutTrue);
  // Spectrum with second cut (cosmic, vertex FV, one muon, and two protons)
  Spectrum sFourthCut("FourthCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFourthCut);
  Spectrum sFourthCutTrue("FourthCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFourthCutTrue);
  // Spectrum with second cut (cosmic, vertex FV, one muon, two protons, and no charged pions)
  Spectrum sFifthCut("FifthCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFifthCut);
  Spectrum sFifthCutTrue("FifthCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFifthCutTrue);
  // Spectrum with second cut (cosmic, vertex FV, one muon, two protons, no charged pions, and no neutral pions)
  Spectrum sSixthCut("SixthCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kSixthCut);
  Spectrum sSixthCutTrue("SixthCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kSixthCutTrue);
  // Spectrum with overall signal definition to sanity check it matches
  Spectrum sRecoSignal("RecoSignal", bEventCount, NuLoader, kEventCount, kNoSpillCut, kRecoIsSignal); 
  Spectrum sRecoTrueSignal("RecoTrueSignal", bEventCount, NuLoader, kEventCount, kNoSpillCut, kRecoIsTrueReco); */
  if(debug) cout<<"l144"<<endl;  
  NuLoader.Go();
  if(runData){
    cout<<"DataLoader.Go"<<endl;
    DataLoader.Go();
  }
  if(debug) cout<<"l150"<<endl;  
  // Loop over variables
  double mcPOT, dataPOT;
  if(debug) cout<<"l153"<<endl;  

  for (std::size_t i = 0; i < Vars.size(); i++) {
    if(debug) cout<<"in loop l156, i="<<i<<endl;  
    string var=VarNames[i];
    auto& [RecoSignals, RecoTrueSignals, RecoBkgSignals] = Spectra.at(i);
    
	if(i==0) mcPOT=RecoSignals->POT();
    /*RecoSignals->OverridePOT(1);
    RecoTrueSignals->OverridePOT(1);
    RecoBkgSignals->OverridePOT(1);
    */

	TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
	TH1D* RecoHisto = RecoSignals->ToTH1(TargetPOT); //setting targetpot to 5e18 ~the data pot 
	TH1D* RecoTrueHisto = RecoTrueSignals->ToTH1(TargetPOT);
	TH1D* RecoBkgHisto = RecoBkgSignals->ToTH1(TargetPOT);
	
    auto& DataSignals=dataSpectra.at(i);
    TH1D* DataHisto;
    vector<TH1D*> DataCutHistos, mcCutHistos;
    if(debug) cout<<"in loop l169, i="<<i<<endl;  
    if(runData){
      //if(i==0) dataPOT=DataSignals->POT();
      dataPOT=DataSignals->POT();
      if(debug) cout<<"making Data histo"<<endl;
      DataSignals->OverridePOT(1);
      //double nspills=DataSignals->Livetime();
      DataHisto= DataSignals->ToTH1(1);//nspills, kLivetime );//TargetPOT);

      for(int c=0; c<ncuts; c++){
        string cutname=cutnames[c];
        //auto dataCutSpectrum=dataCutsSpectra[i][c];
        //TH1D* chanHisto = channelSpectra[i][jChannel]->ToTH1(TargetPOT);
        dataCutsSpectra[i][c]->OverridePOT(1);
        DataCutHistos.push_back(dataCutsSpectra[i][c]->ToTH1(1) );
      }
    }
    //mc cuts
    for(int c=0; c<ncuts; c++){
        string cutname=cutnames[c];
        //auto dataCutSpectrum=dataCutsSpectra[i][c];
        //TH1D* chanHisto = channelSpectra[i][jChannel]->ToTH1(TargetPOT);
        //mcCutsSpectra[i][c]->OverridePOT(1);
        mcCutHistos.push_back(mcCutsSpectra[i][c]->ToTH1(TargetPOT) );
      }


	// Manage under/overflow bins
    //What is this?why do we have to manage the over and underflow? theyre defining overflow as a normal bin and not the overflow I guess?
	RecoHisto->SetBinContent(RecoHisto->GetNbinsX(), RecoHisto->GetBinContent(RecoHisto->GetNbinsX()) + RecoHisto->GetBinContent(RecoHisto->GetNbinsX() + 1));
	RecoTrueHisto->SetBinContent(RecoTrueHisto->GetNbinsX(), RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX()) + RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX() + 1));
	RecoBkgHisto->SetBinContent(RecoBkgHisto->GetNbinsX(), RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX()) + RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX() + 1));
	

	RecoHisto->SetBinContent(1, RecoHisto->GetBinContent(0) + RecoHisto->GetBinContent(1));
	RecoTrueHisto->SetBinContent(1, RecoTrueHisto->GetBinContent(0) + RecoTrueHisto->GetBinContent(1));
	RecoBkgHisto->SetBinContent(1, RecoBkgHisto->GetBinContent(0) + RecoBkgHisto->GetBinContent(1));
	
	PlotCanvas->SetTopMargin(0.13);
	PlotCanvas->SetLeftMargin(0.17);
	PlotCanvas->SetRightMargin(0.05);
	PlotCanvas->SetBottomMargin(0.16);
	
	TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
	leg->SetBorderSize(0);
	leg->SetNColumns(3);
	leg->SetTextSize(TextSize*0.8);
	leg->SetTextFont(FontStyle);
	
	TLegendEntry* legReco = leg->AddEntry(RecoHisto,"Reconstructed","l");
	RecoHisto->SetLineColor(kBlue+2);
	RecoHisto->SetLineWidth(4);
	
    styleHistogram(RecoHisto, VarLabels.at(i)); 
    if(debug) cout<<"in loop l201, i="<<i<<endl;    
    //if(runData) styleHistogram(DataHisto, VarLabels.at(i));
    if(debug) cout<<"in loop l203, i="<<i<<endl;    
	// Style histograms
    /*RecoHisto->GetXaxis()->SetTitleFont(FontStyle);                      
      RecoHisto->GetXaxis()->SetLabelFont(FontStyle);                      
      RecoHisto->GetXaxis()->SetNdivisions(8);                             
      RecoHisto->GetXaxis()->SetLabelSize(TextSize);                       
      RecoHisto->GetXaxis()->SetTitleSize(TextSize);                       
      RecoHisto->GetXaxis()->SetTitleOffset(1.1);                          
      RecoHisto->GetXaxis()->CenterTitle();                                
      RecoHisto->GetXaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());
	
      RecoHisto->GetYaxis()->SetTitleFont(FontStyle);
      RecoHisto->GetYaxis()->SetLabelFont(FontStyle);
      RecoHisto->GetYaxis()->SetNdivisions(6);
      RecoHisto->GetYaxis()->SetLabelSize(TextSize);
      RecoHisto->GetYaxis()->SetTitleSize(TextSize);
      RecoHisto->GetYaxis()->SetTitleOffset(1.3);
      RecoHisto->GetYaxis()->SetTickSize(0);
      RecoHisto->GetYaxis()->CenterTitle();
    
    */
	double imax = RecoHisto->GetMaximum();
	double YAxisRange = 1.3*imax;
    if(debug) cout<<"in loop l226, i="<<i<<endl;    
	RecoHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
	RecoTrueHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
	RecoBkgHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
    if(debug) cout<<"in loop l230, i="<<i<<endl;    

	TLegendEntry* legRecoTrue = leg->AddEntry(RecoTrueHisto,"True","l");
	RecoTrueHisto->SetLineColor(kRed+1);
	RecoTrueHisto->SetLineWidth(4);
    if(debug) cout<<"in loop l235, i="<<i<<endl;    

	TLegendEntry* legRecoBkg = leg->AddEntry(RecoBkgHisto,"Background","l");
	RecoBkgHisto->SetLineColor(kOrange+7);
	RecoBkgHisto->SetLineWidth(4);
    if(debug) cout<<"in loop l240, i="<<i<<endl;    

	PlotCanvas->cd();
	RecoHisto->Draw("hist same");
	RecoTrueHisto->Draw("hist same");
	RecoBkgHisto->Draw("hist same");
    if(runData)DataHisto->Draw("same");
	leg->Draw();
    if(debug) cout<<"in loop l248, i="<<i<<endl;    

	// Save as pdf
	PlotCanvas->SaveAs(dir_figs+Form("/Figs/CAFAna/Selection%s.pdf", var.c_str() ) );

	// Save to root file
	if(runData){
      SaveFile->WriteObject(DataHisto, Form("%s_data", var.c_str() ) ); 
      for(int c=0; c<ncuts; c++){
        SaveFile->WriteObject(DataCutHistos[c], Form("%s_%s_data", var.c_str(), cutnames[c].c_str() ));
      }
      
      DataCutHistos.clear();
    }
      for(int c=0; c<ncuts; c++){
        SaveFile->WriteObject(mcCutHistos[c], Form("%s_%s_mc", var.c_str(), cutnames[c].c_str() ));
      }
      mcCutHistos.clear();

    SaveFile->WriteObject(RecoHisto, Form("%s_reco", var.c_str() ) );
	SaveFile->WriteObject(RecoTrueHisto, Form("%s_reco_true", var.c_str() ) );//PlotNames[i]+"_reco_true");
	SaveFile->WriteObject(RecoBkgHisto, Form("%s_bkg", var.c_str() ) );//PlotNames[i]+"_bkg");
    
    if(debug) cout<<"in loop l259, i="<<i<<endl;    

    for(int jChannel=0; jChannel< channelsRecoSignal.size(); jChannel++){
      TH1D* chanHisto = channelSpectra[i][jChannel]->ToTH1(TargetPOT);
      chanHisto->GetXaxis()->SetTitle(("Reconstructed " + VarLabels.at(i)).c_str());
      SaveFile->WriteObject(chanHisto, Form("%s_truthIs_%s", var.c_str(), chanNames[jChannel].c_str() ) );
    }
    if(debug) cout<<"in loop l266, i="<<i<<endl;    

	delete PlotCanvas;
  }
  // Close file
  WritePOT(SaveFile, mcPOT);
  if(runData) cout<<"DataPOT: "<<dataPOT<<endl;
  SaveFile->Close();
    
  // Get histograms for all cuts
}
