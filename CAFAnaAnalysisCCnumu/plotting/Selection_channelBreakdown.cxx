#include <TROOT.h>
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <string>
#include <assert.h>
#include "TStyle.h"
#include "localColor.h"
#include "TLegend.h"
#include "TLine.h"
#include "TVector2.h"
#include <PlotUtils/MnvPlotter.h>
#include "../Scripts/Utils_noSBNcode.h"
//#include "../Scripts/Definitions.h"

using namespace std;
using namespace PlotUtils;

int main(int args, char* argv[]){
  double gev_to_mev=1.e3;
  double mev_to_gev=1.e-3;

  double sis_minW=1.5;
  double sis_maxW=2.;
  double sis_minQ2=1.;
  
  bool useTrueSIS=true;
  if(args> 1){
    useTrueSIS=atoi(argv[1]); 
  }

  string UserName=std::getenv("USER"); 
  string tag=std::getenv("XSEC_TAG");
  TString RootFilePath = "/exp/sbnd/data/users/" + UserName + "/CAFAnaOutput/"+tag+"/Selection.root";//tag.c_str()+"/SelectionEfficiency.root";
  if(useTrueSIS) RootFilePath = "/exp/sbnd/data/users/" + UserName + "/CAFAnaOutput/"+tag+"/SelectionSIS.root";
  TFile* fin=new TFile(RootFilePath);
  std::vector<std::string> vars=GetSISVarNames();
  std::vector<std::string> channels=getChannelShortNames();
  std::vector<std::string> channel_longnames=getChannelLongNames();
  //std::vector<std::pair<std::string, std::string>> vars2d;//
  //vars2d.push_back({"Q2","W"}); 
  std::vector<std::pair<std::string,std::string>> vars2d={{"pMu", "onebin"}, {"EAvail", "ECompleteness"}, {"Q2", "W"}, {"EAvail", "Enu"}, {"EAvail","Ehad"}, {"pMu", "Enu"}, {"sumE", "Enu"}, {"Ehad", "EhadSummed"}};
  if(useTrueSIS==false) vars2d={{"pMu", "onebin"}};

  string plotTitle= useTrueSIS? "Reco CC #nu_{#mu} SIS" : "Reco  CC #nu_{#mu}" ;
  //string plotTitleDen= useTrueSIS? "True SIS" : "True CC #nu_{#mu}" ;
  Int_t colors[]={30,11,46,6,8,4};
  //{4,8,6,38,11,30};//{154, 153,155,156,195, 190};
  //vars2d.push_back({"EAvail", "ECompleteness"});
  //{{"candPdgs", "ndaughts"}, {"allpdgs", "chi2Ps"}, {"allpdgs", "chi2Mus"}, {"allpdgs", "chi2Pis"}, {"allpdgs", "allpdgs_reco"} };
  //want to plot ccpip_allpdgs_v_chi2Ps in projections 
  /*vector<string> cuts={"kCCNuMu", "kCCPiCharged", "kCCPiProton", "kCCProton"};
  int ncuts=cuts.size();
  bool addrecoCuts=true;
  if (addrecoCuts){
    for(int c=0; c<ncuts; c++){
      cuts.push_back(cuts[c]+"Reco");
    }
    }*/

  //TH1D* cc=(TH1D*)fin->Get("NuEnergy_NuMu");
  //TH1D* npi=(TH1D*) fin->Get("NuEnergy_NPi");
  //TH1D* npinp=(TH1D*) fin->Get("NuEnergy_NpiNp");
  //TH2D* hNuEPE=(TH2D*)fin->Get("NuEnergy_vs_ProtonKE_NPi");

  TVector2* pot=(TVector2*)fin->Get("pot");


  vector<string> names;
  map<string, vector<TH1D*>> hchannel; //hchannel[varname][iChannel]
 map<string, TObjArray*> hchannelArrays; //hchannel[varname][jChannel]
 map<string, TH1D*> hRecoAll;
 map<string, TH1D*> hreco_org;
 //mapTH1D* sumReco=(TH1D*) fin->Get(
 for(int iVar=0; iVar<vars.size(); iVar++){//(string var: vars){
    string var=vars[iVar];
    cout<<"var:"<<var<<endl;
    hchannel[var]={};//just make an empty vector so I can push back later
    hchannelArrays[var]=new TObjArray();
    
    //TH1D* hrecoog=(TH1D*) fin->Get(Form("%_reco", var.c_str() ));
    //assert(hrecoog || Form("Didn't find TH1D %s_reco", var.c_str() ) );

    //TH1D* hreco=(TH1D*) fin->Get(Form("%_truthIs_other", var.c_str() ));
    //hRecoAll[var]=hrecoog;
    //TH1D* htest=hRecoAll[var];
    //htest->Reset();
    //hRecoAll[var]->Reset();
    //assert(hreco || Form("Didn't find TH1D %s_truthIs_other", var.c_str() ) );

    
    //for(int jChan=0; jChan<channels.size(); jChan++){//(string chan :channels){
    for(int jChan=channels.size()-1; jChan>=0; jChan--){//so SIS is on top of legend
      string chan=channels[jChan];
      string name=  var+"_truthIs_"+chan;
      cout<<name<<endl;
      names.push_back(name);
      TH1D* hchan=(TH1D*) fin->Get(name.c_str() );
      assert(hchan || Form("Didn't find TH1D %s", name.c_str() ) );
      if(var=="W") hchan->GetXaxis()->SetRangeUser(1.2, 2.2);
      //hRecoAll[var]->Add(hchan);
      //hchannel.push_back(h);
      cout<<channel_longnames[jChan]<<endl;
      hchan->SetTitle(channel_longnames[jChan].c_str() ); //plotTitle.c_str());
      hchannel[var].push_back(hchan);
      //add some things together before going to channel breakdown 
      if(chan=="lowestW") continue;//gonna add to transition
      if(chan=="trans") {     
        hchan->Add(hchannel[var][kLowestW]);
        hchan->SetTitle("W < 1.5 GeV");
        hchannelArrays[var]->Add(hchan);
      }
      else   hchannelArrays[var]->Add(hchan);
      //if(var=="Q2") cout<<"Channel "<<channel_longnames[jChan]<<": "<<hchan->Integral()/hrecoog->Integral()<<endl;

    }

        /*enum kChannels{
    kSIS=0,
    kDIS,
    kTrans,
    kLowestW,
    kLowQ2,
    kNotCC,
    kOther
    };*/

      //}
    
  }


  //hists.push_back(npinp);
  MnvPlotter plotter(kNukeCCStyle);
  gStyle->SetOptTitle(0);//I want titles even though plotutils defaults to them off
  plotter.print_topdir="/exp/sbnd/data/users/afilkins/SISPlots/"+tag;
  if(useTrueSIS) plotter.print_topdir="/exp/sbnd/data/users/afilkins/SISPlots/TrueSIS/";//+tag.c_str();
 
  plotter.legend_fill_color=-1; // <0 makes it transparent

  TCanvas* c= new TCanvas("c", "c",1280, 800);
  gStyle->SetPalette(kCool);
  c->SetRightMargin(0.15); //for colz plots

  for( string var: vars){
    cout<<var<<endl;
    if(hchannel[var][kSIS]==NULL) cout<<"histogram does not exist"<<endl;

    /*void MnvPlotter::DrawStackedMC(const TObjArray* mcHists,
                                   const Double_t mcScale,
                                   const std::string& legPos,
                                   const Int_t mcBaseColor,               // Color of first stacked histo.                                                                                                                                   
                                   const Int_t mcColorOffset,             // 2nd histo in stack has color mcBaseColor+mcColorOffset, etc.                                                                                                    
                                   const Int_t mcFillStyle,               // Fill style for histograms only.                                                                                                                                 
                                   const char* xaxislabel,
                                   const char* yaxislabel)
    */
    string plotname="selection_chanStack_"+var;
    plotter.DrawStackedMC(hchannelArrays[var],colors,  1. , "TC");//, 152, 1);


    plotter.MultiPrint(c, plotname);



    //hRecoAll[var]->Draw("e2same");
    //plotter.DrawDataStackedMC(hRecoAll[var], hchannelArrays[var], 1. , "TL");//, 152, 1);
    plotter.MultiPrint(c, plotname+"_uncert");

    //DrawColumnNormalizedMigrationHistorgram(h_migration, bool drawAsMatrix, bool coarseContors=false, bool includeUnder/OverFlows=true, bool noText=false)
  

  }
  


  return 0;
}
