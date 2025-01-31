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

  string UserName=std::getenv("USER");
  TString RootFilePath = "/exp/sbnd/data/users/" + UserName + "/CAFAnaOutput/SelectionEfficiency.root";
  TFile* fin=new TFile(RootFilePath);
  std::vector<std::pair<std::string, std::string>> vars2d;//
  vars2d.push_back({"Q2","W"}); 
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
  map<string, TH2D*> h2ds;
  for(pair<string, string> var: vars2d){
    //for(string cut : cuts){
    string nameNum=  "EffNum_"+var.second+"_vs_"+var.first;
    //cout<<name<<endl;
    names.push_back(nameNum);
    TH2D* hnum=(TH2D*) fin->Get(nameNum.c_str() );
    assert(hnum || Form("Didn't find TH2D %s", nameNum.c_str() ) );
    //h2ds.push_back(h);
    h2ds[nameNum]=hnum;

    string nameDen= "EffDenom_"+var.second+"_vs_"+var.first;
    //cout<<name<<endl;
    names.push_back(nameDen);
    TH2D* hden=(TH2D*) fin->Get(nameDen.c_str() );
    if(hden==NULL) cout<<"denom histogram does not exist"<<endl;
    assert(hden || Form("Didn't find TH2D %s", nameDen.c_str() ) );
    //h2ds.push_back(h);
    h2ds[nameDen]=hden;

    //}
  }


  //hists.push_back(npinp);
  MnvPlotter plotter(kNukeCCStyle);
  plotter.print_topdir="/exp/sbnd/data/users/afilkins/SISPlots";

 
  TCanvas* c= new TCanvas("c", "c",1280, 800);
  gStyle->SetPalette(kCool);
  c->SetRightMargin(0.15); //for colz plots
  for( string name: names){
    cout<<name<<endl;
    if(h2ds[name]==NULL) cout<<"histogram does not exist"<<endl;
    //TH2D* h= h2ds[name];//(TH2D*) h2ds[name]->Clone();
    //assert(h || Form("Didn't find TH2D %s", name.c_str() ) ); 
    //cout<<name<<endl;
    //if(name.find("chi2Ps") != string::npos)       h->GetXaxis()->SetTitle("#chi^{2} for Proton Hypothesis");
    // if(name.find("chi2Mus") != string::npos)       h->GetXaxis()->SetTitle("#chi^{2} for Muon Hypothesis");
    //if(name.find("chi2Pis") != string::npos)       h->GetXaxis()->SetTitle("#chi^{2} for Pion Hypothesis");
    
    /*if(name.find("v_allpdgs")!=string::npos){
      vector<string> labels= GetPDGLabels();
      for(int b=1; b<h->GetXaxis()->GetNbins(); b++){
	h->GetXaxis()->SetBinLabel(b, labels[b-1].c_str());
	h->GetYaxis()->SetBinLabel(b, labels[b-1].c_str());
      }
      }*/

    h2ds[name]->Draw("colz");
    plotter.MultiPrint(c, name);
    plotter.DrawNormalizedMigrationHistogram(h2ds[name]);
    plotter.MultiPrint(c, "rownorm_"+name);
    plotter.DrawColumnNormalizedMigrationHistogram(h2ds[name]);
    plotter.MultiPrint(c, "columnnorm_"+name);


  }

  //For Q2 and W specifically
  TH2D* hWQ2=(TH2D*) h2ds["EffNum_W_vs_Q2"]->Clone();
  double recoTot=h2ds["EffNum_W_vs_Q2"]->Integral();
  double recoSIS=h2ds["EffNum_W_vs_Q2"]->Integral( hWQ2->GetXaxis()->FindBin(sis_minQ2), hWQ2->GetXaxis()->GetNbins(), hWQ2->GetYaxis()->FindBin(sis_minW), hWQ2->GetYaxis()->FindBin(sis_maxW) );

  double trueTot=h2ds["EffDenom_W_vs_Q2"]->Integral();
  double trueSIS=h2ds["EffDenom_W_vs_Q2"]->Integral( hWQ2->GetXaxis()->FindBin(sis_minQ2), hWQ2->GetXaxis()->GetNbins(), hWQ2->GetYaxis()->FindBin(sis_minW), hWQ2->GetYaxis()->FindBin(sis_maxW) );
  
  cout<<"Numerator CC: "<< recoTot<< "\t Numerator SIS: "<< recoSIS<<endl;
  cout<<"Denominator CC: "<< trueTot<< "\t Denominator SIS: "<< trueSIS<<endl;
  cout<<"True MC POT: "<< pot->Y()<<" \t Scaled to MC POT: "<<pot->X()<<endl;
  
  
  //I should fix this function to actually work correctly
  //int pbin=GetPDGBinNum(P_PDG)+2;
  //int pibin=GetPDGBinNum(PI_PDG)+2;
  //cout<<"pibin"<<pibin<<endl;
  //int mubin=GetPDGBinNum(MU_PDG)+2;
  
  /* c->SetRightMargin(0.1); //for colz plots
  for(int i=0; i<names.size(); i++){//(string name: names){
    string name=names[i];
    pair<string,string> varxy=vars2d[i/2];
    TH1D* h_x= h2ds[name]->ProjectionX(Form("_%s", varxy.first ));
    TH1D* h_Y= h2ds[name]->ProjectionY(Form("_%s", varxy.first ));

    if( name=="reco_ccpip_allpdgs_v_chi2Ps"){
      cout<<pbin<<endl;
      cout<<"2D integral"<<h2ds[name]->Integral()<<endl;
      h_p->Print("all");
      cout<<h_p->GetMaximum()<<endl;
    }
    if(name.find("chi2") != string::npos){ //it says chi2 in there somewhere
      h_mu->GetXaxis()->SetRangeUser(0, 100);
    }
    h_p->SetLineColor(kBlue+3);
    h_pi->SetLineColor(kMagenta+1);
    h_mu->SetLineColor(kRed+2);

    TLegend *leg=new TLegend( 0.3, 0.8, 0.6, 0.9);
    leg->AddEntry(h_mu, "True Muons", "l");
    leg->AddEntry(h_pi, "True Charged Pions", "l");
    leg->AddEntry(h_p, "True Protons", "l");
    
    //if(name.find("chi2Ps") != string::npos){
    double max=h_p->GetMaximum();
    if(h_mu->GetMaximum() >max) max=h_mu->GetMaximum();
    h_p->SetAxisRange(0, 1.1*max);
      //cout<<name<<" max: "<<max<<endl;
      //h_p->GetYaxis()->SetRangeUser(0, max/4.);
      //}
    h_p->Draw("hist");
    h_mu->Draw("histsame");
    h_pi->Draw("histsame");
    //h_p->Draw("histsame");
    leg->Draw("same");
    plotter.MultiPrint(c, "byParticle1D_"+name);
    }*/


  return 0;
}
