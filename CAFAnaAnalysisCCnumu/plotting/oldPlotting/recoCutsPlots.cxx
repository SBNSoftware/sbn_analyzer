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
#include <PlotUtils/MnvPlotter.h>
#include "CommonIncludes.h"

using namespace std;
using namespace PlotUtils;

int main(int args, char* argv[]){
  double gev_to_mev=1.e3;
  double mev_to_gev=1.e-3;
  TFile* fin=new TFile("Spectra_selectionOpt.root");
  vector<pair<string, string>> vars2d= {{"candPdgs", "ndaughts"}, {"allpdgs", "chi2Ps"}, {"allpdgs", "chi2Mus"}, {"allpdgs", "chi2Pis"}, {"allpdgs", "allpdgs_reco"} };
  //want to plot ccpip_allpdgs_v_chi2Ps in projections 
  vector<string> cuts={"kCCNuMu", "kCCPiCharged", "kCCPiProton", "kCCProton"};
  int ncuts=cuts.size();
  bool addrecoCuts=true;
  if (addrecoCuts){
    for(int c=0; c<ncuts; c++){
      cuts.push_back(cuts[c]+"Reco");
    }
  }

  //TH1D* cc=(TH1D*)fin->Get("NuEnergy_NuMu");
  //TH1D* npi=(TH1D*) fin->Get("NuEnergy_NPi");
  //TH1D* npinp=(TH1D*) fin->Get("NuEnergy_NpiNp");
  //TH2D* hNuEPE=(TH2D*)fin->Get("NuEnergy_vs_ProtonKE_NPi");

  vector<string> names;
  map<string, TH2D*> h2ds;
  for(pair<string, string> var: vars2d){
    for(string cut : cuts){
      string name= GetSpectrumName( var.first+"_v_"+var.second, cut);
      //cout<<name<<endl;
      names.push_back(name);
      TH2D* h=(TH2D*) fin->Get(name.c_str() );
      assert(h || Form("Didn't find TH2D %s", name.c_str() ) );
      //h2ds.push_back(h);
      h2ds[name]=h;
    }
  }


  //hists.push_back(npinp);
  MnvPlotter plotter(kNukeCCStyle);
  plotter.print_topdir="/exp/sbnd/data/users/afilkins/XSec/RecoPlots";

 
  TCanvas* c= new TCanvas("c", "c",1280, 800);
  gStyle->SetPalette(kCool);
  c->SetRightMargin(0.15); //for colz plots
  for( string name: names){
    TH2D* h= h2ds[name];
    assert(h || Form("Didn't find TH2D %s", name.c_str() ) ); 
    //cout<<name<<endl;
    if(name.find("chi2Ps") != string::npos)       h->GetXaxis()->SetTitle("#chi^{2} for Proton Hypothesis");
    if(name.find("chi2Mus") != string::npos)       h->GetXaxis()->SetTitle("#chi^{2} for Muon Hypothesis");
    if(name.find("chi2Pis") != string::npos)       h->GetXaxis()->SetTitle("#chi^{2} for Pion Hypothesis");
    
    if(name.find("v_allpdgs")!=string::npos){
      vector<string> labels= GetPDGLabels();
      for(int b=1; b<h->GetXaxis()->GetNbins(); b++){
	h->GetXaxis()->SetBinLabel(b, labels[b-1].c_str());
	h->GetYaxis()->SetBinLabel(b, labels[b-1].c_str());
      }
    }

    h->Draw("colz");
    plotter.MultiPrint(c, name);
    plotter.DrawNormalizedMigrationHistogram(h);
    plotter.MultiPrint(c, "rownorm_"+name);
    plotter.DrawColumnNormalizedMigrationHistogram(h);
    plotter.MultiPrint(c, "columnnorm_"+name);
  }
  
  //I should fix this function to actually work correctly
  int pbin=GetPDGBinNum(P_PDG)+2;
  int pibin=GetPDGBinNum(PI_PDG)+2;
  cout<<"pibin"<<pibin<<endl;
  int mubin=GetPDGBinNum(MU_PDG)+2;
  
  c->SetRightMargin(0.1); //for colz plots
  for(string name: names){//anti particles are one bin larger
      
    TH1D* h_p= h2ds[name]->ProjectionX("_p", pbin, pbin);
    TH1D* h_pi= h2ds[name]->ProjectionX("_pi", pibin, pibin+1);
    TH1D* h_mu= h2ds[name]->ProjectionX("_mu", mubin, mubin+1);

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
  }


  return 0;
}
