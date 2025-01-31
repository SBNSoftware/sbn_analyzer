#include <TROOT.h>
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <string>
#include "TStyle.h"
#include "localColor.h"
#include "TLegend.h"
#include "TLine.h"
#include <PlotUtils/MnvPlotter.h>

using namespace std;
using namespace PlotUtils;

int main(int args, char* argv[]){
  double gev_to_mev=1.e3;
  double mev_to_gev=1.e-3;
  TFile* fin=new TFile("Spectra.root");
  TH1D* cc=(TH1D*)fin->Get("NuEnergy_NuMu");
  TH1D* npi=(TH1D*) fin->Get("NuEnergy_NPi");
  TH1D* npinp=(TH1D*) fin->Get("NuEnergy_NpiNp");
  TH2D* hNuEPE=(TH2D*)fin->Get("NuEnergy_vs_ProtonKE_NPi");

  vector<TH1D*> hists;
  hists.push_back(cc);
  hists.push_back(npi);
  //hists.push_back(npinp);
  MnvPlotter plotter(kNukeCCStyle);
  plotter.print_topdir="/exp/sbnd/data/users/afilkins/XSec/";

  int nbinsy=hNuEPE->GetYaxis()->GetNbins();
  for(int b=0; b<20; b++){
    double cut=b*10 *mev_to_gev; //10 MeV increments
    int bin=hNuEPE->GetYaxis()->FindBin(cut); 
    cut=hNuEPE->GetYaxis()->GetBinLowEdge(bin); 
    cout<<"bin: "<<bin<<endl;
    TH1D* hpcut=hNuEPE->ProjectionX(Form("NuE_p%f MeV", cut*gev_to_mev), bin, nbinsy+1);
    hpcut->SetTitle(Form("T_{p}>%0.1f MeV",cut*gev_to_mev) ); 
    hists.push_back(hpcut);
  }
  hists[3]->Print();
  hists[5]->Print();
  //hists[0]->SetTitle("CC #nu_#mu \t \t");
  //hists[1]->SetTitle("CC #nu_#mu #geq1#pi^{#pm} \t");
  TCanvas* c= new TCanvas("c", "c",1280, 800);
  cc->Draw();
  vector<int> colors=getColors(2);
  TLegend *leg=new TLegend(0.6, 0.95, 0.6, 0.95);
  leg->SetNColumns(3);
  for(int h=0; h<hists.size(); h++){
    TH1D* hist=hists[h];
    string title=hist->GetTitle();
    
    hist->SetLineColor(colors[h]);
    hist->SetFillColor(colors[h]);
    if(h==0){
      hist->SetLineColor(kBlack);
      hist->SetFillColor(kGray);
    }

    if(h==1){
      hist->SetLineColor(kBlue+2);
      hist->SetFillColor(kBlue-4);
    }
    //if(h==2){
    // hist->SetLineColor(kTeal-5);
    //hist->SetFillColor(kTeal-4);
    //}

    int tot=hist->Integral();
    leg->AddEntry(hist, Form("%s \t %d",title.c_str(), tot ), "lf" );
    hist->Draw("histsame");
  }

  hists[0]->SetTitle("True Events 2.5e17 POT");
  //c->Update();



  //hists[0]->SetFillStyle(0);
  hists[0]->Draw("ap same");
  leg->Draw();
  c->SaveAs("CCTrueEventStack.png");

  TCanvas* c2= new TCanvas("c", "c",1280, 800);
  npi->Draw();
  TLegend *leg2=new TLegend(0.6, 0.95, 0.6, 0.95);
  leg2->SetNColumns(3);

  for(int h=1; h<hists.size(); h++){
    TH1D* hist=hists[h];
    string title=hist->GetTitle();
    
    int tot=hist->Integral();
    leg2->AddEntry(hist, Form("%s \t %d",title.c_str(), tot ), "lf" );
    hist->Draw("histsame");
  }

  hists[1]->SetTitle("True Events 2.5e17 POT");
  //c->Update();



  //hists[0]->SetFillStyle(0);
  hists[1]->Draw("ap same");
  leg2->Draw();
  c2->SaveAs("CCTruePiEventStack.png");

  gStyle->SetPalette(kCool);
  c2->SetRightMargin(0.15);
  hNuEPE->GetYaxis()->SetTitleOffset(0.6);
  hNuEPE->GetZaxis()->SetTitleOffset(0.6);
  hNuEPE->Draw("colz");
  double m_p=0.938;
  double m_pi=0.139;
  //  TLine* l=new TLine(m_p+m_pi, 0, 3, 3-m_p-m_pi);
  //l->Draw("same");
  //TLegend* l3= new TLegend(0.2, 0.8, 0.3, 0.9);
  //l3->AddEntry(l,"T_{p}=E_{#nu}-m_{p}-m{#pi}", "l");
  c2->SaveAs("Proton_v_NeutrinoEnergy_Npi.png");

  return 0;
}
