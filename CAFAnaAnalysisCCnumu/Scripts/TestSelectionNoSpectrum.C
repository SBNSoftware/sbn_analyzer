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
#include "TChain.h"

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


//int main(){
void TestSelectionNoSpectrum(){
  TChain* tree=new TChain("recTree");
  std::ifstream infileDirs("listDirs_Data_MCP2025v3_DevSample.txt");
  string infiledir;
  while (std::getline(infileDirs, infiledir)){
      //tree->Add("/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025Av3/v10_04_06_01/MCP2025Av3_DevSample/flatcaf/bnblight/*/*.flat.caf.root");
    tree->Add((infiledir+"/*.flat.caf.root").c_str() );
  }
  Float_t startx;
  //double* p_startx= &startx;
  int nslc, npfp;
  tree->SetBranchAddress("rec.slc.reco.pfp.trk.start.x", &startx);
  tree->SetBranchAddress("rec.nslc", &nslc);
  //tree->SetBranchAddress("rec.slc.npfp", &npfp);

  int nentries=tree->GetEntries();
  cout<<"nEntries: "<<nentries<<endl;
  TH1D* hStartx=new TH1D("hStartx", ";start.x (cm); Entries", 200, -1000,1000);
  for( int iEntry=0; iEntry<nentries; iEntry++){
    tree->GetEntry(iEntry);
    for(int iSlc=0; iSlc<nslc; iSlc++){
   
      //for(int iPfp=0; iPfp<npfp; iPfp++){
        hStartx->Fill(startx);

      
      }
    
    //}
  


  }

  hStartx->SaveAs("hStartX.root");
}
