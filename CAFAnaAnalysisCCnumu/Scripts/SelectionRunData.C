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
#include <fstream>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionRunData() {
    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    //SpectrumLoader NuLoader(InputFiles);
  //SpectrumLoader NuLoader("/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025Av3/v10_04_06_01/MCP2025Av3_DevSample/flatcaf/bnblight/*/*.flat.caf.root");

  std::ifstream infiles("list_Data_MCP2025v3_DevSample.txt");
  string infile;
  vector<string> fnames; 
  while (std::getline(infiles, infile)){
    fnames.push_back(infile);
    //tree->Add("/pnfs/sbn/data_add/sbn_nd/poms_production/data/MCP2025Av3/v10_04_06_01/MCP2025Av3_DevSample/flatcaf/bnblight/*/*.flat.caf.root");                                                                                          
    //tree->Add((infiledir+"/*.flat.caf.root").c_str() );
  }

  SpectrumLoader NuLoader(fnames);
    // Open csv file to store data
    TString FilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/EventData.csv";
    fstream file; file.open(FilePath, fstream::out); 
    file << "fno,run,subrun,evt,subevt" << std::endl;
    file.close();

    Spectrum sRecoSignals(
        "RecoSignals",
        bEventCount,      //binning
        NuLoader,
        kEventCount,
        kNoSpillCut,
        kNoCut
        //kFirstCut
    );

    NuLoader.Go();

    sRecoSignals.OverridePOT(1);
    cout<<"Livetime: "<< sRecoSignals.Livetime()<<endl;
    TH1D* honebin=sRecoSignals.ToTH1(1);
    honebin->SaveAs("singlebinData.root");
}
