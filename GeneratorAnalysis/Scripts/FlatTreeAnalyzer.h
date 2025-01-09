#ifndef FlatTreeAnalyzer_h
#define FlatTreeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include "../../Utils/Constants.h"

class FlatTreeAnalyzer {

private:
    TFile* fFile;
    TString fInputFile;
    TString fOutputFile;

public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
    Int_t           Mode;
    Char_t          cc;
    Int_t           PDGnu;
    Float_t         Enu_true;
    Int_t           tgt;
    Int_t           tgta;
    Int_t           tgtz;
    Int_t           PDGLep;
    Float_t         ELep;
    Float_t         CosLep;
    Int_t           nfsp;
    Float_t         px[100];   //[nfsp]
    Float_t         py[100];   //[nfsp]
    Float_t         pz[100];   //[nfsp]
    Float_t         E[100];   //[nfsp]
    Int_t           pdg[100];   //[nfsp]
    Float_t         px_init[10];   //[ninitp]
    Float_t         py_init[10];   //[ninitp]
    Float_t         pz_init[10];   //[ninitp]
    Float_t         E_init[10];   //[ninitp]
    Int_t           pdg_init[10];   //[ninitp]
    Int_t           nvertp;
    Float_t         px_vert[100];   //[nvertp]
    Float_t         py_vert[100];   //[nvertp]
    Float_t         pz_vert[100];   //[nvertp]
    Float_t         E_vert[100];   //[nvertp]
    Int_t           pdg_vert[100];   //[nvertp]
    Float_t         Weight;
    Double_t        fScaleFactor;

    // List of branches
    TBranch        *b_Mode;   //!
    TBranch        *b_cc;   //!
    TBranch        *b_PDGnu;   //!
    TBranch        *b_Enu_true;   //!
    TBranch        *b_tgt;   //!
    TBranch        *b_tgta;   //!
    TBranch        *b_tgtz;   //!
    TBranch        *b_PDGLep;   //!
    TBranch        *b_ELep;   //!
    TBranch        *b_CosLep;   //!
    TBranch        *b_nfsp;   //!
    TBranch        *b_px;   //!
    TBranch        *b_py;   //!
    TBranch        *b_pz;   //!
    TBranch        *b_E;   //!
    TBranch        *b_pdg;   //!
    TBranch        *b_ninitp;   //!
    TBranch        *b_px_init;   //!
    TBranch        *b_py_init;   //!
    TBranch        *b_pz_init;   //!
    TBranch        *b_E_init;   //!
    TBranch        *b_pdg_init;   //!
    TBranch        *b_nvertp;   //!
    TBranch        *b_px_vert;   //!
    TBranch        *b_py_vert;   //!
    TBranch        *b_pz_vert;   //!
    TBranch        *b_E_vert;   //!
    TBranch        *b_pdg_vert;   //!
    TBranch        *b_Weight;   //!
    TBranch        *b_fScaleFactor;   //!

    FlatTreeAnalyzer(TString in, TString out, TTree *tree=0);
    virtual ~FlatTreeAnalyzer();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FlatTreeAnalyzer_cxx
FlatTreeAnalyzer::FlatTreeAnalyzer(TString InputFile, TString OutputFile, TTree *tree) : fChain(0) 
{
    fInputFile = InputFile;
    fOutputFile = OutputFile;
    TString Directory = "/pnfs/sbnd/persistent/users/" + (TString)Constants::UserName + "/HighSamples/";
    TString FullName = Directory + fInputFile;

    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(FullName+".root");
        if (!f || !f->IsOpen()) {
            f = new TFile(FullName+".root");
        }
        f->GetObject("FlatTree_VARS",tree);
        fFile = f;

    }
    Init(tree);
}

FlatTreeAnalyzer::~FlatTreeAnalyzer()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t FlatTreeAnalyzer::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t FlatTreeAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void FlatTreeAnalyzer::Init(TTree *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("Mode", &Mode, &b_Mode);
    fChain->SetBranchAddress("cc", &cc, &b_cc);
    fChain->SetBranchAddress("PDGnu", &PDGnu, &b_PDGnu);
    fChain->SetBranchAddress("Enu_true", &Enu_true, &b_Enu_true);
    fChain->SetBranchAddress("tgt", &tgt, &b_tgt);
    fChain->SetBranchAddress("tgta", &tgta, &b_tgta);
    fChain->SetBranchAddress("tgtz", &tgtz, &b_tgtz);
    fChain->SetBranchAddress("PDGLep", &PDGLep, &b_PDGLep);
    fChain->SetBranchAddress("ELep", &ELep, &b_ELep);
    fChain->SetBranchAddress("CosLep", &CosLep, &b_CosLep);
    fChain->SetBranchAddress("nfsp", &nfsp, &b_nfsp);
    fChain->SetBranchAddress("px", px, &b_px);
    fChain->SetBranchAddress("py", py, &b_py);
    fChain->SetBranchAddress("pz", pz, &b_pz);
    fChain->SetBranchAddress("E", E, &b_E);
    fChain->SetBranchAddress("pdg", pdg, &b_pdg);
    fChain->SetBranchAddress("px_init", px_init, &b_px_init);
    fChain->SetBranchAddress("py_init", py_init, &b_py_init);
    fChain->SetBranchAddress("pz_init", pz_init, &b_pz_init);
    fChain->SetBranchAddress("E_init", E_init, &b_E_init);
    fChain->SetBranchAddress("pdg_init", pdg_init, &b_pdg_init);
    fChain->SetBranchAddress("nvertp", &nvertp, &b_nvertp);
    fChain->SetBranchAddress("px_vert", px_vert, &b_px_vert);
    fChain->SetBranchAddress("py_vert", py_vert, &b_py_vert);
    fChain->SetBranchAddress("pz_vert", pz_vert, &b_pz_vert);
    fChain->SetBranchAddress("E_vert", E_vert, &b_E_vert);
    fChain->SetBranchAddress("pdg_vert", pdg_vert, &b_pdg_vert);
    fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
    fChain->SetBranchAddress("fScaleFactor", &fScaleFactor, &b_fScaleFactor);

    Notify();
}

Bool_t FlatTreeAnalyzer::Notify()
{
    return kTRUE;
}

void FlatTreeAnalyzer::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t FlatTreeAnalyzer::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef FlatTreeAnalyzer_cxx
