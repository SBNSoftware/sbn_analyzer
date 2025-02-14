#include "Definitions.h"
#include "Utils_noSBNcode.h"
using namespace ana;




const Binning VarBinning(string var, bool isTrue = true)
{
  if (var=="oneBin" || var=="onebin")
    return bEventCount;
  else if (var=="pMu")
    return Binning::Custom({0.1, 0.283, 0.466, 0.649, 0.832, 1.015, 1.2, 1.8, 2.5, 4., 6.});
    //bMuonMomentumBins;
  else if (var == "Emu" || var == "Enu" )//|| var == "pMu")
    return Binning::Simple(50, 0, 3.5);
  else if (var == "EAvail" || var == "Ehad" )//|| var == "pMu")
    return Binning::Simple(15, 0, 4);
  else if (var == "ECompleteness" )//|| var == "pMu")
    return Binning::Simple(25, 0, 1);
  else if (var == "W")
    return Binning::Simple(50, 0, 3);
  else if (var == "Q2")
    return Binning::Simple(50, 0, 2);
  else if (var == "trkLen")
    return Binning::Simple(40, 0, 200);
  else if (var == "nParticles")
    return Binning::Simple(100, 0, 100);
  else if (var == "chi2Ps" || var == "chi2Mus" || var == "chi2Pis" || var == "chi2Ks")
    return Binning::Simple(80, 0, 400);
  else if (var == "ndaughts")
    return Binning::Simple(20, 0, 20);
  else if (var == "Ep" || var == "Epi" || var == "Tpi" || var == "pP" || var == "pPi")
    return Binning::Simple(50, 0, 3.5);
  else if (var == "Tp")
    return Binning::Simple(90, 0, 3.);
  else if (var == "Npi" || var == "Np")
    return Binning::Simple(15, 0, 15);
  else if (var == "pdg" || var == "pdgMu" || var == "pdgPi" || var == "pdgP" || var == "allpdgs" || "candPdgs")
    return Binning::Simple(10, -1, 9, {"None", "e^{-}", "e^{+}", "#mu^{-}", "#mu^{+}", "#pi^{+}", "#pi^{-}", "p", "#gamma", "Other"});
  return Binning::Simple(100, 0, 10);
}


std::tuple<Var,Var, TruthVar> GetVarTuple(std::string var){
  if (var=="onebin"){
    return {kEventCount, kEventCount, kTrueEventCount};
  }
  else if (var=="pMu"){
    return {kMuonMomentum, kRecoTruthMuonMomentum, kTruthMuonMomentum}; 
  }
  else if (var=="W"){
    return {kRecoTruthInvariantMass, kRecoTruthInvariantMass, kTruthInvariantMass};//tmp fix
    //return {kInvariantMass, kRecoTruthInvariantMass, kTruthInvariantMass};  
  }   
  else if (var=="Q2"){
    return {kRecoTruthQ2, kRecoTruthQ2, kTruthQ2};//tmp fix
    //return {kQ2, kRecoTruthQ2, kTruthQ2};  
  }   
  else if (var=="EAvail"){
    return {kEAvail, kRecoTruthEAvail, kTruthEAvail};//tmp fix
    //return {kQ2, kRecoTruthQ2, kTruthQ2};  
  }   
  else if (var=="ECompleteness"){
    return {kEnergyCompleteness, kRecoTruthEnergyCompleteness, kTruthEnergyCompleteness};//tmp fix
    //return {kQ2, kRecoTruthQ2, kTruthQ2};  
  }   
  else if (var=="cosTheta"){
    return {kMuonCosTheta, kRecoTruthMuonCosTheta, kTruthMuonCosTheta};  
  }   

  else{
    cout<<"Variable string "<<var.c_str()<<" not found. Returning one bin event count variable"<<endl;
    return {kEventCount, kEventCount, kTrueEventCount};
  }

}

static const std::vector<std::tuple<Var, Var, TruthVar>> GetVars(std::vector<std::string> varnames)
{
  std::vector<std::tuple<Var, Var, TruthVar>> Vars;
  for (string varname : varnames){
      std::tuple VarTuple = GetVarTuple(varname);
      Vars.push_back(VarTuple);
    }
  return Vars;
}



static const std::vector<std::tuple<Var, Var, TruthVar>> GetSISVars()
{
  std::vector<std::string> varnames = GetSISVarNames();
  std::vector<std::tuple<Var, Var, TruthVar>> Vars = GetVars(varnames);
  //
  //for (string varname : varnames)
  //  {
   //   std::tuple VarTuple = GetVarTuple(varname);
   //   Vars.push_back(VarTuple);
    //}
  return Vars;
}


static const std::vector<Binning> GetSISBins(){
  std::vector<std::string> varnames = GetSISVarNames();
  std::vector<Binning> VarBins;
  for (string varname : varnames){
    VarBins.push_back( VarBinning(varname) );

  }
  return VarBins;
}

void WritePOT(TFile *file, double mcPOT){
  TVector2* pot=new TVector2(TargetPOT, mcPOT);
  file->WriteTObject(pot,"pot");
}


/*SpillVar GetSpillVar(string var, bool isTrue = true, bool isTrkTrue = false)
  {
  if (isTrkTrue)
  return GetSpillVarTruthTree(var);
  if (isTrue)
  {

  if (var == "pMu")
  return kMuTrkTrueP; // This is whats the true energy of the track that is IDed as the best muon track. Not necessarily the highest energy muon.
  else if (var == "pPi")
  return kPiTrkTrueP;
  else if (var == "pP")
  return kPTrkTrueP;

  else if (var == "pdg")
  return kLongestTrkPDG;
  ///==== Fix this ======
  else if (var == "Emu")
  return kMuE; // This isnt really the variable I want, this is the true energy of the true highest energy muon, I want the true energy of the thing that I think is a muon

  else if (var == "bestSlice")
  return kBestNuID;
  // else if(var=="bestSlice") return kTrueBestSlice;

  else if (var == "NPi")
  return kNPi;
  else if (var == "NPiPlus")
  return kNPiPlus;
  else if (var == "NPiMinus")
  return kNPiMinus;
  else if (var == "NPi0")
  return kNPiZero;
  else if (var == "NPiChargedPi0")
  return kNPi;
  else if (var == "Np")
  return kNProton;
  else if (var == "Np50MeV")
  return kN50MeVProtons;

  else if (var == "event")
  return kEvent;
  else if (var == "run")
  return kRun;
  else if (var == "subrun")
  return kSubrun;
  else if (var == "nuIndex")
  return kNuIndex;
  else if (var == "pdg")
  return kLongestTrkPDG;
  else if (var == "pdgMu")
  return kMuTrkPDG;
  else if (var == "pdgP")
  return kPTrkPDG;
  else if (var == "pdgPi")
  return kPiTrkPDG;

  else if (var == "W")
  return kTrueW; // idk that this is really what I want, I think this really is "truth tree" or whatever this distinction I was making
  else if (var == "Q2")
  return kTrueQ2;
  }
  else
  { // reco variables
  // if(var=="Enu") return kNuTrueEnergy;
  if (var == "pdg")
  return kLongestTrkPDGReco;
  else if (var == "pMu")
  return kMuRecoP;
  else if (var == "pPi")
  return kPiRecoP;
  else if (var == "pP")
  return kPRecoP;

  else if (var == "pdgMu")
  return kMuTrkPDG;
  else if (var == "pdgP")
  return kMuTrkPDG;
  else if (var == "pdgPi")
  return kPiTrkPDG;
  else if(var=="Epi") return kChargedPiE;
  else if(var=="Ep") return kProtonE;
  else if(var=="Ep") return kProtonE;
  else if(var=="Tp") return kProtonKinE;


  else if(var=="bestSlice") return kBestNuID;
  //else if(var=="bestSlice") return kTrueBestSlice;

  else if(var=="NPi") return kNPi;
  else if(var=="NPiPlus") return kNPiPlus;
  else if(var=="NPiMinus") return kNPiMinus;
  else if(var=="NPi0") return kNPiZero;
  else if (var=="NPiChargedPi0") return kNPi;
  else if(var=="Np") return kNProton;
  else if(var=="Np50MeV") return kN50MeVProtons;


  else if(var=="event") return kEvent;
  else if(var=="run") return kRun;
  else if(var=="subrun") return kSubrun;
  else if(var=="nuIndex") return kNuIndex;
        
  }
  // else{ //reco variables
  //  }
  // else return NULL;
  // else{
  cout << "Utils::GetSpillVar could not find variable " << var << "." << endl;
  assert(false);
  // throw std::invalid_argument("Utils::GetSpillVarTruthTree does not recognize variable");
  return kRun;
  //}

  // else{

  //}
  }
  */

