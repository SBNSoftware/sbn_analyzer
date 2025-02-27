
std::string tag=std::getenv("XSEC_TAG");

std::string AxisTitle(std::string var)
{
  if (var == "onebin")
    return "single bin";
  else if (var == "Emu")
    return "Muon Energy (GeV)";
  else if (var == "Enu")
    return "Neutrino Energy (GeV)";
  else if (var == "Ep")
    return "Proton Energy (GeV)";
  else if (var == "Epi")
    return "Pion Energy (GeV)";
  else if (var=="EAvail" || var=="Eavail")
    return "Available Energy (Gev)";
  else if (var=="EHad" || var== "Ehad")
    return ("Hadronic Energy");
  else if (var == "Tp")
    return "Proton Kinetic Energy (GeV)";
  else if (var == "Tpi")
    return "Pion Kinetic Energy (GeV)";
  else if (var == "pMu")
    return "Muon Momentum (GeV)";
  else if (var == "pPi")
    return "#pi^{#pm} Momentum (GeV)";
  else if (var == "pMu")
    return "Proton Momentum (GeV)";
  else if (var == "trkLen")
    return "Primary Track Length";
  else if (var == "Ehad")
    return "Hadronic Energy (GeV)";
  else if (var == "VtxEnergy")
    return "Vertex Energy (GeV)";
  else if (var == "Q2")
    return "Q^{2} (GeV/c)^{2}";
  else if (var == "W")
    return "W (GeV/c^{2})";
  else if (var == "y")
    return "Inelasticity";
  else if (var == "x")
    return "Bjorken x";
  else if (var == "xi")
    return "Nachtmann";
  else if (var == "PhiMu")
    return "Muon Phi wrt Beam (deg)";
  else if (var == "CosThetaMu")
    return "Cos(Muon Theta wrt Beam)";
  else if (var == "ThetaMu")
    return "Muon Theta wrt Beam (deg)";
  else if (var == "ThetaMuX")
    return "Muon Theta_{X} wrt Beam (deg)";
  else if (var == "ThetaMuY")
    return "Muon Theta_{Y} wrt Beam (deg)";
  else if (var == "ETheta")
    return "E_{\\mu}*(1-cos(\\theta_{\\mu}))";
  else if (var == "nParticles")
    return "# True Particles";
  else if (var == "pdg")
    return "Longest Track's PDG";
  else if (var == "allpdgs")
    return "All Tracks PDG";
  else if (var == "pdgMu")
    return "Reconstructed Muon Track's True PDG";
  else if (var == "pdgP")
    return "Reconstructed Proton Track's True PDG";
  else if (var == "pdgPi")
    return "Reconstructed #pi^{#pm} Track's True PDG";
  else if (var == "chi2Ps")
    return "#chi^{2} for Proton Hypothesis";
  else if (var == "chi2Pis")
    return "#chi^{2} for Pion Hypothesis";
  else if (var == "chi2Mus")
    return "#chi^{2} for Muon Hypothesis";
  //else if (var == "W")
    //return "W (GeV)";
  //else if (var == "Q2")
   // return "Q2 (GeV)";
  else
    return var;
}

std::vector<std::string> GetSISVarNames()
{
  std::vector<std::string> varnames = {"onebin", "pMu", "EAvail","Q2", "W"};
  return varnames;
}


static const std::vector<std::string> GetVarLabels(std::vector<std::string> varnames){
  std::vector<std::string> VarLabels;
  for (std::string varname : varnames){
    VarLabels.push_back( AxisTitle(varname) );

  }
  return VarLabels;
}


static const std::vector<std::string> GetSISVarLabels(){
  std::vector<std::string> varnames = GetSISVarNames();
  std::vector<std::string> VarLabels=GetVarLabels(varnames);
  return VarLabels;
}




std::pair<std::string, std::string> GetEffSpectrumNames(std::string var){
    std::string specnameA="TrueSignal_"+var;
    std::string specnameB="RecoTrueSignal_"+var;
    return{specnameA, specnameB};
}


//make a function to check if a directory exists and make it if not
/*void checkDir(TString histdir){
  if( 0 != system( Form( "test -d %s", histdir.Data() ) ) )
    {
      int madedir = system( Form( "mkdir -m 755 -p %s", histdir.Data() ) );
            
      if( 0 != madedir )
        Error( "checkDir", Form("Could not make directory '%s'", histdir.Data() ) );
    }
    }*/

