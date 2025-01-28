#ifndef CaloEhadTools_h
#define CaloEhadTools_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

#include "TMath.h"
#include <TVector3.h>
#include <TLorentzVector.h>

class CaloEhadTools {
 private: 
  double fMuonMomentum;
  double fMuonCosTheta;
  double fHadronMomentum;
  double fHadronCosTheta;
  //These really need to go but tmp hacks.
  double fLeadingProtonMomentum;
  double fLeadingProtonCosTheta;
  double fRecoilProtonMomentum;
  double fRecoilProtonCosTheta;
  double fCosOpeningAngleProtons;
  double fCosOpeningAngleMuonTotalProton;
  double fTransverseMomentum;
  double fDeltaAlphaT;
  double fMissingMomentum;
  double fAlphaThreeD;
  double fCosOpeningAngleMomentumTransferTotalProton;


  double fInvariantMass;
  double fCosOpeningAngleLProtonMuon;
  double fCosOpeningAngleRProtonMuon;

 public:
  // Default constructor
  CaloEhadTools(TVector3 MuonVector, vector<TVector3> HadronVectors, vector<int> HadronPIDs);
  // Default destructor
  ~CaloEhadTools() = default;

  // Getter functions
  double ReturnMuonMomentum();
  double ReturnMuonCosTheta();
  double ReturnLeadingProtonMomentum();
  double ReturnLeadingProtonCosTheta();
  double ReturnRecoilProtonMomentum();
  double ReturnRecoilProtonCosTheta();
  double ReturnCosOpeningAngleProtons();
  double ReturnCosOpeningAngleMuonTotalProton();
  double ReturnTransverseMomentum();
  double ReturnDeltaAlphaT();
  double ReturnMissingMomentum();
  double ReturnAlphaThreeD();
  double ReturnCosOpeningAngleMomentumTransferTotalProton();
  double ReturnInvariantMass();
  double ReturnCosOpeningAngleLProtonMuon();
  double ReturnCosOpeningAngleRProtonMuon();
};

#endif
