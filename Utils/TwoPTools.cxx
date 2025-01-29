
#ifndef TwoPTools_cxx
#define TwoPTools_cxx

#include <TMath.h>
#include <TVector3.h>

#include "TwoPTools.h"

using namespace std;

TwoPTools::TwoPTools(TVector3 MuonVector, TVector3 LeadingProtonVector, TVector3 RecoilProtonVector) {
    // Swap lead and recoil protons to make sure lead has higher magnitude
    if (RecoilProtonVector.Mag() > LeadingProtonVector.Mag()) std::swap(LeadingProtonVector, RecoilProtonVector); 

    // Momenta and theta variables
    fMuonMomentum = MuonVector.Mag();
    fMuonCosTheta = MuonVector.CosTheta();
    fLeadingProtonMomentum = LeadingProtonVector.Mag();
    fLeadingProtonCosTheta = LeadingProtonVector.CosTheta();
    fRecoilProtonMomentum = RecoilProtonVector.Mag();
    fRecoilProtonCosTheta = RecoilProtonVector.CosTheta();

    // Opening angle variables
    TVector3 TotalProtonVector = LeadingProtonVector + RecoilProtonVector;
    fCosOpeningAngleProtons = std::cos(LeadingProtonVector.Angle(RecoilProtonVector));
    fCosOpeningAngleMuonTotalProton = std::cos(MuonVector.Angle(TotalProtonVector));
    fCosOpeningAngleLProtonMuon = std::cos(LeadingProtonVector.Angle(MuonVector));
    fCosOpeningAngleRProtonMuon = std::cos(RecoilProtonVector.Angle(MuonVector));    
    
    // Transverse vectors
    TVector3 MuonVectorTrans(MuonVector.X(), MuonVector.Y(), 0);
    TVector3 TotalProtonVectorTrans(TotalProtonVector.X(), TotalProtonVector.Y(), 0);
    TVector3 DeltaPVectorTrans = MuonVectorTrans + TotalProtonVectorTrans;

    // For GKI calculations
    double TotalProtonMomentum = TotalProtonVector.Mag();
    double Ecal = (TMath::Sqrt((fMuonMomentum*fMuonMomentum) + (0.106*0.106)) + 
    TMath::Sqrt((TotalProtonMomentum*TotalProtonMomentum) + (0.938272*0.938272)) - 0.938272 + 0.0309);
    TVector3 EcalLongVector(0, 0, Ecal);

    //Momentum transfer vector
    TVector3 QVector = EcalLongVector - MuonVector;

    //Missing momentum vector
    TVector3 PVectorN = MuonVector + TotalProtonVector - EcalLongVector;

    // GKI opening angle
    fCosOpeningAngleMomentumTransferTotalProton = std::cos(QVector.Angle(TotalProtonVector));

    // Transverse momentum variable
    fTransverseMomentum = DeltaPVectorTrans.Mag();

    //Missing momentum variable
    fMissingMomentum = PVectorN.Mag();

    //Proton energies
    double ProtonMass = 0.9383; //in Gev/c^2
    double LProtonEnergy = sqrt(pow(ProtonMass, 2) +  pow(LeadingProtonVector.Mag(),2));
    double RProtonEnergy = sqrt(pow(ProtonMass, 2) +  pow(RecoilProtonVector.Mag(),2));
    
    //Invariant Mass of the Hadronic System variable
    fInvariantMass = sqrt(pow((LProtonEnergy + RProtonEnergy), 2) - pow(TotalProtonVector.Mag(), 2));

    // Transverse momentum angular orientation with respect to transverse muon
    fDeltaAlphaT = TMath::ACos(
        (-MuonVectorTrans).Dot(DeltaPVectorTrans) / 
        (MuonVectorTrans.Mag() * DeltaPVectorTrans.Mag())
    ) * (180. / TMath::Pi()) ;
    if (fDeltaAlphaT > 180.) { fDeltaAlphaT -= 180.; }
	if (fDeltaAlphaT < 0.) { fDeltaAlphaT += 180.; }

    // Angular orientation with respect to momentum transfer vector
    fAlphaThreeD = TMath::ACos(
        (QVector).Dot(PVectorN) / 
        (QVector.Mag() * PVectorN.Mag())
    ) * (180. / TMath::Pi()) ;
    if (fAlphaThreeD > 180.) { fAlphaThreeD -= 180.; }
	if (fAlphaThreeD < 0.) { fAlphaThreeD += 180.; }
}

double TwoPTools::ReturnMuonMomentum() {
    return fMuonMomentum;
}

double TwoPTools::ReturnMuonCosTheta() {
    return fMuonCosTheta;
}

double TwoPTools::ReturnLeadingProtonMomentum() {
    return fLeadingProtonMomentum;
}

double TwoPTools::ReturnLeadingProtonCosTheta() {
    return fLeadingProtonCosTheta;
}

double TwoPTools::ReturnRecoilProtonMomentum() {
    return fRecoilProtonMomentum;
}

double TwoPTools::ReturnRecoilProtonCosTheta() {
    return fRecoilProtonCosTheta;
}

double TwoPTools::ReturnCosOpeningAngleProtons() {
    return fCosOpeningAngleProtons;
}

double TwoPTools::ReturnCosOpeningAngleMuonTotalProton() {
    return fCosOpeningAngleMuonTotalProton;
}

double TwoPTools::ReturnTransverseMomentum() {
    return fTransverseMomentum;
}

double TwoPTools::ReturnDeltaAlphaT() {
    return fDeltaAlphaT;
}

double TwoPTools::ReturnMissingMomentum() {
    return fMissingMomentum;
}

double TwoPTools::ReturnInvariantMass(){
    return fInvariantMass;
}

double TwoPTools::ReturnCosOpeningAngleLProtonMuon(){
    return fCosOpeningAngleLProtonMuon;
}

double TwoPTools::ReturnCosOpeningAngleRProtonMuon(){
    return fCosOpeningAngleRProtonMuon;
}

double TwoPTools::ReturnAlphaThreeD() {
    return fAlphaThreeD;
}

double TwoPTools::ReturnCosOpeningAngleMomentumTransferTotalProton() {
    return fCosOpeningAngleMomentumTransferTotalProton;
}

#endif
