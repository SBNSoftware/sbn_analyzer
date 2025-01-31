#ifndef Definitions_h
#define Definitions_h

// SBNAna includes.
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrack.h"
#include "sbnanaobj/StandardRecord/SRVector3D.h"

// Root includes.
#include "TMath.h"
#include "TVector3.h"

// std includes.
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

// Utils includes.
#include "../../Utils/CaloEhadTools.cxx"
#include "../../Utils/Constants.h"
#include "../../Utils/Tools.cxx"
#include "../../Utils/TwoPTools.cxx"

using namespace Constants;

namespace ana {
  Tools tools;

  // Files with samples
  const std::string TargetPath = "/pnfs/sbnd/persistent/users/twester/sbnd/v09_78_04/cv";
  const std::vector<std::string> InputFiles = tools.GetInputFiles(TargetPath);

  // Constants
  const float fFVXMax = 180.f;
  const float fFVXMin = 5.f;
  const float fFVYMax = 180.f;
  const float fFVYMin = 0.f;
  const float fFVZMax = 450.f;
  const float fFVZMin = 10.f;

  const float fMuCutMuScore = 30.0f;
  const float fMuCutPrScore = 60.0f;
  const float fMuCutLength = 50.0f;
  const float fPrCutPrScore = 100.0f;

  const std::map<int, std::tuple<float, float>> PDGToThreshold = {
    {13, {0.1f, 8.f}},                                 // Muon
    {2212, {0.3f, 1.0f}},                               // Proton
    {211, {0.07f, std::numeric_limits<float>::max()}},  // Pi plus
    {-211, {0.07f, std::numeric_limits<float>::max()}}, // Pi minus
    {111, {0.0f, std::numeric_limits<float>::max()}}    // Pi zero
  };

  bool debug=true;//false;

  ///////////
  // Binning
  ///////////

  // Create the binning schemes for the Vars we wish to plot.
  const Binning bEventCount = Binning::Custom(ArrayNBinsEventCount);
  const Binning bVertexX = Binning::Custom(ArrayNBinsVertexX);
  const Binning bVertexY = Binning::Custom(ArrayNBinsVertexY);
  const Binning bVertexZ = Binning::Custom(ArrayNBinsVertexZ);
  const Binning bAngleBins = Binning::Custom(ArrayNBinsAngle);
  const Binning bDeltaAlphaBins = Binning::Custom(ArrayNBinsDeltaAlphaT);
  const Binning bTransverseMomentumBins = Binning::Custom(ArrayNBinsTransverseMomentum);
  const Binning bMuonMomentumBins = Binning::Custom(ArrayNBinsMuonMomentum);
  const Binning bLeadingProtonMomentumBins = Binning::Custom(ArrayNBinsLeadingProtonMomentum);
  const Binning bRecoilProtonMomentumBins = Binning::Custom(ArrayNBinsRecoilProtonMomentum);
  const Binning bInvariantMassBins = Binning::Custom(ArrayNBinsInvariantMass);
  const Binning bCosOpeningAngleLProtonMuonBins = Binning::Custom(ArrayNBinsCosAngleLPMu);
  const Binning bCosOpeningAngleRProtonMuonBins = Binning::Custom(ArrayNBinsCosAngleRPMu);

  // GKI
  const Binning bAlphaThreeDBins = Binning::Custom(ArrayNBinsAlphaThreeD);
  const Binning bMissingMomentumBins = Binning::Custom(ArrayNBinsMissingMomentum);

  // Bins for cut plots
  const Binning bNuScore = Binning::Simple(30, 0, 1.0);
  const Binning bFMatchScore = Binning::Simple(40, 0, 40.0);
  const Binning bFMatchTime = Binning::Simple(20, 0, 5.0);
  const Binning bMuChi2 = Binning::Simple(30, 0, 60.0);
  const Binning bProtonChi2 = Binning::Simple(100, 0, 200.0);

  // Resolution bins
  const Binning bActualResolution = Binning::Simple(51, -1., 1.);

  // Double differential bins
  const Binning bTransverseMomentumInMuonCosTheta = Binning::Custom(
								    tools.Return2DBinIndices(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices));
  const Binning bDeltaAlphaTInMuonCosTheta = Binning::Custom(
							     tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices));
  const Binning bCosOpeningAngleProtonsInMuonCosTheta = Binning::Custom(
									tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices));
  const Binning bCosOpeningAngleMuonTotalProtonInMuonCosTheta = Binning::Custom(
										tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices));
  // GKI
  const Binning bMissingMomentumInMuonCosTheta = Binning::Custom(
								 tools.Return2DBinIndices(TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices));
  const Binning bAlphaThreeDInMuonCosTheta = Binning::Custom(
							     tools.Return2DBinIndices(TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices));
  const Binning bCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta = Binning::Custom(tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices));

  //////////////
  // Functions
  //////////////

  // Check vector is in fiducial volume
  bool bIsInFV(const caf::Proxy<caf::SRVector3D> *data) {
    if (std::isnan(data->x) || std::isnan(data->y) || std::isnan(data->z))
      return false;
    return (
	    (TMath::Abs(data->x) > fFVXMin) && (TMath::Abs(data->x) < fFVXMax) &&
	    (TMath::Abs(data->y) > fFVYMin) && (TMath::Abs(data->y) < fFVYMax) &&
	    (data->z > fFVZMin) && (data->z < fFVZMax));
  }

  // Check multiplicity of particles in a given pdg
  int iCountMultParticle(const caf::SRTrueInteractionProxy *nu, int pdg, float lb, float up) {
    int count = 0;
    for (auto const &prim : nu->prim) {
      float totp = std::sqrt(std::pow(prim.genp.x, 2) + std::pow(prim.genp.y, 2) + std::pow(prim.genp.z, 2));
      if (prim.pdg == pdg && totp >= lb && totp < up) {
	count++;
      }
    }
    return count;
  }

  // Looks for a muon track and gets its PID
  std::tuple<bool, int> bOneMuon(const caf::SRSliceProxy *slc) {
    // Momentum range for muons
    float lb = std::get<0>(PDGToThreshold.at(13));
    float ub = std::get<1>(PDGToThreshold.at(13));

    std::vector<int> CandidateMuons;
    std::vector<int> CandidateMuonsTrkLen;

    for (auto const &pfp : slc->reco.pfp) {
      bool bSkipPFP = false;
      float fMuAverage = 0.0f;
      float fPrAverage = 0.0f;
      for (int i = 0; i < 3; i++) {
	// Skip events with Nan's or 0's in chi squared values
	if (
	    std::isnan(pfp.trk.chi2pid[i].chi2_muon) ||
	    pfp.trk.chi2pid[i].chi2_muon == 0. ||
	    std::isnan(pfp.trk.chi2pid[i].chi2_proton) ||
	    pfp.trk.chi2pid[i].chi2_proton == 0.) {
	  bSkipPFP = true;
	  break;
	}
	fMuAverage += pfp.trk.chi2pid[i].chi2_muon / 3;
	fPrAverage += pfp.trk.chi2pid[i].chi2_proton / 3;
      }
      if (bSkipPFP)
	continue;

      // Check start point is in FV and assign momentum based on end point
      if (!bIsInFV(&pfp.trk.start))
	continue;

      float fMomentum;
      if (!bIsInFV(&pfp.trk.end)) {
	if (std::isnan(pfp.trk.mcsP.fwdP_muon))
	  continue;
	fMomentum = pfp.trk.mcsP.fwdP_muon;
      } else {
	fMomentum = pfp.trk.rangeP.p_muon;
      }

      if (std::isnan(pfp.trk.len))
	continue;
      if (
	  fMuAverage < fMuCutMuScore &&
	  fPrAverage > fMuCutPrScore &&
	  pfp.trk.len > fMuCutLength &&
	  fMomentum >= lb &&
	  fMomentum < ub) {
	CandidateMuons.push_back(pfp.id);
	CandidateMuonsTrkLen.push_back(pfp.trk.len);
      }
    }

    // Choose candidate muon with longest track
    if (CandidateMuons.size() == 0)
      return {false, -1};

    float fLongestTrack = CandidateMuonsTrkLen.at(0);
    int iMuonIndex = CandidateMuons.at(0);
    for (std::size_t i = 1; i < CandidateMuons.size(); i++) {
      if (CandidateMuonsTrkLen.at(i) > fLongestTrack) {
	fLongestTrack = CandidateMuonsTrkLen.at(i);
	iMuonIndex = CandidateMuons.at(i);
      }
    }
    return {true, iMuonIndex};
  }

  // Look for protons and get their PIDs
  std::tuple<bool, std::vector<int>> bTwoProtons(const caf::SRSliceProxy *slc, int MuonID) {
    // Momentum range for protons
    float lb = std::get<0>(PDGToThreshold.at(2212));
    float ub = std::get<1>(PDGToThreshold.at(2212));

    std::vector<int> ProtonIDs;
    for (auto const &pfp : slc->reco.pfp) {
      if (pfp.id == MuonID)
	continue; // skip pfp tagged as muon
      
      float fPrAverage = 0.0f;
      bool bSkipPFP = false;

      for (int i = 0; i < 3; i++) {
	// Skip events with Nan's or 0's in chi squared values
	if (
	    std::isnan(pfp.trk.chi2pid[i].chi2_muon) ||
	    pfp.trk.chi2pid[i].chi2_muon == 0.) {
	  bSkipPFP = true;
	  break;
	}
      	fPrAverage += pfp.trk.chi2pid[i].chi2_proton / 3;
      }
      if (bSkipPFP)
	continue;

      // Check full track is in FV
      if (!(bIsInFV(&pfp.trk.start) && bIsInFV(&pfp.trk.end)))
	continue;
      float fMomentum = pfp.trk.rangeP.p_proton;

      if (
	  fPrAverage < fPrCutPrScore &&
	  fMomentum >= lb &&
	  fMomentum < ub)
	ProtonIDs.push_back(pfp.id);
    }
    return {ProtonIDs.size() == 2, ProtonIDs};
  }

  bool bNoChargedPions(const caf::SRSliceProxy *slc, std::vector<int> TaggedIDs) {
    // Momentum range for charged pions
    float lb = std::get<0>(PDGToThreshold.at(211));
    float ub = std::get<1>(PDGToThreshold.at(211));

    for (auto const &pfp : slc->reco.pfp) {
      if (std::find(TaggedIDs.begin(), TaggedIDs.end(), pfp.id) != TaggedIDs.end())
	continue;
      if (!(bIsInFV(&pfp.trk.start) && bIsInFV(&pfp.trk.end)))
	continue;

      float fMomentum = pfp.trk.rangeP.p_pion;
      if (fMomentum >= lb && fMomentum < ub)
	return false; // tag pion
    }
    return true;
  }

  bool bNoShowers(const caf::SRSliceProxy *slc, std::vector<int> TaggedIDs) {
    for (auto const &pfp : slc->reco.pfp) {
      if (std::find(TaggedIDs.begin(), TaggedIDs.end(), pfp.id) != TaggedIDs.end())
	continue;
      if (pfp.trackScore == -5.f)
	continue; // clear cosmic
      if (pfp.trackScore < 0.5)
	return false;
    }
    return true;
  }

  // Gets reconstructed momentum vectors given the particle IDs
  std::tuple<TVector3, TVector3, TVector3> GetVectors(const caf::SRSliceProxy *slc, int MuonID, int ProtonID1, int ProtonID2) {
    TVector3 Muon(1, 1, 1);
    TVector3 LeadingProton(1, 1, 1);
    TVector3 RecoilProton(1, 1, 1);

    for (auto const &pfp : slc->reco.pfp) {
      if (pfp.id == MuonID) {
	Muon.SetTheta(TMath::ACos(pfp.trk.costh));
	Muon.SetPhi(pfp.trk.phi);
	if (bIsInFV(&pfp.trk.end)) {
	  Muon.SetMag(pfp.trk.rangeP.p_muon);
	} else {
	  Muon.SetMag(pfp.trk.mcsP.fwdP_muon);
	}
      } else if (pfp.id == ProtonID1) {
	LeadingProton.SetTheta(TMath::ACos(pfp.trk.costh));
	LeadingProton.SetPhi(pfp.trk.phi);
	LeadingProton.SetMag(pfp.trk.rangeP.p_proton);
      } else if (pfp.id == ProtonID2) {
	RecoilProton.SetTheta(TMath::ACos(pfp.trk.costh));
	RecoilProton.SetPhi(pfp.trk.phi);
	RecoilProton.SetMag(pfp.trk.rangeP.p_proton);
      }
    }
    return {Muon, LeadingProton, RecoilProton};
  }

  TVector3 GetParticleVector(const caf::SRSliceProxy *slc, int PID, int partPDG) {
    TVector3 partVec(-999, -999, -999);
    // This is really only gonna work for track based particles and needs new caveats for like electrons etc
    for (auto const &pfp : slc->reco.pfp) {
      if (pfp.id == PID) {
	//If its a track
        //if(pfp.trackScore>0.5){//this def needs to be optimized... 

	  partVec.SetTheta(TMath::ACos(pfp.trk.costh));
	  partVec.SetPhi(pfp.trk.phi);
	  if (abs(partPDG) == 13) {
            if (bIsInFV(&pfp.trk.end)) {
	      partVec.SetMag(pfp.trk.rangeP.p_muon);
            } else {
	      partVec.SetMag(pfp.trk.mcsP.fwdP_muon);
            }
	  } else if (abs(partPDG) == 2212) {
            partVec.SetMag(pfp.trk.rangeP.p_proton);
	  } else if (abs(partPDG) == 211) {
            partVec.SetMag(pfp.trk.rangeP.p_pion);
	  } else {
            cout << "This PDG doesn't have a well defined momentum, using rangeP.p_muon. PDG: " << partPDG << endl;
            partVec.SetMag(pfp.trk.rangeP.p_muon);
	  }
        }
	/*      }
      else{//shower like -- I think we really need a case for looking like a cosmic though?
	pfp.shw.bestplane_energy
	}*/
      return partVec;
    }//end loop over pfps 
    cout<<"PFP with PID ="<< PID<<"not found."<<endl;
    return partVec;//return -999 vec if i didnt find the PID
  }

  double GetShowerEnergy(const caf::SRSliceProxy *slc, int PID){
    double energy=-999.;
    for (auto const &pfp: slc->reco.pfp){
      if(pfp.id==PID){
	energy=pfp.shw.bestplane_energy;
	//seriously don't know about this
	break;
      }
    }
      
    return energy;
  }

  // Get true momentum vectors
  std::tuple<TVector3, vector<TVector3>, vector<int>> GetTrueVector(const caf::SRTrueInteractionProxy *nu) {
    if(debug) cout<<"entering GetTrueVector"<<endl;
    TVector3 Muon(1, 1, 1);
    vector<TVector3> HadronVectors; // LeadingProton(1, 1, 1);
    vector<int> HadronPIDs;
    // TVector3 RecoilProton(1, 1, 1);
    bool FirstMuon = false;
    bool FirstProton = false;
    bool SecondProton = false;
    if(debug) cout<<"GetTrueVector:l342"<<endl;   
    for (auto const &prim : nu->prim) {
      if(debug) cout<<"GetTrueVector:l344"<<endl;   
      float totp = std::sqrt(std::pow(prim.genp.x, 2) + std::pow(prim.genp.y, 2) + std::pow(prim.genp.z, 2));
      int pdg = prim.pdg;
      if(debug) cout<<"GetTrueVector:l347 pdg: "<<prim.pdg<<"\t totp: "<<totp<<endl;   
      if (
	  (prim.pdg == 13) &&
	  (totp >= std::get<0>(PDGToThreshold.at(13))) &&
	  (totp < std::get<1>(PDGToThreshold.at(13)))) {

	    Muon.SetXYZ(prim.genp.x, prim.genp.y, prim.genp.z);
	    FirstMuon = true;
      } else if (
		 (prim.pdg == 2212) &&
		 (totp >= std::get<0>(PDGToThreshold.at(2212))) &&
		 (totp < std::get<1>(PDGToThreshold.at(2212)))) {
                  if(debug) cout<<"GetTrueVector:l359"<<endl;   

	    HadronVectors.push_back({1, 1, 1});
	    HadronVectors.back().SetXYZ(prim.genp.x, prim.genp.y, prim.genp.z);
	    FirstProton = true;
	    HadronPIDs.push_back(pdg);
      } else  {//need some sort of if for other pdgs but that involves setting a threshold for every possible particle type? or i just need to make a blanket statement for other particles or something
            if(debug) cout<<"GetTrueVector:l369"<<endl;   
	    HadronVectors.push_back({1, 1, 1});
	    HadronVectors.back().SetXYZ(prim.genp.x, prim.genp.y, prim.genp.z);
	    HadronPIDs.push_back(pdg);
      }
            if(debug) cout<<"GetTrueVector:l370"<<endl;    
    }

    if(debug) cout<<"l371"<<endl;
    if (!(FirstMuon)) {
      std::cout << "GetTrueVector: Didn't find a muon" << std::endl;
      exit(-1);
    }
    if(debug) cout<<"Returning TrueVector"<<endl;

    return {Muon, HadronVectors, HadronPIDs};
  }

  // 2p version
  std::tuple<TVector3, TVector3, TVector3> GetTrueVectorCC2p(const caf::SRTrueInteractionProxy *nu) {
    TVector3 Muon(1, 1, 1);
    TVector3 LeadingProton(1, 1, 1);
    TVector3 RecoilProton(1, 1, 1);
    bool FirstMuon = false;
    bool FirstProton = false;
    bool SecondProton = false;

    for (auto const &prim : nu->prim) {
      float totp = std::sqrt(std::pow(prim.genp.x, 2) + std::pow(prim.genp.y, 2) + std::pow(prim.genp.z, 2));
      if ( (prim.pdg == 13) &&  (totp >= std::get<0>(PDGToThreshold.at(13))) && (totp < std::get<1>(PDGToThreshold.at(13)))) {
        Muon.SetXYZ(prim.genp.x, prim.genp.y, prim.genp.z);
        FirstMuon = true;
      } else if ( (prim.pdg == 2212) && (totp >= std::get<0>(PDGToThreshold.at(2212))) && (totp < std::get<1>(PDGToThreshold.at(2212))) && !FirstProton) {
        LeadingProton.SetXYZ(prim.genp.x, prim.genp.y, prim.genp.z);
        FirstProton = true;
      } else if ( (prim.pdg == 2212) && (totp >= std::get<0>(PDGToThreshold.at(2212))) && (totp < std::get<1>(PDGToThreshold.at(2212)))) {
        RecoilProton.SetXYZ(prim.genp.x, prim.genp.y, prim.genp.z);
        SecondProton = true;
      }
    }

    if (!(FirstMuon && FirstProton && SecondProton)) {
      std::cout << "GetTrueVectorCCp: All particles not found" << std::endl;
      exit(-1);
    }

    return {Muon, LeadingProton, RecoilProton};
  }

  // Get all vars from TwoPTools helper
  std::vector<double> GetVarsFromHelper(TwoPTools Helper) {
    std::vector<double> vars;

    double MuonCosTheta = Helper.ReturnMuonCosTheta();
    vars.push_back(MuonCosTheta);

    double LeadingProtonCosTheta = Helper.ReturnLeadingProtonCosTheta();
    vars.push_back(LeadingProtonCosTheta);

    double RecoilProtonCosTheta = Helper.ReturnRecoilProtonCosTheta();
    vars.push_back(RecoilProtonCosTheta);

    double CosOpeningAngleProtons = Helper.ReturnCosOpeningAngleProtons();
    vars.push_back(CosOpeningAngleProtons);

    double CosOpeningAngleMuonTotalProton = Helper.ReturnCosOpeningAngleMuonTotalProton();
    vars.push_back(CosOpeningAngleMuonTotalProton);

    double DeltaAlphaT = Helper.ReturnDeltaAlphaT();
    vars.push_back(DeltaAlphaT);

    double TransverseMomentum = Helper.ReturnTransverseMomentum();
    vars.push_back(TransverseMomentum);

    double MuonMomentum = Helper.ReturnMuonMomentum();
    vars.push_back(MuonMomentum);

    double LeadingProtonMomentum = Helper.ReturnLeadingProtonMomentum();
    vars.push_back(LeadingProtonMomentum);

    double RecoilProtonMomentum = Helper.ReturnRecoilProtonMomentum();
    vars.push_back(RecoilProtonMomentum);

    // GKI
    double CosOpeningAngleMomentumTransferTotalProton = Helper.ReturnCosOpeningAngleMomentumTransferTotalProton();
    vars.push_back(CosOpeningAngleMomentumTransferTotalProton);

    double AlphaThreeD = Helper.ReturnAlphaThreeD();
    vars.push_back(AlphaThreeD);

    double MissingMomentum = Helper.ReturnMissingMomentum();
    vars.push_back(MissingMomentum);

    // Additional variables
    double InvariantMass = Helper.ReturnInvariantMass();
    vars.push_back(InvariantMass);

    double CosOpeningAngleLProtonMuon = Helper.ReturnCosOpeningAngleLProtonMuon();
    vars.push_back(CosOpeningAngleLProtonMuon);

    double CosOpeningAngleRProtonMuon = Helper.ReturnCosOpeningAngleRProtonMuon();
    vars.push_back(CosOpeningAngleRProtonMuon);

    return vars;
  }

  // Get all vars from TwoPTools helper
  /*std::vector<double> GetVarsFromHelper(CaloEhadTools Helper)
    {
    std::vector<double> vars;

    double MuonCosTheta = Helper.ReturnMuonCosTheta();
    vars.push_back(MuonCosTheta);

    double LeadingProtonCosTheta = Helper.ReturnLeadingProtonCosTheta();
    vars.push_back(LeadingProtonCosTheta);

    double RecoilProtonCosTheta = Helper.ReturnRecoilProtonCosTheta();
    vars.push_back(RecoilProtonCosTheta);

    double CosOpeningAngleProtons = Helper.ReturnCosOpeningAngleProtons();
    vars.push_back(CosOpeningAngleProtons);

    double CosOpeningAngleMuonTotalProton = Helper.ReturnCosOpeningAngleMuonTotalProton();
    vars.push_back(CosOpeningAngleMuonTotalProton);

    double DeltaAlphaT = Helper.ReturnDeltaAlphaT();
    vars.push_back(DeltaAlphaT);

    double TransverseMomentum = Helper.ReturnTransverseMomentum();
    vars.push_back(TransverseMomentum);

    double MuonMomentum = Helper.ReturnMuonMomentum();
    vars.push_back(MuonMomentum);

    double LeadingProtonMomentum = Helper.ReturnLeadingProtonMomentum();
    vars.push_back(LeadingProtonMomentum);

    double RecoilProtonMomentum = Helper.ReturnRecoilProtonMomentum();
    vars.push_back(RecoilProtonMomentum);

    // GKI
    double CosOpeningAngleMomentumTransferTotalProton = Helper.ReturnCosOpeningAngleMomentumTransferTotalProton();
    vars.push_back(CosOpeningAngleMomentumTransferTotalProton);

    double AlphaThreeD = Helper.ReturnAlphaThreeD();
    vars.push_back(AlphaThreeD);

    double MissingMomentum = Helper.ReturnMissingMomentum();
    vars.push_back(MissingMomentum);

    // Additional variables
    double InvariantMass = Helper.ReturnInvariantMass();
    vars.push_back(InvariantMass);

    double CosOpeningAngleLProtonMuon = Helper.ReturnCosOpeningAngleLProtonMuon();
    vars.push_back(CosOpeningAngleLProtonMuon);

    double CosOpeningAngleRProtonMuon = Helper.ReturnCosOpeningAngleRProtonMuon();
    vars.push_back(CosOpeningAngleRProtonMuon);

    return vars;
    }*/

  ////////////
  // MultiVars
  ////////////

  const TruthMultiVar kTruthVars([](const caf::SRTrueInteractionProxy *nu) {
      // Get vectors
      if(debug) cout<<"Getting kTruthVars"<<endl;
      auto [Muon, LeadingProton, RecoilProton] = GetTrueVectorCC2p(nu);

      // Extract vars from helper
      TwoPTools Helper(Muon, LeadingProton, RecoilProton);
      return GetVarsFromHelper(Helper); });

  // Contains all variables we are interested in
  const MultiVar kVars([](const caf::SRSliceProxy *slc) -> std::vector<double> {
      // Get IDs for tagged particles
      std::vector<int> TaggedIDs;
      auto [OneMuon, MuonID] = bOneMuon(slc);
      auto [TwoProtons, ProtonIDs] = bTwoProtons(slc, MuonID);
      auto [Muon, LeadingProton, RecoilProton] = GetVectors(slc, MuonID, ProtonIDs.at(0), ProtonIDs.at(1));

      // Extract vars from helper
      TwoPTools Helper(Muon, LeadingProton, RecoilProton);
      std::vector<double> vars = GetVarsFromHelper(Helper);
      return vars;
    });

  //////////////
  // Vars
  //////////////

  // Dummy variables to keep track of events
  const Var kEventCount([](const caf::SRSliceProxy *slc) -> double { return 0.5; });
  const TruthVar kTrueEventCount([](const caf::SRTrueInteractionProxy *nu) -> double { return 0.5; });

  // Cosmic cut variables
  const Var kNuScore([](const caf::SRSliceProxy *slc) -> double { return slc->nu_score; });
  const Var kFMatchScore([](const caf::SRSliceProxy *slc) -> double { return slc->fmatch.score; });
  const Var kFMatchTime([](const caf::SRSliceProxy *slc) -> double { return slc->fmatch.time; });

  // Mu chi2 for muon
  const Var kMuMuChi2([](const caf::SRSliceProxy *slc) -> double {
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 13) {
	  bool bSkipPFP = false;
	  double fMuAverage = 0.;
	  for (int i = 0; i < 3; i++) { 
	    if (
		std::isnan(pfp.trk.chi2pid[i].chi2_muon) ||
		pfp.trk.chi2pid[i].chi2_muon == 0.
		) {
	      bSkipPFP = true;
	      break;
	    } else {
	      fMuAverage += pfp.trk.chi2pid[i].chi2_muon / 3;
	    }
	  }
	  if (!bSkipPFP) return fMuAverage;
	}
      }
      return 0; });

  // Proton chi2 for muon
  const Var kMuProtonChi2([](const caf::SRSliceProxy *slc) -> double {
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 13) {
	  bool bSkipPFP = false;
	  double fPrAverage = 0.;
	  for (int i = 0; i < 3; i++) { 
	    if (
		std::isnan(pfp.trk.chi2pid[i].chi2_proton) ||
		pfp.trk.chi2pid[i].chi2_proton == 0.
		) {
	      bSkipPFP = true;
	      break;
	    } else {
	      fPrAverage += pfp.trk.chi2pid[i].chi2_proton / 3;
	    }
	  }
	  if (!bSkipPFP) return fPrAverage;
	}
      }
      return 0; });

  // Mu chi2 for proton
  const Var kProtonMuChi2([](const caf::SRSliceProxy *slc) -> double {
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 2212) {
	  bool bSkipPFP = false;
	  double fMuAverage = 0.;
	  for (int i = 0; i < 3; i++) { 
	    if (
		std::isnan(pfp.trk.chi2pid[i].chi2_muon) ||
		pfp.trk.chi2pid[i].chi2_muon == 0.
		) {
	      bSkipPFP = true;
	      break;
	    } else {
	      fMuAverage += pfp.trk.chi2pid[i].chi2_muon / 3;
	    }
	  }
	  if (!bSkipPFP) return fMuAverage;
	}
      }
      return 0; });

  // Proton chi2 for proton
  const Var kProtonProtonChi2([](const caf::SRSliceProxy *slc) -> double {
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 2212) {
	  bool bSkipPFP = false;
	  double fPrAverage = 0.;
	  for (int i = 0; i < 3; i++) { 
	    if (
		std::isnan(pfp.trk.chi2pid[i].chi2_proton) ||
		pfp.trk.chi2pid[i].chi2_proton == 0.
		) {
	      bSkipPFP = true;
	      break;
	    } else {
	      fPrAverage += pfp.trk.chi2pid[i].chi2_proton / 3;
	    }
	  }
	  if (!bSkipPFP) return fPrAverage;
	}
      }
      return 0; });

  // Mu chi2 for second proton
  const Var kSecondProtonMuChi2([](const caf::SRSliceProxy *slc) -> double {
      bool firstProton = false;
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 2212) {
	  bool bSkipPFP = false;
	  double fMuAverage = 0.;
	  for (int i = 0; i < 3; i++) { 
	    if (
		std::isnan(pfp.trk.chi2pid[i].chi2_muon) ||
		pfp.trk.chi2pid[i].chi2_muon == 0.
		) {
	      bSkipPFP = true;
	      break;
	    } else {
	      fMuAverage += pfp.trk.chi2pid[i].chi2_muon / 3;
	    }
	  }
	  if (!firstProton && !bSkipPFP) firstProton = true;
	  else if (firstProton && !bSkipPFP) return fMuAverage;
	}
      }
      return 0; });

  // Proton chi2 for second proton
  const Var kSecondProtonProtonChi2([](const caf::SRSliceProxy *slc) -> double {
      bool firstProton = false;
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 2212) {
	  bool bSkipPFP = false;
	  double fPrAverage = 0.;
	  for (int i = 0; i < 3; i++) { 
	    if (
		std::isnan(pfp.trk.chi2pid[i].chi2_proton) ||
		pfp.trk.chi2pid[i].chi2_proton == 0.
		) {
	      bSkipPFP = true;
	      break;
	    } else {
	      fPrAverage += pfp.trk.chi2pid[i].chi2_proton / 3;
	    }
	  }
	  if (!firstProton && !bSkipPFP) firstProton = true;
	  else if (firstProton && !bSkipPFP) return fPrAverage;
	}
      }
      return 0; });

  // Mu chi2 for pion
  const Var kPionMuChi2([](const caf::SRSliceProxy *slc) -> double {
      bool firstProton = false;
      for (auto const &pfp : slc->reco.pfp) {
        if (pfp.trk.truth.p.pdg == 211 || pfp.trk.truth.p.pdg == -211) {
	  bool bSkipPFP = false;
	  double fMuAverage = 0.;
	  for (int i = 0; i < 3; i++) {
	    if (
		std::isnan(pfp.trk.chi2pid[i].chi2_muon) ||
		pfp.trk.chi2pid[i].chi2_muon == 0.) {
	      bSkipPFP = true;
	      break;
	    } else {
	      fMuAverage += pfp.trk.chi2pid[i].chi2_muon / 3;
	    }
	  }
	  if (!firstProton && !bSkipPFP)
	    firstProton = true;
	  else if (firstProton && !bSkipPFP)
	    return fMuAverage;
        }
      }
      return 0;
    });

  // Proton chi2 for pion
  const Var kPionProtonChi2([](const caf::SRSliceProxy *slc) -> double {
      bool firstProton = false;
      for (auto const &pfp : slc->reco.pfp) {
        if (pfp.trk.truth.p.pdg == 211 || pfp.trk.truth.p.pdg == -211) {
	  bool bSkipPFP = false;
	  double fPrAverage = 0.;
	  for (int i = 0; i < 3; i++) {
	    if (
		std::isnan(pfp.trk.chi2pid[i].chi2_proton) ||
		pfp.trk.chi2pid[i].chi2_proton == 0.) {
	      bSkipPFP = true;
	      break;
	    } else {
	      fPrAverage += pfp.trk.chi2pid[i].chi2_proton / 3;
	    }
	  }
	  if (!firstProton && !bSkipPFP)
	    firstProton = true;
	  else if (firstProton && !bSkipPFP)
	    return fPrAverage;
        }
      }
      return 0;
    });

  // For all the following vars, there are three variations:
  //     1. Reconstructed value `Var`
  //     2. True value `TruthVar`
  //     3. True value accessed through truth information of
  //        a `SRSlice`, and stored as `Var`
  // The difference between 1 and 3 is that the latter can be
  // used in a `Spectrum` that uses a `Cut` instead of a `TruthCut`.

  // Vertex distribution
  const Var kVertexX([](const caf::SRSliceProxy *slc) -> float {
      if (std::isnan(slc->vertex.x))
        return 160.;
      return slc->vertex.x;
    });

  const TruthVar kTruthVertexX([](const caf::SRTrueInteractionProxy *nu) -> float {
      if (std::isnan(nu->position.x))
        return 160.;
      return nu->position.x;
    });

  const Var kRecoTruthVertexX([](const caf::SRSliceProxy *slc) -> float {
      return kTruthVertexX(&slc->truth);
    });

  const Var kVertexY([](const caf::SRSliceProxy *slc) -> float {
      if (std::isnan(slc->vertex.y))
        return 0.;
      return slc->vertex.y;
    });

  const TruthVar kTruthVertexY([](const caf::SRTrueInteractionProxy *nu) -> float {
      if (std::isnan(nu->position.y)) return 0.;
      return nu->position.y; });

  const Var kRecoTruthVertexY([](const caf::SRSliceProxy *slc) -> float { return kTruthVertexY(&slc->truth); });

  const Var kVertexZ([](const caf::SRSliceProxy *slc) -> float {
      if (std::isnan(slc->vertex.z)) return 230.;
      return slc->vertex.z; });

  const TruthVar kTruthVertexZ([](const caf::SRTrueInteractionProxy *nu) -> float {
      if (std::isnan(nu->position.z)) return 230.;
      return nu->position.z; });

  const Var kRecoTruthVertexZ([](const caf::SRSliceProxy *slc) -> float { return kTruthVertexZ(&slc->truth); });

  // Muon angle
  const Var kMuonCosTheta([](const caf::SRSliceProxy *slc) -> double { 
      auto [OneMuon, MuonID] = bOneMuon(slc);
      TVector3 MuonVec=GetParticleVector(slc, MuonID, 13);
      return MuonVec.CosTheta();
    });

  const TruthVar kTruthMuonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double { 
      auto [MuonVec, HadronVector, HadronPDG] = GetTrueVector(nu);
      return MuonVec.CosTheta();
    });

  const Var kRecoTruthMuonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthMuonCosTheta(&slc->truth); });

  // Muon angle
  const Var kMuonTheta([](const caf::SRSliceProxy *slc) -> double { 
      auto [OneMuon, MuonID] = bOneMuon(slc);
      TVector3 MuonVec=GetParticleVector(slc, MuonID, 13);
      return MuonVec.Theta();
    });

  const TruthVar kTruthMuonTheta([](const caf::SRTrueInteractionProxy *nu) -> double { 
      auto [MuonVec, HadronVector, HadronPDG] = GetTrueVector(nu);
      return MuonVec.Theta();
    });

  const Var kRecoTruthMuonTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthMuonTheta(&slc->truth); });


  // Leading proton angle
  const Var kLeadingProtonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(1); });

  const TruthVar kTruthLeadingProtonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(1); });

  const Var kRecoTruthLeadingProtonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthLeadingProtonCosTheta(&slc->truth); });

  // Recoil proton angle
  const Var kRecoilProtonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(2); });

  const TruthVar kTruthRecoilProtonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(2); });

  const Var kRecoTruthRecoilProtonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthRecoilProtonCosTheta(&slc->truth); });

  // Opening angle between protons
  const Var kCosOpeningAngleProtons([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(3); });

  const TruthVar kTruthCosOpeningAngleProtons([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(3); });

  const Var kRecoTruthCosOpeningAngleProtons([](const caf::SRSliceProxy *slc) -> double { return kTruthCosOpeningAngleProtons(&slc->truth); });

  // Opening angle between muon and total proton
  const Var kCosOpeningAngleMuonTotalProton([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(4); });
  const TruthVar kTruthCosOpeningAngleMuonTotalProton([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(4); });
  const Var kRecoTruthCosOpeningAngleMuonTotalProton([](const caf::SRSliceProxy *slc) -> double { return kTruthCosOpeningAngleMuonTotalProton(&slc->truth); });

  // Delta alpha transverse
  const Var kDeltaAlphaT([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(5); });
  const TruthVar kTruthDeltaAlphaT([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(5); });
  const Var kRecoTruthDeltaAlphaT([](const caf::SRSliceProxy *slc) -> double { return kTruthDeltaAlphaT(&slc->truth); });

  // Transverse momentum
  const Var kTransverseMomentum([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(6); });
  const TruthVar kTruthTransverseMomentum([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(6); });
  const Var kRecoTruthTransverseMomentum([](const caf::SRSliceProxy *slc) -> double { return kTruthTransverseMomentum(&slc->truth); });

  // Muon momentum
  const Var kMuonMomentum([](const caf::SRSliceProxy *slc) -> double {
      std::vector<int> TaggedIDs;
      auto [OneMuon, MuonID] = bOneMuon(slc);
      return GetParticleVector(slc, MuonID, 13).Mag();
    });

  const TruthVar kTruthMuonMomentum([](const caf::SRTrueInteractionProxy *nu) -> double {
      if(debug) cout<<"entering kTruthMuonMomentum"<<endl;
      auto [Muon, HadronVector, HadronPDG] = GetTrueVector(nu); 
      if(debug) cout<<"l852"<<endl;
      double muonP=Muon.Mag();
      if(debug) cout<<"l854"<<endl;
      if (muonP > 0) double a=0; 
      else if(debug) cout<<"Muon mag is "<<muonP<<endl;
      return Muon.Mag();
    });

  const Var kRecoTruthMuonMomentum([](const caf::SRSliceProxy *slc) -> double {
      return kTruthMuonMomentum(&slc->truth);
    });


  //Neutrino energy
  const TruthVar kTruthNeutrinoEnergy([](const caf::SRTrueInteractionProxy *nu) -> double {
      if(debug) cout<<"entering kTruthNeutrinoEnergy"<<endl;
      double Enu=nu->E;
      return Enu;
    });


  //Ehad
  const TruthVar kTruthHadronicEnergy([](const caf::SRTrueInteractionProxy *nu) -> double {
      if(debug) cout<<"entering kTruthHadronicEnergy"<<endl;
      double Enu=kTruthNeutrinoEnergy(nu);
      double Emu=std::hypot(kTruthMuonMomentum(nu), M_MU);
      return (Enu-Emu);
    });
  
  const Var kRecoTruthHadronicEnergy([](const caf::SRSliceProxy *slc) -> double{
      return kTruthHadronicEnergy(&slc->truth);
    });

  //EAvaialable
  const Var kEAvail([](const caf::SRSliceProxy *slc) -> double { 
      double eavail=0;
      double p=-999; 
      auto [OneMuon, MuonID] = bOneMuon(slc);
      for (auto const &pfp: slc->reco.pfp){
	if (pfp.id==MuonID) continue;
	//if(pfp.trackScore<0.5 && pfp.shw.bestplane_energy>0){
    if(pfp.shw.bestplane_energy>0)
	  eavail+=pfp.shw.bestplane_energy;
    //lets just look at treating everything as a shower first even though theyre obviously not 
	/*}//end if shower
	else{
	  //GetParticleVector(slc, pfp.id, pfp.pdg);
	  int bestplane=pfp.trk.bestplane;
	  if (bestplane==-1){ cout<<"No calorimetry found for PID "<<pfp.id<<endl;
	    continue;
	  }
	  int partPDG=pfp.trk.chi2pid[bestplane].pdg;//no idea if this is relieabe but minimum particle chi2 seems reasonable enough to me 
	  double caloT=pfp.trk.calo[bestplane].ke*mev_to_gev; //these numbers are huge like 200... must be in MeV? or is really not what I think it is. 
      if(std::is_nan(caloT)){ 
        cout<<"caloT is a nan. partpdg"<<partPDG<<endl;
        continue;
      }
	  eavail+=caloT;
	  //Get a momentum here too 
	  
	  if (abs(partPDG) == 13) {
            if (bIsInFV(&pfp.trk.end)) {
	      p+=pfp.trk.rangeP.p_muon;
            } else {
	      p+=pfp.trk.mcsP.fwdP_muon;
            }
	  } else if (abs(partPDG) == 2212) {
            p+=pfp.trk.rangeP.p_proton;
	  } else if (abs(partPDG) == 211) {
            p+=pfp.trk.rangeP.p_pion;
	  } else {
            cout << "This PDG doesn't have a well defined momentum, using rangeP.p_muon. PDG: " << partPDG << endl;
            p+=pfp.trk.rangeP.p_muon;
	  }
	    cout<<"calo T: "<<caloT<<" p: "<<p<< " PDG: "<<partPDG<< endl;
	}//end if track*/
      }//end loop over pfps
      return eavail; 
    });

  const TruthVar kTruthEAvail([](const caf::SRTrueInteractionProxy *nu) -> double { 
      //auto [Muon, HadronVector, HadronPDG] = GetTrueVector(nu); 
      //return kTruthVars(nu).at(6); 
      return kTruthHadronicEnergy(nu);
    });
  const Var kRecoTruthEAvail([](const caf::SRSliceProxy *slc) -> double { return kTruthEAvail(&slc->truth); });
  

  //Q2
  
  const TruthVar kTruthQ2([](const caf::SRTrueInteractionProxy *nu) -> double {
      double Q2= CalcQ2( kTruthNeutrinoEnergy(nu), kTruthMuonMomentum(nu), kTruthMuonTheta(nu) );
      return Q2;
    });

  const Var kRecoTruthQ2([](const caf::SRSliceProxy *slc) -> double {
      //double Q2= CalcQ2( kRecoTruthNeutrinoEnergy(slc), kRecoTruthMuonEnergy(slc), kRecoTruthMuonTheta(slc) );
      return kTruthQ2(&slc->truth); });
  
  
  //W -- InvariantMass
  const Var kInvariantMass([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(13); });// FIX ME!
  
  const TruthVar kTruthInvariantMass([](const caf::SRTrueInteractionProxy *nu) -> double { 
      //auto [Muon, HadronVector, HadronPDG] = GetTrueVector(nu); 
      double W= CalcW( kTruthNeutrinoEnergy(nu), kTruthMuonMomentum(nu), kTruthMuonTheta(nu) );
      return W;
    });
  
  const Var kRecoTruthInvariantMass([](const caf::SRSliceProxy *slc) -> double { return kTruthInvariantMass(&slc->truth); });
  
  //same thing but calling it W because I know thats what how I'm going to want call this function
  const Var kW([](const caf::SRSliceProxy *slc) -> double { return kInvariantMass(slc);});
  const TruthVar kTruthW([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthInvariantMass(nu);});
  const Var kRecoTruthW([](const caf::SRSliceProxy *slc) -> double { return kTruthInvariantMass(&slc->truth); });



  // Leading proton momentum
  const Var kLeadingProtonMomentum([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(8); });
  const TruthVar kTruthLeadingProtonMomentum([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(8); });
  const Var kRecoTruthLeadingProtonMomentum([](const caf::SRSliceProxy *slc) -> double { return kTruthLeadingProtonMomentum(&slc->truth); });

  // Recoil proton momentum
  const Var kRecoilProtonMomentum([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(9); });
  const TruthVar kTruthRecoilProtonMomentum([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(9); });
  const Var kRecoTruthRecoilProtonMomentum([](const caf::SRSliceProxy *slc) -> double { return kTruthRecoilProtonMomentum(&slc->truth); });

  // Opening angle between momentum transfer and total proton
  const Var kCosOpeningAngleMomentumTransferTotalProton([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(10); });
  const TruthVar kTruthCosOpeningAngleMomentumTransferTotalProton([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(10); });
  const Var kRecoTruthCosOpeningAngleMomentumTransferTotalProton([](const caf::SRSliceProxy *slc) -> double { return kTruthCosOpeningAngleMomentumTransferTotalProton(&slc->truth); });

  // Alpha three dimensional
  const Var kAlphaThreeD([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(11); });
  const TruthVar kTruthAlphaThreeD([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(11); });
  const Var kRecoTruthAlphaThreeD([](const caf::SRSliceProxy *slc) -> double { return kTruthAlphaThreeD(&slc->truth); });

  // Missing momentum
  const Var kMissingMomentum([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(12); });
  const TruthVar kTruthMissingMomentum([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(12); });
  const Var kRecoTruthMissingMomentum([](const caf::SRSliceProxy *slc) -> double { return kTruthMissingMomentum(&slc->truth); });


  // Opening Angle between Leading Proton and Muon
  const Var kCosOpeningAngleLProtonMuon([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(14); });

  const TruthVar kTruthCosOpeningAngleLProtonMuon([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(14); });
  const Var kRecoTruthCosOpeningAngleLProtonMuon([](const caf::SRSliceProxy *slc) -> double { return kTruthCosOpeningAngleLProtonMuon(&slc->truth); });

  // Opening Angle between Recoil Proton and Muon
  const Var kCosOpeningAngleRProtonMuon([](const caf::SRSliceProxy *slc) -> double { return kVars(slc).at(15); });
  const TruthVar kTruthCosOpeningAngleRProtonMuon([](const caf::SRTrueInteractionProxy *nu) -> double { return kTruthVars(nu).at(15); });
  const Var kRecoTruthCosOpeningAngleRProtonMuon([](const caf::SRSliceProxy *slc) -> double { return kTruthCosOpeningAngleRProtonMuon(&slc->truth); });

  ////////////////////////////////
  // Double differential variables
  ////////////////////////////////

  // All variables are in muon cos theta

  // Transverse momentum
  const Var kTransverseMomentumInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
      float fTransverseMomentum = kTransverseMomentum(slc);

      if (fTransverseMomentum < TwoDArrayTransverseMomentum[0]) {
        fTransverseMomentum = (TwoDArrayTransverseMomentum[0] + TwoDArrayTransverseMomentum[1]) / 2.;
      } else if (fTransverseMomentum > TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum]) {
        fTransverseMomentum = (TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum] + TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum - 1]) / 2.;
      }

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kMuonCosTheta(slc), TwoDArrayNBinsMuonCosTheta);
      int SerialTransverseMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
										  TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices,
										  MuonCosThetaTwoDIndex,
										  fTransverseMomentum);
      return SerialTransverseMomentumInMuonCosThetaIndex;
    });


  const TruthVar kTruthTransverseMomentumInMuonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
      float fTransverseMomentum = kTruthTransverseMomentum(nu);
      if (fTransverseMomentum < TwoDArrayTransverseMomentum[0]) { fTransverseMomentum = (TwoDArrayTransverseMomentum[0] + TwoDArrayTransverseMomentum[1]) / 2.; }
      else if (fTransverseMomentum > TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum]) { fTransverseMomentum = (TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum] + TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum - 1]) / 2.; }

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kTruthMuonCosTheta(nu), TwoDArrayNBinsMuonCosTheta);
      int SerialTransverseMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
										  TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices,
										  MuonCosThetaTwoDIndex,
										  fTransverseMomentum
										  );
      return SerialTransverseMomentumInMuonCosThetaIndex; 
    });


  const Var kRecoTruthTransverseMomentumInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthTransverseMomentumInMuonCosTheta(&slc->truth); });

  // Delta alpha transverse
  const Var kDeltaAlphaTInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
      float fDeltaAlphaT = kDeltaAlphaT(slc);

      if (fDeltaAlphaT < TwoDArrayDeltaAlphaT[0]) { fDeltaAlphaT = (TwoDArrayDeltaAlphaT[0] + TwoDArrayDeltaAlphaT[1]) / 2.; }
      else if (fDeltaAlphaT > TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT]) { fDeltaAlphaT = (TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT] + TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT - 1]) / 2.; }

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kMuonCosTheta(slc), TwoDArrayNBinsMuonCosTheta);
      int SerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
									   TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,
									   MuonCosThetaTwoDIndex,
									   fDeltaAlphaT
									   );
      return SerialDeltaAlphaTInMuonCosThetaIndex; });
  const TruthVar kTruthDeltaAlphaTInMuonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
      float fDeltaAlphaT = kTruthDeltaAlphaT(nu);

      if (fDeltaAlphaT < TwoDArrayDeltaAlphaT[0]) { fDeltaAlphaT = (TwoDArrayDeltaAlphaT[0] + TwoDArrayDeltaAlphaT[1]) / 2.; }
      else if (fDeltaAlphaT > TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT]) { fDeltaAlphaT = (TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT] + TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT - 1]) / 2.; }

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kTruthMuonCosTheta(nu), TwoDArrayNBinsMuonCosTheta);
      int SerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
									   TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,
									   MuonCosThetaTwoDIndex,
									   fDeltaAlphaT
									   );
      return SerialDeltaAlphaTInMuonCosThetaIndex;
    });

  const Var kRecoTruthDeltaAlphaTInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthDeltaAlphaTInMuonCosTheta(&slc->truth); });

  // Opening angle between protons
  const Var kCosOpeningAngleProtonsInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
      float fCosOpeningAngleProtons = kCosOpeningAngleProtons(slc);

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kMuonCosTheta(slc), TwoDArrayNBinsMuonCosTheta);
      int SerialCosOpeningAngleProtonsInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
										      TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices,
										      MuonCosThetaTwoDIndex,
										      fCosOpeningAngleProtons
										      );
      return SerialCosOpeningAngleProtonsInMuonCosThetaIndex; 
    });

  const TruthVar kTruthCosOpeningAngleProtonsInMuonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
      float fCosOpeningAngleProtons = kTruthCosOpeningAngleProtons(nu);

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kTruthMuonCosTheta(nu), TwoDArrayNBinsMuonCosTheta);
      int SerialCosOpeningAngleProtonsInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
										      TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices,
										      MuonCosThetaTwoDIndex,
										      fCosOpeningAngleProtons
										      );
      return SerialCosOpeningAngleProtonsInMuonCosThetaIndex; 
    });
 
  const Var kRecoTruthCosOpeningAngleProtonsInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthCosOpeningAngleProtonsInMuonCosTheta(&slc->truth); });

  // Opening angle betwen muon and total proton
  const Var kCosOpeningAngleMuonTotalProtonInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
      float fCosOpeningAngleMuonTotalProton = kCosOpeningAngleMuonTotalProton(slc);

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kMuonCosTheta(slc), TwoDArrayNBinsMuonCosTheta);
      int SerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
											      TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices,
											      MuonCosThetaTwoDIndex,
											      fCosOpeningAngleMuonTotalProton
											      );
      return SerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex; 
    });


  const TruthVar kTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
      float fCosOpeningAngleMuonTotalProton = kTruthCosOpeningAngleMuonTotalProton(nu);

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kTruthMuonCosTheta(nu), TwoDArrayNBinsMuonCosTheta);
      int SerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
											      TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices,
											      MuonCosThetaTwoDIndex,
											      fCosOpeningAngleMuonTotalProton
											      );
      return SerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex; });


  const Var kRecoTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta(&slc->truth); });

  // Missing momentum
  const Var kMissingMomentumInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
      float fMissingMomentum = kMissingMomentum(slc);

      if (fMissingMomentum < TwoDArrayMissingMomentum[0]) { fMissingMomentum = (TwoDArrayMissingMomentum[0] + TwoDArrayMissingMomentum[1]) / 2.; }
      else if (fMissingMomentum > TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum]) { fMissingMomentum = (TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum] + TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum - 1]) / 2.; }

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kMuonCosTheta(slc), TwoDArrayNBinsMuonCosTheta);
      int SerialMissingMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
									       TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices,
									       MuonCosThetaTwoDIndex,
									       fMissingMomentum
									       );
      return SerialMissingMomentumInMuonCosThetaIndex; });

  const TruthVar kTruthMissingMomentumInMuonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
      float fMissingMomentum = kTruthMissingMomentum(nu);

      if (fMissingMomentum < TwoDArrayMissingMomentum[0]) { fMissingMomentum = (TwoDArrayMissingMomentum[0] + TwoDArrayMissingMomentum[1]) / 2.; }
      else if (fMissingMomentum > TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum]) { fMissingMomentum = (TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum] + TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum - 1]) / 2.; }

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kTruthMuonCosTheta(nu), TwoDArrayNBinsMuonCosTheta);
      int SerialMissingMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
									       TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices,
									       MuonCosThetaTwoDIndex,
									       fMissingMomentum
									       );
      return SerialMissingMomentumInMuonCosThetaIndex; 
    });


  const Var kRecoTruthMissingMomentumInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthMissingMomentumInMuonCosTheta(&slc->truth); });

  // Alpha three dimenional
  const Var kAlphaThreeDInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
      float fAlphaThreeD = kAlphaThreeD(slc);

      if (fAlphaThreeD < TwoDArrayAlphaThreeD[0]) { fAlphaThreeD = (TwoDArrayAlphaThreeD[0] + TwoDArrayAlphaThreeD[1]) / 2.; }
      else if (fAlphaThreeD > TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD]) { fAlphaThreeD = (TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD] + TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD - 1]) / 2.; }

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kMuonCosTheta(slc), TwoDArrayNBinsMuonCosTheta);
      int SerialAlphaThreeDInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
									   TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices,
									   MuonCosThetaTwoDIndex,
									   fAlphaThreeD
									   );
      return SerialAlphaThreeDInMuonCosThetaIndex; 
    });

  const TruthVar kTruthAlphaThreeDInMuonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
      float fAlphaThreeD = kTruthAlphaThreeD(nu);

      if (fAlphaThreeD < TwoDArrayAlphaThreeD[0]) { fAlphaThreeD = (TwoDArrayAlphaThreeD[0] + TwoDArrayAlphaThreeD[1]) / 2.; }
      else if (fAlphaThreeD > TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD]) { fAlphaThreeD = (TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD] + TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD - 1]) / 2.; }

      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kTruthMuonCosTheta(nu), TwoDArrayNBinsMuonCosTheta);
      int SerialAlphaThreeDInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
									   TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices,
									   MuonCosThetaTwoDIndex,
									   fAlphaThreeD
									   );
      return SerialAlphaThreeDInMuonCosThetaIndex; 
    });

  const Var kRecoTruthAlphaThreeDInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double { return kTruthAlphaThreeDInMuonCosTheta(&slc->truth); });

  // Opening angle betwen momentum transfer vector and total proton
  const Var kCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
      float fCosOpeningAngleMomentumTransferTotalProton = kCosOpeningAngleMomentumTransferTotalProton(slc);
      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kMuonCosTheta(slc), TwoDArrayNBinsMuonCosTheta);
      int SerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
													  TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices,
													  MuonCosThetaTwoDIndex,
													  fCosOpeningAngleMomentumTransferTotalProton
													  );
      return SerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex;
    });


  const TruthVar kTruthCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
      float fCosOpeningAngleMomentumTransferTotalProton = kTruthCosOpeningAngleMomentumTransferTotalProton(nu);
      int MuonCosThetaTwoDIndex = tools.ReturnIndex(kTruthMuonCosTheta(nu), TwoDArrayNBinsMuonCosTheta);
      int SerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
													  TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices,
													  MuonCosThetaTwoDIndex,
													  fCosOpeningAngleMomentumTransferTotalProton
													  );
      return SerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex; 
    });
  
  const Var kRecoTruthCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
      return kTruthCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta(&slc->truth); 
    });

  //////////////
  // Truth Cuts
  //////////////

  const TruthCut kTruthIsSignal([](const caf::SRTrueInteractionProxy *nu) {
      return(
	     bIsInFV(&nu->position) &&                                                                               // check position is in fiducial volume
	     nu->iscc &&                                                                                             // check it is charged current interaction
	     nu->pdg == 14 &&                                                                                        // check neutrino is muon neutrino
	     iCountMultParticle(nu, 13, std::get<0>(PDGToThreshold.at(13)), std::get<1>(PDGToThreshold.at(13))) == 1 //&&       // check for one muon
	     // iCountMultParticle(nu, 2212, std::get<0>(PDGToThreshold.at(2212)), std::get<1>(PDGToThreshold.at(2212))) == 2 && // check for two protons
	     // iCountMultParticle(nu, 211, std::get<0>(PDGToThreshold.at(211)), std::get<1>(PDGToThreshold.at(211))) == 0 &&    // no positively charged pions
	     // iCountMultParticle(nu, -211, std::get<0>(PDGToThreshold.at(-211)), std::get<1>(PDGToThreshold.at(-211))) == 0 && // no negatively charged pions
	     // iCountMultParticle(nu, 111, std::get<0>(PDGToThreshold.at(111)), std::get<1>(PDGToThreshold.at(111))) == 0       // no neutral pions
	     ); 
    });
  
  const TruthCut kTruthNoSignal([](const caf::SRTrueInteractionProxy *nu) { return !kTruthIsSignal(nu); });

  // Truth cuts for other topologies

  // CCnumu 1pi=/-, 0pi0,  (no idea what Ngt means?)
  const TruthCut kCCNgt0p1pi([](const caf::SRTrueInteractionProxy *nu) {
      int nChargedPions = iCountMultParticle(nu, 211, std::get<0>(PDGToThreshold.at(211)), std::get<1>(PDGToThreshold.at(211))) + iCountMultParticle(nu, -211, std::get<0>(PDGToThreshold.at(-211)), std::get<1>(PDGToThreshold.at(-211)));
      return (
	      bIsInFV(&nu->position) &&
	      nu->iscc &&
	      nu->pdg == 14 &&
	      iCountMultParticle(nu, 13, std::get<0>(PDGToThreshold.at(13)), std::get<1>(PDGToThreshold.at(13))) == 1 && // check for one muon
	      nChargedPions == 1 &&                                                                                      // one charged pion
	      iCountMultParticle(nu, 111, std::get<0>(PDGToThreshold.at(111)), std::get<1>(PDGToThreshold.at(111))) == 0 // no neutral pions
	      ); });

  const TruthCut kCCNg2p0pi([](const caf::SRTrueInteractionProxy *nu) { 
      return (
	      bIsInFV(&nu->position) &&
	      nu->iscc &&
	      nu->pdg == 14 &&
	      iCountMultParticle(nu, 13, std::get<0>(PDGToThreshold.at(13)), std::get<1>(PDGToThreshold.at(13))) == 1 &&       // check for one muon
	      iCountMultParticle(nu, 2212, std::get<0>(PDGToThreshold.at(2212)), std::get<1>(PDGToThreshold.at(2212))) > 2 &&  // check for more than two protons
	      iCountMultParticle(nu, 211, std::get<0>(PDGToThreshold.at(211)), std::get<1>(PDGToThreshold.at(211))) == 0 &&    // no positively charged pion
	      iCountMultParticle(nu, -211, std::get<0>(PDGToThreshold.at(-211)), std::get<1>(PDGToThreshold.at(-211))) == 0 && // no negatively charged pion
	      iCountMultParticle(nu, 111, std::get<0>(PDGToThreshold.at(111)), std::get<1>(PDGToThreshold.at(111))) == 0       // no neutral pions
	      ); 
    });

  const TruthCut kCC1p0pi([](const caf::SRTrueInteractionProxy *nu) { 
      return (
	      bIsInFV(&nu->position) &&
	      nu->iscc &&
	      nu->pdg == 14 &&
	      iCountMultParticle(nu, 13, std::get<0>(PDGToThreshold.at(13)), std::get<1>(PDGToThreshold.at(13))) == 1 &&       // check for one muon
	      iCountMultParticle(nu, 2212, std::get<0>(PDGToThreshold.at(2212)), std::get<1>(PDGToThreshold.at(2212))) == 1 && // check for one proton
	      iCountMultParticle(nu, 211, std::get<0>(PDGToThreshold.at(211)), std::get<1>(PDGToThreshold.at(211))) == 0 &&    // no positively charged pion
	      iCountMultParticle(nu, -211, std::get<0>(PDGToThreshold.at(-211)), std::get<1>(PDGToThreshold.at(-211))) == 0 && // no negatively charged pion
	      iCountMultParticle(nu, 111, std::get<0>(PDGToThreshold.at(111)), std::get<1>(PDGToThreshold.at(111))) == 0       // no neutral pions
	      ); });

  const TruthCut kCC0p0pi([](const caf::SRTrueInteractionProxy *nu) { 
      return (
	      bIsInFV(&nu->position) &&
	      nu->iscc &&
	      nu->pdg == 14 &&
	      iCountMultParticle(nu, 13, std::get<0>(PDGToThreshold.at(13)), std::get<1>(PDGToThreshold.at(13))) == 1 &&       // check for one muon
	      iCountMultParticle(nu, 2212, std::get<0>(PDGToThreshold.at(2212)), std::get<1>(PDGToThreshold.at(2212))) == 0 && // no protons
	      iCountMultParticle(nu, 211, std::get<0>(PDGToThreshold.at(211)), std::get<1>(PDGToThreshold.at(211))) == 0 &&    // no positively charged pion
	      iCountMultParticle(nu, -211, std::get<0>(PDGToThreshold.at(-211)), std::get<1>(PDGToThreshold.at(-211))) == 0 && // no negatively charged pion
	      iCountMultParticle(nu, 111, std::get<0>(PDGToThreshold.at(111)), std::get<1>(PDGToThreshold.at(111))) == 0       // no neutral pions
	      ); 
    });

  const TruthCut kCCNgt0pNg1pi([](const caf::SRTrueInteractionProxy *nu) {
      int nChargedPions = iCountMultParticle(nu, 211, std::get<0>(PDGToThreshold.at(211)), std::get<1>(PDGToThreshold.at(211))) + iCountMultParticle(nu, -211, std::get<0>(PDGToThreshold.at(-211)), std::get<1>(PDGToThreshold.at(-211)));
      return (
	      bIsInFV(&nu->position) &&
	      nu->iscc &&
	      nu->pdg == 14 &&
	      iCountMultParticle(nu, 13, std::get<0>(PDGToThreshold.at(13)), std::get<1>(PDGToThreshold.at(13))) == 1 && // check for one muon
	      nChargedPions > 1 &&                                                                                       // more than one charged pion
	      iCountMultParticle(nu, 111, std::get<0>(PDGToThreshold.at(111)), std::get<1>(PDGToThreshold.at(111))) == 0 // no neutral pions
	      ); 
    });

  const TruthCut kOtherTopology([](const caf::SRTrueInteractionProxy *nu) { 
      return !(
	       kTruthIsSignal(nu) ||
	       kCC1p0pi(nu) ||
	       kCCNg2p0pi(nu) ||
	       kCCNgt0p1pi(nu) ||
	       kCC0p0pi(nu)
	       // kCCNgt0pNg1pi
	       ); 
    });

  //////////////
  // Cuts
  //////////////

  // Check fmatch is in beam
  const Cut kIsInBeam([](const caf::SRSliceProxy *slc) { return ((slc->fmatch.time > 0.) && (slc->fmatch.time < 1.800)); });

  // Check for cosmics (not const for data script)
  Cut kCosmicCut([](const caf::SRSliceProxy *slc) { 
      return (
	      slc->nu_score > 0.4 &&     // check how neutrino like slice is
	      slc->fmatch.score < 7.0 && // check flash match score
	      slc->fmatch.time > 0. &&   // check flash is in beam
	      slc->fmatch.time < 1.8); 
    });
  
  const Cut kIsCosmic([](const caf::SRSliceProxy *slc) { return (slc->truth.genie_mode == -1); });
  const Cut kIsNotCosmic([](const caf::SRSliceProxy *slc) { return (slc->truth.genie_mode != -1); });

  const Cut kHasMuon([](const caf::SRSliceProxy *slc) {
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 13) {
	  return true;
	}
      }
      return false;
    });
  
  const Cut kHasProton([](const caf::SRSliceProxy *slc) {
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 2212) {
	  return true;
	}
      }
      return false; 
    });

  const Cut kHasSecondProton([](const caf::SRSliceProxy *slc) {
      bool firstProton = false;
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 2212) {
	  if (!firstProton) firstProton = true;
	  else if (firstProton) return true;
	}
      }
      return false; 
    });

  const Cut kHasPion([](const caf::SRSliceProxy *slc) {
      for (auto const& pfp : slc -> reco.pfp) {
	if (pfp.trk.truth.p.pdg == 211 || pfp.trk.truth.p.pdg == -211) {
	  return true;
	}
      }
      return false; 
    });

  // Check reconstructed event is signal
  const Cut kRecoIsSignal([](const caf::SRSliceProxy *slc) {
      std::vector<int> TaggedIDs;

      // Reject cosmic events
      if (!kCosmicCut(slc)) return false; 

      // Check neutrino vertex is in fiducial volume
      if (!bIsInFV(&slc->vertex)) return false;

      // Check there is one muon in signal
      auto [OneMuon, MuonID] = bOneMuon(slc);
      if (!OneMuon) return false;
      TaggedIDs.push_back(MuonID);

      // Check there are two protons in signal
      //auto [TwoProtons, ProtonIDs] = bTwoProtons(slc, MuonID);
      //if (!TwoProtons) return false;
      //TaggedIDs.insert(TaggedIDs.end(), ProtonIDs.begin(), ProtonIDs.end());

      // Check there are no charged pions
      //if (!bNoChargedPions(slc, TaggedIDs)) return false;

      // // Check there are no shower-like objects (neutral pions)
      //if (!bNoShowers(slc, TaggedIDs)) return false;

      // Signal definition satisifed
      return true; 
    });

  const Cut kRecoIsTrueReco([](const caf::SRSliceProxy *slc) { return (kRecoIsSignal(slc) && kTruthIsSignal(&slc->truth)); });

  const Cut kRecoIsBackground([](const caf::SRSliceProxy *slc) { return (kRecoIsSignal(slc) && kTruthNoSignal(&slc->truth)); });

  const Cut kNoInvalidVariables([](const caf::SRSliceProxy *slc) {
      if (std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z)) return false;
      // for (auto const& pfp : slc -> reco.pfp) {
      //     if (std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.start.y) || std::isnan(pfp.trk.start.z)) return false;
      //     if (std::isnan(pfp.trk.end.x)   || std::isnan(pfp.trk.end.y)   || std::isnan(pfp.trk.end.z)) return false;
      //     if (std::isnan(pfp.trk.len)) return false;
      //     if (std::isnan(pfp.trk.mcsP.fwdP_muon) || std::isnan(pfp.trk.rangeP.p_muon)) return false;

      //     for (int i = 0; i < 3; i++) {
      //         if (
      //             std::isnan(pfp.trk.chi2pid[i].chi2_muon) ||
      //             pfp.trk.chi2pid[i].chi2_muon == 0. ||
      //             std::isnan(pfp.trk.chi2pid[i].chi2_proton) ||
      //             pfp.trk.chi2pid[i].chi2_proton == 0.
      //         ) return false;
      //     }
      // }
      // TODO: figure out a better way to implement this, we don't really care about ALL
      //       the daughter particles, but only about those that we tag successfully and 
      //       we will use to compute our reconstructed variables
      return true; });

  const Cut kTrueSignal([](const caf::SRSliceProxy *slc) { return (kNoInvalidVariables(slc) && kTruthIsSignal(&slc->truth)); });

  // These cuts cumulatively reconstruct kRecoIsSignal, to compute cut
  // efficiencies and purity. The cuts are applied in the following order:
  //     1. Cosmic cut
  //     2. Neutrino vertex in fiducial volume cut
  //     3. One muon cut
  //     4. Two protons cut
  //     5. No charged pions cut
  //     6. No neutral pions cut

  const Cut kFirstCut([](const caf::SRSliceProxy *slc) { return kCosmicCut(slc); });
  const Cut kFirstCutTrue([](const caf::SRSliceProxy *slc) { return (kFirstCut(slc) && kTruthIsSignal(&slc->truth)); });

  const Cut kSecondCut([](const caf::SRSliceProxy *slc) { return (kFirstCut(slc) && bIsInFV(&slc->vertex)); });
  const Cut kSecondCutTrue([](const caf::SRSliceProxy *slc) { return (kSecondCut(slc) && kTruthIsSignal(&slc->truth)); });

  const Cut kThirdCut([](const caf::SRSliceProxy *slc) {
      auto [OneMuon, MuonID] = bOneMuon(slc);
      return (OneMuon && kSecondCut(slc)); });

  const Cut kThirdCutTrue([](const caf::SRSliceProxy *slc) { return (kThirdCut(slc) && kTruthIsSignal(&slc->truth)); });

  const Cut kFourthCut([](const caf::SRSliceProxy *slc) {
      auto [OneMuon, MuonID] = bOneMuon(slc);
      if (!OneMuon) return false;
      auto [TwoProtons, ProtonIDs] = bTwoProtons(slc, MuonID);
          
      return (TwoProtons && kSecondCut(slc)); 
    });

  const Cut kFourthCutTrue([](const caf::SRSliceProxy *slc) { return (kFourthCut(slc) && kTruthIsSignal(&slc->truth)); });

  const Cut kFifthCut([](const caf::SRSliceProxy *slc) {
      std::vector<int> TaggedIDs;

      auto [OneMuon, MuonID] = bOneMuon(slc);
      if (!OneMuon) return false;
      TaggedIDs.push_back(MuonID);

      auto [TwoProtons, ProtonIDs] = bTwoProtons(slc, MuonID);
      if (!TwoProtons) return false;
      TaggedIDs.insert(TaggedIDs.end(), ProtonIDs.begin(), ProtonIDs.end());

      return (bNoChargedPions(slc, TaggedIDs) && kSecondCut(slc)); 
    });

  const Cut kFifthCutTrue([](const caf::SRSliceProxy *slc) { return (kFifthCut(slc) && kTruthIsSignal(&slc->truth)); });

  const Cut kSixthCut([](const caf::SRSliceProxy *slc) {
      std::vector<int> TaggedIDs;

      auto [OneMuon, MuonID] = bOneMuon(slc);
      if (!OneMuon) return false;
      TaggedIDs.push_back(MuonID);

      auto [TwoProtons, ProtonIDs] = bTwoProtons(slc, MuonID);
      if (!TwoProtons) return false;
      TaggedIDs.insert(TaggedIDs.end(), ProtonIDs.begin(), ProtonIDs.end());

      return (bNoChargedPions(slc, TaggedIDs) && bNoShowers(slc, TaggedIDs) && kSecondCut(slc)); });
			
  const Cut kSixthCutTrue([](const caf::SRSliceProxy *slc) { return (kSixthCut(slc) && kTruthIsSignal(&slc->truth)); });

  ///////////
  // SpillVar
  ///////////

  const SpillCut kSpillPrintFile([](const caf::StandardRecordProxy *sr) {
      std::cout << (std::string)sr->hdr.sourceName << "  " << sr->hdr.sourceIndex << "  " << sr->hdr.fno << std::endl; 
      return true; });

  const SpillVar kSpillData([](const caf::StandardRecordProxy *sr) {
      fstream file;
      std::string FileName = "/exp/sbnd/data/users/" + UserName + "/CAFAnaOutput/EventData.csv";
      file.open(FileName, fstream::out | fstream::app);
      for (auto const& slc : sr->slc) {
	if (kRecoIsSignal(&slc)) {
	  file << sr->hdr.fno << ",";
	  file << sr->hdr.run << ",";
	  file << sr->hdr.subrun << ",";
	  file << sr->hdr.evt << ",";
	  file << sr->hdr.subevt << std::endl;       
	}
      }
      return 0.5; 
    });

  /////////////////////////////
  // Vector with varialbes/bins
  /////////////////////////////

  static const std::vector<std::tuple<Var, Var, TruthVar>> VarsOld = {
    {kEventCount, kEventCount, kTrueEventCount},
    {kVertexX, kRecoTruthVertexX, kTruthVertexX},
    //{kVertexY, kRecoTruthVertexY, kTruthVertexY},
    //{kVertexZ, kRecoTruthVertexZ, kTruthVertexZ},
    //{kMuonCosTheta, kRecoTruthMuonCosTheta, kTruthMuonCosTheta},
    //{kLeadingProtonCosTheta, kRecoTruthLeadingProtonCosTheta, kTruthLeadingProtonCosTheta},
    //{kRecoilProtonCosTheta, kRecoTruthRecoilProtonCosTheta, kTruthRecoilProtonCosTheta},
    //{kCosOpeningAngleProtons, kRecoTruthCosOpeningAngleProtons, kTruthCosOpeningAngleProtons},
    //{kCosOpeningAngleMuonTotalProton, kRecoTruthCosOpeningAngleMuonTotalProton, kTruthCosOpeningAngleMuonTotalProton},
    //{kDeltaAlphaT, kRecoTruthDeltaAlphaT, kTruthDeltaAlphaT},
    //{kTransverseMomentum, kRecoTruthTransverseMomentum, kTruthTransverseMomentum},
    {kMuonMomentum, kRecoTruthMuonMomentum, kTruthMuonMomentum} //,
    //{kLeadingProtonMomentum, kRecoTruthLeadingProtonMomentum, kTruthLeadingProtonMomentum},
    //{kRecoilProtonMomentum, kRecoTruthRecoilProtonMomentum, kTruthRecoilProtonMomentum},
    //{kCosOpeningAngleMomentumTransferTotalProton, kRecoTruthCosOpeningAngleMomentumTransferTotalProton, kTruthCosOpeningAngleMomentumTransferTotalProton},
    //{kAlphaThreeD, kRecoTruthAlphaThreeD, kTruthAlphaThreeD},
    //{kMissingMomentum, kRecoTruthMissingMomentum, kTruthMissingMomentum},
    //{kInvariantMass, kRecoTruthInvariantMass, kTruthInvariantMass},
    //{kCosOpeningAngleLProtonMuon, kRecoTruthCosOpeningAngleLProtonMuon, kTruthCosOpeningAngleLProtonMuon},
    //{kCosOpeningAngleRProtonMuon, kRecoTruthCosOpeningAngleRProtonMuon, kTruthCosOpeningAngleRProtonMuon},
    //{kTransverseMomentumInMuonCosTheta, kRecoTruthTransverseMomentumInMuonCosTheta, kTruthTransverseMomentumInMuonCosTheta},
    //{kDeltaAlphaTInMuonCosTheta, kRecoTruthDeltaAlphaTInMuonCosTheta, kTruthDeltaAlphaTInMuonCosTheta},
    //{kCosOpeningAngleProtonsInMuonCosTheta, kRecoTruthCosOpeningAngleProtonsInMuonCosTheta, kTruthCosOpeningAngleProtonsInMuonCosTheta},
    //{kCosOpeningAngleMuonTotalProtonInMuonCosTheta, kRecoTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta, kTruthCosOpeningAngleMuonTotalProtonInMuonCosTheta},
    //{kMissingMomentumInMuonCosTheta, kRecoTruthMissingMomentumInMuonCosTheta, kTruthMissingMomentumInMuonCosTheta},
    //{kAlphaThreeDInMuonCosTheta, kRecoTruthAlphaThreeDInMuonCosTheta, kTruthAlphaThreeDInMuonCosTheta},
    //{kCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta, kRecoTruthCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta, kTruthCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta}
  };

  static const std::vector<Binning> VarBinsOld = {
    bEventCount,
    bVertexX,
    bVertexY,
    bVertexZ,
    bAngleBins,
    bAngleBins,
    bAngleBins,
    bAngleBins,
    bAngleBins,
    bDeltaAlphaBins,
    bTransverseMomentumBins,
    bMuonMomentumBins,
    bLeadingProtonMomentumBins,
    bRecoilProtonMomentumBins,
    bAngleBins,
    bAlphaThreeDBins,
    bMissingMomentumBins,
    bInvariantMassBins,
    bCosOpeningAngleLProtonMuonBins,
    bCosOpeningAngleRProtonMuonBins,
    bTransverseMomentumInMuonCosTheta,
    bDeltaAlphaTInMuonCosTheta,
    bCosOpeningAngleProtonsInMuonCosTheta,
    bCosOpeningAngleMuonTotalProtonInMuonCosTheta,
    bMissingMomentumInMuonCosTheta,
    bAlphaThreeDInMuonCosTheta,
    bCosOpeningAngleMomentumTransferTotalProtonInMuonCosTheta
  };

  ////////////////////////////
  // Weights for our fake data
  ////////////////////////////

  const Var kDoubleMECWeight([](const caf::SRSliceProxy *slc) -> double {
      if (slc->truth.genie_mode == 10) return 2;
      return 1; });

  const TruthVar kTrueDoubleMECWeight([](const caf::SRTrueInteractionProxy *nu) -> double {
      if (nu->genie_mode == 10) return 2;
      return 1; });

  const Var kDoubleQEWeight([](const caf::SRSliceProxy *slc) -> double {
      if (slc->truth.genie_mode == 0) return 2;
      return 1; });

  const TruthVar kTrueDoubleQEWeight([](const caf::SRTrueInteractionProxy *nu) -> double {
      if (nu->genie_mode == 0) return 2;
      return 1; });

  const Var kCombinedWeight([](const caf::SRSliceProxy *slc) -> double {
      if (slc->truth.genie_mode == 0) return 0.5;
      if (slc->truth.genie_mode == 10) return 1.5;
      return 1; });

  const TruthVar kTrueCombinedWeight([](const caf::SRTrueInteractionProxy *nu) -> double {
      if (nu->genie_mode == 0) return 0.5;
      if (nu->genie_mode == 10) return 1.5;
      return 1; });

  static const std::vector<std::tuple<Var, TruthVar>> FakeWeights = {
    {kDoubleMECWeight, kTrueDoubleMECWeight},
    {kDoubleQEWeight, kTrueDoubleQEWeight},
    {kCombinedWeight, kTrueCombinedWeight}
  };
} // namespace ana

#endif
