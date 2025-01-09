#define FlatTreeAnalyzer_cxx
#include "FlatTreeAnalyzer.h"
#include "../../Utils/TwoPTools.h"
#include "../../Utils/Tools.h"
#include "../../Utils/Constants.h"

#include <TH1D.h>
#include <TRandom.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>

using namespace std;
using namespace Constants;

// Function to divide by the bin width and to get xsecs
void Reweight(TH1D* h);

//----------------------------------------//

void FlatTreeAnalyzer::Loop() {

    //----------------------------------------//	

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
    double A = 40.; // so that we can have xsecs per nucleus

    int NInte = 6; // Interaction processes: All, QE, MEC, RES, DIS, COH
    std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};

    //----------------------------------------//	

    Tools tools;

    //----------------------------------------//	

    // Output file
    TString Directory = "/pnfs/sbnd/persistent/users/" + (TString)UserName + "/HighSamples/";
    TString FileNameAndPath = Directory+"FlatTree/FlatTreeAnalyzerOutput_"+fOutputFile+".root";
    TFile* file = new TFile(FileNameAndPath,"recreate");

    std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
    std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
    
    //----------------------------------------//

    // Plot declaration
    
    // Single differential
    TH1D* TrueVertexXPlot[NInte];
    TH1D* TrueVertexYPlot[NInte];
    TH1D* TrueVertexZPlot[NInte];
    TH1D* TrueMuonCosThetaPlot[NInte];
    TH1D* TrueLeadingProtonCosThetaPlot[NInte];
    TH1D* TrueRecoilProtonCosThetaPlot[NInte];
    TH1D* TrueLeadingProtonMomentumPlot[NInte];
    TH1D* TrueRecoilProtonMomentumPlot[NInte];
    TH1D* TrueMuonMomentumPlot[NInte];
    TH1D* TrueCosOpeningAngleProtonsPlot[NInte];
    TH1D* TrueCosOpeningAngleMuonTotalProtonPlot[NInte];
    TH1D* TrueTransverseMomentumPlot[NInte];
    TH1D* TrueDeltaAlphaTPlot[NInte];
    TH1D* TrueInvariantMassPlot[NInte];
    TH1D* TrueCosOpeningAngleLProtonMuonPlot[NInte]; 
    TH1D* TrueCosOpeningAngleRProtonMuonPlot[NInte];   
	
    TH1D* TrueNoFSILeadingProtonCosThetaPlot[NInte];
    TH1D* TrueNoFSIRecoilProtonCosThetaPlot[NInte];
    TH1D* TrueNoFSILeadingProtonMomentumPlot[NInte];
    TH1D* TrueNoFSIRecoilProtonMomentumPlot[NInte];
    TH1D* TrueNoFSIMuonMomentumPlot[NInte];
    TH1D* TrueNoFSICosOpeningAngleProtonsPlot[NInte];
    TH1D* TrueNoFSICosOpeningAngleMuonTotalProtonPlot[NInte];
    TH1D* TrueNoFSITransverseMomentumPlot[NInte];
    TH1D* TrueNoFSIDeltaAlphaTPlot[NInte];
    TH1D* TrueNoFSIInvariantMassPlot[NInte];
    TH1D* TrueNoFSICosOpeningAngleLProtonMuonPlot[NInte];
    TH1D* TrueNoFSICosOpeningAngleRProtonMuonPlot[NInte];

    // Single differential - GKI
    TH1D* TrueCosOpeningAngleMomentumTransferTotalProtonPlot[NInte];
    TH1D* TrueMissingMomentumPlot[NInte];
    TH1D* TrueAlphaThreeDPlot[NInte];
    
    TH1D* TrueNoFSICosOpeningAngleMomentumTransferTotalProtonPlot[NInte];
    TH1D* TrueNoFSIMissingMomentumPlot[NInte];
    TH1D* TrueNoFSIAlphaThreeDPlot[NInte];

    // Double differential
    TH1D* TrueSerialTransverseMomentum_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialDeltaAlphaT_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot[NInte];

    TH1D* TrueSerialNoFSITransverseMomentum_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialNoFSIDeltaAlphaT_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialNoFSICosOpeningAngleProtons_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialNoFSICosOpeningAngleMuonTotalProton_InMuonCosThetaPlot[NInte];

    // Double differential - GKI
    TH1D* TrueSerialMissingMomentum_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialAlphaThreeD_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialCosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot[NInte];

    TH1D* TrueSerialNoFSIMissingMomentum_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialNoFSIAlphaThreeD_InMuonCosThetaPlot[NInte];
    TH1D* TrueSerialNoFSICosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot[NInte];

    // Loop over the interaction processes

    for (int inte = 0; inte < NInte; inte++) {
        //--------------------------------------------------//

        // Final state
        TrueVertexXPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexXPlot",";#vec{v}_{x} [cm]",NBinsVertexX,ArrayNBinsVertexX.data());
        TrueVertexYPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexYPlot",";#vec{v}_{y} [cm]",NBinsVertexY,ArrayNBinsVertexY.data());
        TrueVertexZPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueVertexZPlot",";#vec{v}_{z} [cm]",NBinsVertexZ,ArrayNBinsVertexZ.data());
        TrueMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonCosThetaPlot",";cos(#theta_{#vec{p}_{#mu}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueLeadingProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueLeadingProtonCosThetaPlot",";cos(#theta_{#vec{p}_{L}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueRecoilProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueRecoilProtonCosThetaPlot",";cos(#theta_{#vec{p}_{R}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueLeadingProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueLeadingProtonMomentumPlot",";|#vec{p}_{L}| [GeV/c]",NBinsLeadingProtonMomentum,ArrayNBinsLeadingProtonMomentum.data());
        TrueRecoilProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueRecoilProtonMomentumPlot",";|#vec{p}_{R}| [GeV/c]",NBinsRecoilProtonMomentum,ArrayNBinsRecoilProtonMomentum.data());
        TrueMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMuonMomentumPlot",";|#vec{p}_{#mu}| [GeV/c]",NBinsMuonMomentum,ArrayNBinsMuonMomentum.data());
        TrueCosOpeningAngleProtonsPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosOpeningAngleProtonsPlot",";cos(#theta_{#vec{p}_{L},#vec{p}_{R}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueCosOpeningAngleMuonTotalProtonPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosOpeningAngleMuonTotalProtonPlot",";cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueTransverseMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueTransverseMomentumPlot",";#delta P_{T} [GeV/c]",NBinsTransverseMomentum,ArrayNBinsTransverseMomentum.data());
        TrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",";#delta #alpha_{T} [deg]",NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT.data());
        
	// Final state - GKI
        TrueCosOpeningAngleMomentumTransferTotalProtonPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosOpeningAngleMomentumTransferTotalProtonPlot",";cos(#theta_{#vec{q},#vec{p}_{sum}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueMissingMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueMissingMomentumPlot",";p_{n} [GeV/c]",NBinsMissingMomentum,ArrayNBinsMissingMomentum.data());
        TrueAlphaThreeDPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueAlphaThreeDPlot",";#alpha_{3D} [deg]",NBinsAlphaThreeD,ArrayNBinsAlphaThreeD.data());
        
	//Final state - additional var
	TrueInvariantMassPlot[inte] = new TH1D(InteractionLabels[inte] + "TrueInvariantMassPlot", ";W [GeV]", NBinsInvariantMass, ArrayNBinsInvariantMass.data());
        TrueCosOpeningAngleLProtonMuonPlot[inte] = new TH1D(InteractionLabels[inte] + "TrueCosOpeningAngleLProtonMuonPlot", ";cos(#theta_{#vec{p}_{L},#vec{p}_{#mu}})", NBinsCosAngleLPMu, ArrayNBinsCosAngleLPMu.data());
        TrueCosOpeningAngleRProtonMuonPlot[inte] = new TH1D(InteractionLabels[inte] + "TrueCosOpeningAngleRProtonMuonPlot", ";cos(#theta_{#vec{p}_{R},#vec{p}_{#mu}})", NBinsCosAngleLPMu, ArrayNBinsCosAngleLPMu.data());


	// Before final state interactions
        TrueNoFSILeadingProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSILeadingProtonCosThetaPlot",";cos(#theta_{#vec{p}_{L}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueNoFSIRecoilProtonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIRecoilProtonCosThetaPlot",";cos(#theta_{#vec{p}_{R}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueNoFSILeadingProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSILeadingProtonMomentumPlot",";|#vec{p}_{L}| [GeV/c]",NBinsLeadingProtonMomentum,ArrayNBinsLeadingProtonMomentum.data());
        TrueNoFSIRecoilProtonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIRecoilProtonMomentumPlot",";|#vec{p}_{R}| [GeV/c]",NBinsRecoilProtonMomentum,ArrayNBinsRecoilProtonMomentum.data());
        TrueNoFSIMuonMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIMuonMomentumPlot",";|#vec{p}_{#mu}| [GeV/c]",NBinsMuonMomentum,ArrayNBinsMuonMomentum.data());
        TrueNoFSICosOpeningAngleProtonsPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSICosOpeningAngleProtonsPlot",";cos(#theta_{#vec{p}_{L},#vec{p}_{R}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueNoFSICosOpeningAngleMuonTotalProtonPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSICosOpeningAngleMuonTotalProtonPlot",";cos(#theta_{#vec{p}_{#mu},#vec{p}_{sum}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueNoFSITransverseMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSITransverseMomentumPlot",";#delta P_{T} [GeV/c]",NBinsTransverseMomentum,ArrayNBinsTransverseMomentum.data());
        TrueNoFSIDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIDeltaAlphaTPlot",";#delta #alpha_{T} [deg]",NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT.data());

        // Before final state - GKI
        TrueNoFSICosOpeningAngleMomentumTransferTotalProtonPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSICosOpeningAngleMomentumTransferTotalProtonPlot",";cos(#theta_{#vec{q},#vec{p}_{sum}})",NBinsAngle,ArrayNBinsAngle.data());
        TrueNoFSIMissingMomentumPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIMissingMomentumPlot",";p_{n} [GeV/c]",NBinsMissingMomentum,ArrayNBinsMissingMomentum.data());
        TrueNoFSIAlphaThreeDPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueNoFSIAlphaThreeDPlot",";#alpha_{3D} [deg]",NBinsAlphaThreeD,ArrayNBinsAlphaThreeD.data());

	//Before final state - additional var
 	TrueNoFSIInvariantMassPlot[inte] = new TH1D(InteractionLabels[inte] + "TrueNoFSIInvariantMassPlot", ";W [GeV]", NBinsInvariantMass, ArrayNBinsInvariantMass.data());
        TrueNoFSICosOpeningAngleLProtonMuonPlot[inte] = new TH1D(InteractionLabels[inte] + "TrueNoFSICosOpeningAngleLProtonMuonPlot", ";cos(#theta_{#vec{p}_{L},#vec{p}_{#mu}})", NBinsCosAngleLPMu, ArrayNBinsCosAngleLPMu.data());
        TrueNoFSICosOpeningAngleRProtonMuonPlot[inte] = new TH1D(InteractionLabels[inte] + "TrueNoFSICosOpeningAngleRProtonMuonPlot", ";cos(#theta_{#vec{p}_{R},#vec{p}_{#mu}})", NBinsCosAngleLPMu, ArrayNBinsCosAngleLPMu.data());
        
	// Double differential final state
        TrueSerialTransverseMomentum_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialTransverseMomentum_InMuonCosThetaPlot",
            LabelXAxisTwoDTransverseMomentumInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices)[0]
        );
        TrueSerialDeltaAlphaT_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialDeltaAlphaT_InMuonCosThetaPlot",
            LabelXAxisTwoDDeltaAlphaTInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]
        );
        TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot",
            LabelXAxisTwoDCosOpeningAngleProtonsInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices)[0]
        );
        TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot",
            LabelXAxisTwoDCosOpeningMuonTotalProtonInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices)[0]
        );

        // Double differential final state - GKI
        TrueSerialMissingMomentum_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialMissingMomentum_InMuonCosThetaPlot",
            LabelXAxisTwoDMissingMomentumInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices)[0]
        );
        TrueSerialAlphaThreeD_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialAlphaThreeD_InMuonCosThetaPlot",
            LabelXAxisTwoDAlphaThreeDInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices)[0]
        );
        TrueSerialCosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialCosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot",
            LabelXAxisTwoDCosOpeningMomentumTransferTotalProtonInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices)[0]
        );

        // Double differential pre-FSI
        TrueSerialNoFSITransverseMomentum_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialNoFSITransverseMomentum_InMuonCosThetaPlot",
            LabelXAxisTwoDTransverseMomentumInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices)[0]
        );
        TrueSerialNoFSIDeltaAlphaT_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialNoFSIDeltaAlphaT_InMuonCosThetaPlot",
            LabelXAxisTwoDDeltaAlphaTInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices)[0]
        );
        TrueSerialNoFSICosOpeningAngleProtons_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialNoFSICosOpeningAngleProtons_InMuonCosThetaPlot",
            LabelXAxisTwoDCosOpeningAngleProtonsInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices)[0]
        );
        TrueSerialNoFSICosOpeningAngleMuonTotalProton_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialNoFSICosOpeningAngleMuonTotalProton_InMuonCosThetaPlot",
            LabelXAxisTwoDCosOpeningMuonTotalProtonInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices)[0]
        );

        // Double differential pre-FSI - GKI
        TrueSerialNoFSIMissingMomentum_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialNoFSIMissingMomentum_InMuonCosThetaPlot",
            LabelXAxisTwoDMissingMomentumInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices)[0]
        );
        TrueSerialNoFSIAlphaThreeD_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialNoFSIAlphaThreeD_InMuonCosThetaPlot",
            LabelXAxisTwoDAlphaThreeDInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices)[0]
        );
        TrueSerialNoFSICosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot[inte] = new TH1D(
            InteractionLabels[inte]+"TrueSerialNoFSICosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot",
            LabelXAxisTwoDCosOpeningMomentumTransferTotalProtonInMuonCosTheta,
            tools.Return2DNBins(TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices),
            &tools.Return2DBinIndices(TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices)[0]
        );

    //--------------------------------------------------//
    } // End of the loop over the interaction processes							

    //----------------------------------------//

    // Counters

    int CounterEventsPassedSelection = 0;
    int CounterQEEventsPassedSelection = 0;
    int CounterMECEventsPassedSelection = 0;
    int CounterRESEventsPassedSelection = 0;
    int CounterDISEventsPassedSelection = 0;
    int CounterCOHEventsPassedSelection = 0;

    // Counters for pre-FSI RES events
    int RESMode[24];
    for (int i = 0; i < 24; i++) {
        RESMode[i] = 0;
    }
    int NoFSICounterRESEventsPassedSelection = 0;

    //----------------------------------------//
    
    // Loop over the events

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        //----------------------------------------//	
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
        if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

        //----------------------------------------//	
            
        double weight = fScaleFactor*Units*A*Weight;
        if (fOutputFile.Contains("GiBUU")) { weight = weight/105.; } // To increase the stats, the GiBUU sample has been produced in 105 samples

        //----------------------------------------//	

        // Signal definition
        if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
        if (cc != 1) { continue; } // make sure that we have only CC interactions		

        // CC2p0pi event selection

        // Loop over final state particles
        int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0;
        int MuonTagging = 0, ElectronTagging = 0, PhotonTagging = 0;
        vector <int> ProtonID; ProtonID.clear();
        vector <int> MuonID; MuonID.clear();

        for (int i = 0; i < nfsp; i++) {
            double pf = TMath::Sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
            if (pdg[i] == 13 && (pf > 0.1 && pf < 1.2)) {
                MuonTagging++;
                MuonID.push_back(i);
            }
            if (pdg[i] == 2212 && (pf > 0.3 && pf < 1.)) {
                ProtonTagging++;
                ProtonID.push_back(i);
            }
            if (fabs(pdg[i]) == 211 && pf > 0.07)  {
                ChargedPionTagging++;
            }
            if (pdg[i] == 111)  {
                NeutralPionTagging++;
            }
            if (fabs(pdg[i]) == 11)  {
                ElectronTagging++;
            }
            if (fabs(pdg[i]) == 22)  {
                PhotonTagging++;
            }
        } // End of the loop over the final state particles

        // Tags for pre-FSI particles
        int NoFSIProtonTagging = 0, NoFSIChargedPionTagging = 0, NoFSINeutralPionTagging = 0;
        int NoFSIMuonTagging = 0, NoFSIElectronTagging = 0, NoFSIPhotonTagging = 0;
        vector <int> NoFSIProtonID; NoFSIProtonID.clear();
        vector <int> NoFSIMuonID; NoFSIMuonID.clear();

        // Loop over pre FSI variables for non-MEC files
        if (!fOutputFile.Contains("MEC")) {
            // Loop over pre FSI variables
            for (int i = 0; i < nvertp; i++) {
                double NoFSIpf = TMath::Sqrt(px_vert[i]*px_vert[i] + py_vert[i]*py_vert[i] + pz_vert[i]*pz_vert[i]);
                if (pdg_vert[i] == 13 && (NoFSIpf > 0.1 && NoFSIpf < 1.2)) {
                    NoFSIMuonTagging++;
                    NoFSIMuonID.push_back(i);
                }
                if (pdg_vert[i] == 2212 && (NoFSIpf > 0.3 && NoFSIpf < 1.)) {
                    NoFSIProtonTagging++;
                    NoFSIProtonID.push_back(i);
                }
                if (fabs(pdg_vert[i]) == 211 && NoFSIpf > 0.07)  {
                    NoFSIChargedPionTagging++;
                }
                if (pdg_vert[i] == 111)  {
                    NoFSINeutralPionTagging++;
                }
                if (fabs(pdg_vert[i]) == 11)  {
                    NoFSIElectronTagging++;
                }
                if (fabs(pdg_vert[i]) == 22)  {
                    NoFSIPhotonTagging++;
                }
            } // End of loop over pre-FSI particles
        }


        // Check if signal definition for final state particles is satisfied
        if (
            ProtonTagging == 2 &&
            ChargedPionTagging == 0 &&
            NeutralPionTagging == 0 && 
            MuonTagging == 1
        ) {
            //----------------------------------------//	

            // https://arxiv.org/pdf/2106.15809.pdf

            CounterEventsPassedSelection++;
            
            // Classify the events based on the interaction type

            int genie_mode = -1.;
            if (TMath::Abs(Mode) == 1) { CounterQEEventsPassedSelection++; genie_mode = 1; } // QE
            else if (TMath::Abs(Mode) == 2) { CounterMECEventsPassedSelection++; genie_mode = 2; } // MEC
            else if (
                TMath::Abs(Mode) == 10 ||
                TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
                TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
                ) { CounterRESEventsPassedSelection++; genie_mode = 3; } // RES
            else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { CounterDISEventsPassedSelection++; genie_mode = 4; } // DIS
            else if (TMath::Abs(Mode) == 16) { CounterCOHEventsPassedSelection++; genie_mode = 5;} // COH
            else { continue; } 

            // Feb 8 2022: Only case that is not covered is 15 = diffractive

            //----------------------------------------//

            // Create momentum vectors and helper
            TVector3 Muon(px[MuonID[0]], py[MuonID[0]], pz[MuonID[0]]);
            TVector3 LeadingProton(px[ProtonID[0]], py[ProtonID[0]], pz[ProtonID[0]]);
            TVector3 RecoilProton(px[ProtonID[1]], py[ProtonID[1]], pz[ProtonID[1]]);
            TwoPTools Helper(Muon, LeadingProton, RecoilProton);

            // Get variables of interest
            double LeadingProtonCosTheta = Helper.ReturnLeadingProtonCosTheta();
            double RecoilProtonCosTheta = Helper.ReturnRecoilProtonCosTheta();
            double LeadingProtonMomentum = Helper.ReturnLeadingProtonMomentum();
            double RecoilProtonMomentum = Helper.ReturnRecoilProtonMomentum();
            double MuonMomentum = Helper.ReturnMuonMomentum();
            double CosOpeningAngleProtons = Helper.ReturnCosOpeningAngleProtons();
            double CosOpeningAngleMuonTotalProton = Helper.ReturnCosOpeningAngleMuonTotalProton();
            double TransverseMomentum = Helper.ReturnTransverseMomentum();
            double DeltaAlphaT = Helper.ReturnDeltaAlphaT();
            double CosOpeningAngleMomentumTransferTotalProton = Helper.ReturnCosOpeningAngleMomentumTransferTotalProton();
            double MissingMomentum = Helper.ReturnMissingMomentum();
            double AlphaThreeD = Helper.ReturnAlphaThreeD();
	    double InvariantMass = Helper.ReturnInvariantMass();
	    double CosOpeningAngleLProtonMuon = Helper.ReturnCosOpeningAngleLProtonMuon();
	    double CosOpeningAngleRProtonMuon = Helper.ReturnCosOpeningAngleRProtonMuon();
            //----------------------------------------//

            // Filling in the histo regardless of interaction mode

            TRandom* rd = new TRandom();
            rd->SetSeed(ientry);				
            double Vx = rd->Uniform(ArrayNBinsVertexX[0],ArrayNBinsVertexX[NBinsVertexX]);
            double Vy = rd->Uniform(ArrayNBinsVertexY[0],ArrayNBinsVertexY[NBinsVertexY]);
            double Vz = rd->Uniform(ArrayNBinsVertexZ[0],ArrayNBinsVertexZ[NBinsVertexZ]);
            TrueVertexXPlot[0]->Fill(Vx,weight);
            TrueVertexYPlot[0]->Fill(Vy,weight);
            TrueVertexZPlot[0]->Fill(Vz,weight);

            TrueMuonCosThetaPlot[0]->Fill(CosLep,weight);
            TrueLeadingProtonCosThetaPlot[0]->Fill(LeadingProtonCosTheta,weight);
            TrueRecoilProtonCosThetaPlot[0]->Fill(RecoilProtonCosTheta,weight);
            TrueLeadingProtonMomentumPlot[0]->Fill(LeadingProtonMomentum,weight);
            TrueRecoilProtonMomentumPlot[0]->Fill(RecoilProtonMomentum,weight);
            TrueMuonMomentumPlot[0]->Fill(MuonMomentum,weight);
            TrueCosOpeningAngleProtonsPlot[0]->Fill(CosOpeningAngleProtons,weight);
            TrueCosOpeningAngleMuonTotalProtonPlot[0]->Fill(CosOpeningAngleMuonTotalProton,weight);
            TrueTransverseMomentumPlot[0]->Fill(TransverseMomentum,weight);
            TrueDeltaAlphaTPlot[0]->Fill(DeltaAlphaT,weight);
            TrueMissingMomentumPlot[0]->Fill(MissingMomentum,weight);
            TrueAlphaThreeDPlot[0]->Fill(AlphaThreeD,weight);
            TrueCosOpeningAngleMomentumTransferTotalProtonPlot[0]->Fill(CosOpeningAngleMomentumTransferTotalProton,weight);
	    TrueInvariantMassPlot[0]->Fill(InvariantMass,weight);
	    TrueCosOpeningAngleLProtonMuonPlot[0]->Fill(CosOpeningAngleLProtonMuon,weight);
	    TrueCosOpeningAngleRProtonMuonPlot[0]->Fill(CosOpeningAngleRProtonMuon,weight);

            //----------------------------------------//

            // Filling in the histo based on the interaction mode

            TrueVertexXPlot[genie_mode]->Fill(Vx,weight);
            TrueVertexYPlot[genie_mode]->Fill(Vy,weight);
            TrueVertexZPlot[genie_mode]->Fill(Vz,weight);
            TrueMuonCosThetaPlot[genie_mode]->Fill(CosLep,weight);
            TrueLeadingProtonCosThetaPlot[genie_mode]->Fill(LeadingProtonCosTheta,weight);
            TrueRecoilProtonCosThetaPlot[genie_mode]->Fill(RecoilProtonCosTheta,weight);
            TrueLeadingProtonMomentumPlot[genie_mode]->Fill(LeadingProtonMomentum,weight);
            TrueRecoilProtonMomentumPlot[genie_mode]->Fill(RecoilProtonMomentum,weight);
            TrueMuonMomentumPlot[genie_mode]->Fill(MuonMomentum,weight);
            TrueCosOpeningAngleProtonsPlot[genie_mode]->Fill(CosOpeningAngleProtons,weight);
            TrueCosOpeningAngleMuonTotalProtonPlot[genie_mode]->Fill(CosOpeningAngleMuonTotalProton,weight);
            TrueTransverseMomentumPlot[genie_mode]->Fill(TransverseMomentum,weight);
            TrueDeltaAlphaTPlot[genie_mode]->Fill(DeltaAlphaT,weight);
            TrueMissingMomentumPlot[genie_mode]->Fill(MissingMomentum,weight);
            TrueAlphaThreeDPlot[genie_mode]->Fill(AlphaThreeD,weight);
            TrueCosOpeningAngleMomentumTransferTotalProtonPlot[genie_mode]->Fill(CosOpeningAngleMomentumTransferTotalProton,weight);
	    TrueInvariantMassPlot[genie_mode]->Fill(InvariantMass,weight);
            TrueCosOpeningAngleLProtonMuonPlot[genie_mode]->Fill(CosOpeningAngleLProtonMuon,weight);
            TrueCosOpeningAngleRProtonMuonPlot[genie_mode]->Fill(CosOpeningAngleRProtonMuon,weight);
            //----------------------------------------//

            // Check for underflow/overflow for double differential plots
            if (TransverseMomentum < TwoDArrayTransverseMomentum[0]) { TransverseMomentum = (TwoDArrayTransverseMomentum[0] + TwoDArrayTransverseMomentum[1]) / 2.; }
            else if (TransverseMomentum > TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum]) { TransverseMomentum = (TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum] + TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum - 1]) / 2.; }

            if (DeltaAlphaT < TwoDArrayDeltaAlphaT[0]) { DeltaAlphaT = (TwoDArrayDeltaAlphaT[0] + TwoDArrayDeltaAlphaT[1]) / 2.; }
            else if (DeltaAlphaT > TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT]) { DeltaAlphaT = (TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT] + TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT - 1]) / 2.; }

            // GKI
            if (MissingMomentum < TwoDArrayMissingMomentum[0]) { MissingMomentum = (TwoDArrayMissingMomentum[0] + TwoDArrayMissingMomentum[1]) / 2.; }
            else if (MissingMomentum > TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum]) { MissingMomentum = (TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum] + TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum - 1]) / 2.; }

            if (AlphaThreeD < TwoDArrayAlphaThreeD[0]) { AlphaThreeD = (TwoDArrayAlphaThreeD[0] + TwoDArrayAlphaThreeD[1]) / 2.; }
            else if (AlphaThreeD > TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD]) { AlphaThreeD = (TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD] + TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD - 1]) / 2.; }

            // Double differential indices
            int MuonCosThetaTwoDIndex = tools.ReturnIndex(CosLep, TwoDArrayNBinsMuonCosTheta);
            int SerialTransverseMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                TransverseMomentum
            );
            int SerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                DeltaAlphaT
            );
            int SerialCosOpeningAngleProtonsInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                CosOpeningAngleProtons
            );
            int SerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                CosOpeningAngleMuonTotalProton
            );

            // Double differential indices - GKI
            int SerialMissingMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                MissingMomentum
            );
            int SerialAlphaThreeDInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                AlphaThreeD
            );
            int SerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices,
                MuonCosThetaTwoDIndex,
                CosOpeningAngleMomentumTransferTotalProton
            );

            // Fill in the histograms for double differential plots
            TrueSerialTransverseMomentum_InMuonCosThetaPlot[0]->Fill(SerialTransverseMomentumInMuonCosThetaIndex,weight);
            TrueSerialDeltaAlphaT_InMuonCosThetaPlot[0]->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
            TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot[0]->Fill(SerialCosOpeningAngleProtonsInMuonCosThetaIndex,weight);
            TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot[0]->Fill(SerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex,weight);

            TrueSerialTransverseMomentum_InMuonCosThetaPlot[genie_mode]->Fill(SerialTransverseMomentumInMuonCosThetaIndex,weight);
            TrueSerialDeltaAlphaT_InMuonCosThetaPlot[genie_mode]->Fill(SerialDeltaAlphaTInMuonCosThetaIndex,weight);
            TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot[genie_mode]->Fill(SerialCosOpeningAngleProtonsInMuonCosThetaIndex,weight);
            TrueSerialCosOpeningAngleMuonTotalProton_InMuonCosThetaPlot[genie_mode]->Fill(SerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex,weight);

            // GKI
            TrueSerialMissingMomentum_InMuonCosThetaPlot[0]->Fill(SerialMissingMomentumInMuonCosThetaIndex,weight);
            TrueSerialAlphaThreeD_InMuonCosThetaPlot[0]->Fill(SerialAlphaThreeDInMuonCosThetaIndex,weight);
            //TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot[0]->Fill(SerialCosOpeningAngleProtonsInMuonCosThetaIndex,weight);
            TrueSerialCosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot[0]->Fill(SerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex,weight);

            TrueSerialMissingMomentum_InMuonCosThetaPlot[genie_mode]->Fill(SerialMissingMomentumInMuonCosThetaIndex,weight);
            TrueSerialAlphaThreeD_InMuonCosThetaPlot[genie_mode]->Fill(SerialAlphaThreeDInMuonCosThetaIndex,weight);
            //TrueSerialCosOpeningAngleProtons_InMuonCosThetaPlot[genie_mode]->Fill(SerialCosOpeningAngleProtonsInMuonCosThetaIndex,weight);
            TrueSerialCosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot[genie_mode]->Fill(SerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex,weight);

            //----------------------------------------//
        } // End final state particles check

        // Check if signal definition for pre-FSI particles is satisfied
        else if (
            NoFSIProtonTagging == 2 &&
            NoFSIChargedPionTagging == 0 &&
            NoFSINeutralPionTagging == 0 && 
            NoFSIMuonTagging == 1 &&
            !fOutputFile.Contains("MEC") // pre-FSI variables only for non-MEC files
        ) {
            // Classify events based on interaction type
            int NoFSIgenie_mode = -1.;
            if (TMath::Abs(Mode) == 1) { NoFSIgenie_mode = 1; } // QE
            else if (TMath::Abs(Mode) == 2) { NoFSIgenie_mode = 2; } // MEC
            else if (
                TMath::Abs(Mode) == 10 ||
                TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
                TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
                ) { 
                    NoFSIgenie_mode = 3; 
                    if (fOutputFile.Contains("GiBUU")) { RESMode[Mode]++; NoFSICounterRESEventsPassedSelection++; }
                } // RES
            else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { NoFSIgenie_mode = 4; } // DIS
            else if (TMath::Abs(Mode) == 16) { NoFSIgenie_mode = 5;} // COH
            else { continue; } 

            // Create momentum vectors and helper
            TVector3 NoFSIMuon(px_vert[NoFSIMuonID[0]], py_vert[NoFSIMuonID[0]], pz_vert[NoFSIMuonID[0]]);
            TVector3 NoFSILeadingProton(px_vert[NoFSIProtonID[0]], py_vert[NoFSIProtonID[0]], pz_vert[NoFSIProtonID[0]]);
            TVector3 NoFSIRecoilProton(px_vert[NoFSIProtonID[1]], py_vert[NoFSIProtonID[1]], pz_vert[NoFSIProtonID[1]]);
            TwoPTools NoFSIHelper(NoFSIMuon, NoFSILeadingProton, NoFSIRecoilProton);

            // Get variables of interest
            double NoFSILeadingProtonCosTheta = NoFSIHelper.ReturnLeadingProtonCosTheta();
            double NoFSIRecoilProtonCosTheta = NoFSIHelper.ReturnRecoilProtonCosTheta();
            double NoFSILeadingProtonMomentum = NoFSIHelper.ReturnLeadingProtonMomentum();
            double NoFSIRecoilProtonMomentum = NoFSIHelper.ReturnRecoilProtonMomentum();
            double NoFSIMuonMomentum = NoFSIHelper.ReturnMuonMomentum();
            double NoFSICosOpeningAngleProtons = NoFSIHelper.ReturnCosOpeningAngleProtons();
            double NoFSICosOpeningAngleMuonTotalProton = NoFSIHelper.ReturnCosOpeningAngleMuonTotalProton();
            double NoFSITransverseMomentum = NoFSIHelper.ReturnTransverseMomentum();
            double NoFSIDeltaAlphaT = NoFSIHelper.ReturnDeltaAlphaT();
            double NoFSICosOpeningAngleMomentumTransferTotalProton = NoFSIHelper.ReturnCosOpeningAngleMomentumTransferTotalProton();
            double NoFSIMissingMomentum = NoFSIHelper.ReturnMissingMomentum();
            double NoFSIAlphaThreeD = NoFSIHelper.ReturnAlphaThreeD();
	    double NoFSIInvariantMass = NoFSIHelper.ReturnInvariantMass();
            double NoFSICosOpeningAngleLProtonMuon = NoFSIHelper.ReturnCosOpeningAngleLProtonMuon();
            double NoFSICosOpeningAngleRProtonMuon = NoFSIHelper.ReturnCosOpeningAngleRProtonMuon();
            
	    // Filling in the histo regardless of interaction mode
            TrueNoFSILeadingProtonCosThetaPlot[0]->Fill(NoFSILeadingProtonCosTheta,weight);
            TrueNoFSIRecoilProtonCosThetaPlot[0]->Fill(NoFSIRecoilProtonCosTheta,weight);
            TrueNoFSILeadingProtonMomentumPlot[0]->Fill(NoFSILeadingProtonMomentum,weight);
            TrueNoFSIRecoilProtonMomentumPlot[0]->Fill(NoFSIRecoilProtonMomentum,weight);
            TrueNoFSIMuonMomentumPlot[0]->Fill(NoFSIMuonMomentum,weight);
            TrueNoFSICosOpeningAngleProtonsPlot[0]->Fill(NoFSICosOpeningAngleProtons,weight);
            TrueNoFSICosOpeningAngleMuonTotalProtonPlot[0]->Fill(NoFSICosOpeningAngleMuonTotalProton,weight);
            TrueNoFSITransverseMomentumPlot[0]->Fill(NoFSITransverseMomentum,weight);
            TrueNoFSIDeltaAlphaTPlot[0]->Fill(NoFSIDeltaAlphaT,weight);
            TrueNoFSICosOpeningAngleMomentumTransferTotalProtonPlot[0]->Fill(NoFSICosOpeningAngleMomentumTransferTotalProton,weight);
            TrueNoFSIMissingMomentumPlot[0]->Fill(NoFSIMissingMomentum,weight);
            TrueNoFSIAlphaThreeDPlot[0]->Fill(NoFSIAlphaThreeD,weight);
	    TrueNoFSIInvariantMassPlot[0]->Fill(NoFSIInvariantMass,weight);
            TrueNoFSICosOpeningAngleLProtonMuonPlot[0]->Fill(NoFSICosOpeningAngleLProtonMuon,weight);
            TrueNoFSICosOpeningAngleRProtonMuonPlot[0]->Fill(NoFSICosOpeningAngleRProtonMuon,weight);
            
	    // Filling in the histo based on the interaction mode
            TrueNoFSILeadingProtonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSILeadingProtonCosTheta,weight);
            TrueNoFSIRecoilProtonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSIRecoilProtonCosTheta,weight);
            TrueNoFSILeadingProtonMomentumPlot[NoFSIgenie_mode]->Fill(NoFSILeadingProtonMomentum,weight);
            TrueNoFSIRecoilProtonMomentumPlot[NoFSIgenie_mode]->Fill(NoFSIRecoilProtonMomentum,weight);
            TrueNoFSIMuonMomentumPlot[NoFSIgenie_mode]->Fill(NoFSIMuonMomentum,weight);
            TrueNoFSICosOpeningAngleProtonsPlot[NoFSIgenie_mode]->Fill(NoFSICosOpeningAngleProtons,weight);
            TrueNoFSICosOpeningAngleMuonTotalProtonPlot[NoFSIgenie_mode]->Fill(NoFSICosOpeningAngleMuonTotalProton,weight);
            TrueNoFSITransverseMomentumPlot[NoFSIgenie_mode]->Fill(NoFSITransverseMomentum,weight);
            TrueNoFSIDeltaAlphaTPlot[NoFSIgenie_mode]->Fill(NoFSIDeltaAlphaT,weight);
            TrueNoFSICosOpeningAngleMomentumTransferTotalProtonPlot[NoFSIgenie_mode]->Fill(NoFSICosOpeningAngleMomentumTransferTotalProton,weight);
            TrueNoFSIMissingMomentumPlot[NoFSIgenie_mode]->Fill(NoFSIMissingMomentum,weight);
            TrueNoFSIAlphaThreeDPlot[NoFSIgenie_mode]->Fill(NoFSIAlphaThreeD,weight);
	    TrueNoFSIInvariantMassPlot[NoFSIgenie_mode]->Fill(NoFSIInvariantMass,weight);
            TrueNoFSICosOpeningAngleLProtonMuonPlot[NoFSIgenie_mode]->Fill(NoFSICosOpeningAngleLProtonMuon,weight);
            TrueNoFSICosOpeningAngleRProtonMuonPlot[NoFSIgenie_mode]->Fill(NoFSICosOpeningAngleRProtonMuon,weight);	    

            // Check for underflow/overflow for double differential plots
            if (NoFSITransverseMomentum < TwoDArrayTransverseMomentum[0]) { NoFSITransverseMomentum = (TwoDArrayTransverseMomentum[0] + TwoDArrayTransverseMomentum[1]) / 2.; }
            else if (NoFSITransverseMomentum > TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum]) { NoFSITransverseMomentum = (TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum] + TwoDArrayTransverseMomentum[TwoDNBinsTransverseMomentum - 1]) / 2.; }

            if (NoFSIDeltaAlphaT < TwoDArrayDeltaAlphaT[0]) { NoFSIDeltaAlphaT = (TwoDArrayDeltaAlphaT[0] + TwoDArrayDeltaAlphaT[1]) / 2.; }
            else if (NoFSIDeltaAlphaT > TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT]) { NoFSIDeltaAlphaT = (TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT] + TwoDArrayDeltaAlphaT[TwoDNBinsDeltaAlphaT - 1]) / 2.; }
            
            if (NoFSIMissingMomentum < TwoDArrayMissingMomentum[0]) { NoFSIMissingMomentum = (TwoDArrayMissingMomentum[0] + TwoDArrayMissingMomentum[1]) / 2.; }
            else if (NoFSIMissingMomentum > TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum]) { NoFSIMissingMomentum = (TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum] + TwoDArrayMissingMomentum[TwoDNBinsMissingMomentum - 1]) / 2.; }

            if (NoFSIAlphaThreeD < TwoDArrayAlphaThreeD[0]) { NoFSIAlphaThreeD = (TwoDArrayAlphaThreeD[0] + TwoDArrayAlphaThreeD[1]) / 2.; }
            else if (NoFSIAlphaThreeD > TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD]) { NoFSIAlphaThreeD = (TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD] + TwoDArrayAlphaThreeD[TwoDNBinsAlphaThreeD - 1]) / 2.; }

            // Double differential indices
            int NoFSIMuonCosThetaTwoDIndex = tools.ReturnIndex(CosLep, TwoDArrayNBinsMuonCosTheta);
            int NoFSISerialTransverseMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsTransverseMomentumInMuonCosThetaSlices,
                NoFSIMuonCosThetaTwoDIndex,
                NoFSITransverseMomentum
            );
            int NoFSISerialDeltaAlphaTInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsDeltaAlphaTInMuonCosThetaSlices,
                NoFSIMuonCosThetaTwoDIndex,
                NoFSIDeltaAlphaT
            );
            int NoFSISerialCosOpeningAngleProtonsInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsCosOpeningAngleProtonsInMuonCosThetaSlices,
                NoFSIMuonCosThetaTwoDIndex,
                NoFSICosOpeningAngleProtons
            );
            int NoFSISerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsCosOpeningAngleMuonTotalProtonInMuonCosThetaSlices,
                NoFSIMuonCosThetaTwoDIndex,
                NoFSICosOpeningAngleMuonTotalProton
            );
            // GKI
            int NoFSIMuonCosThetaThreeDIndex = tools.ReturnIndex(CosLep, TwoDArrayNBinsMuonCosTheta);
            int NoFSISerialMissingMomentumInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsMissingMomentumInMuonCosThetaSlices,
                NoFSIMuonCosThetaThreeDIndex,
                NoFSIMissingMomentum
            );
            int NoFSISerialAlphaThreeDInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsAlphaThreeDInMuonCosThetaSlices,
                NoFSIMuonCosThetaThreeDIndex,
                NoFSIAlphaThreeD
            );
            int NoFSISerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex = tools.ReturnIndexIn2DList(
                TwoDArrayNBinsCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaSlices,
                NoFSIMuonCosThetaThreeDIndex,
                NoFSICosOpeningAngleMomentumTransferTotalProton
            );

            // Fill in the histograms for double differential plots
            TrueSerialNoFSITransverseMomentum_InMuonCosThetaPlot[0]->Fill(NoFSISerialTransverseMomentumInMuonCosThetaIndex,weight);
            TrueSerialNoFSIDeltaAlphaT_InMuonCosThetaPlot[0]->Fill(NoFSISerialDeltaAlphaTInMuonCosThetaIndex,weight);
            TrueSerialNoFSICosOpeningAngleProtons_InMuonCosThetaPlot[0]->Fill(NoFSISerialCosOpeningAngleProtonsInMuonCosThetaIndex,weight);
            TrueSerialNoFSICosOpeningAngleMuonTotalProton_InMuonCosThetaPlot[0]->Fill(NoFSISerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex,weight);

            TrueSerialNoFSITransverseMomentum_InMuonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSISerialTransverseMomentumInMuonCosThetaIndex,weight);
            TrueSerialNoFSIDeltaAlphaT_InMuonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSISerialDeltaAlphaTInMuonCosThetaIndex,weight);
            TrueSerialNoFSICosOpeningAngleProtons_InMuonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSISerialCosOpeningAngleProtonsInMuonCosThetaIndex,weight);
            TrueSerialNoFSICosOpeningAngleMuonTotalProton_InMuonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSISerialCosOpeningAngleMuonTotalProtonInMuonCosThetaIndex,weight);
            TrueSerialNoFSIMissingMomentum_InMuonCosThetaPlot[0]->Fill(NoFSISerialMissingMomentumInMuonCosThetaIndex,weight);
            TrueSerialNoFSIAlphaThreeD_InMuonCosThetaPlot[0]->Fill(NoFSISerialAlphaThreeDInMuonCosThetaIndex,weight);
            TrueSerialNoFSICosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot[0]->Fill(NoFSISerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex,weight);

            TrueSerialNoFSIMissingMomentum_InMuonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSISerialMissingMomentumInMuonCosThetaIndex,weight);
            TrueSerialNoFSIAlphaThreeD_InMuonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSISerialAlphaThreeDInMuonCosThetaIndex,weight);
            TrueSerialNoFSICosOpeningAngleMomentumTransferTotalProton_InMuonCosThetaPlot[NoFSIgenie_mode]->Fill(NoFSISerialCosOpeningAngleMomentumTransferTotalProtonInMuonCosThetaIndex,weight);

        } // End pre-FSI particles check
        else { continue; }
    } // End of the loop over the events

    //----------------------------------------//

    std::cout << "Percentage of events passing the selection cuts = " << 
    double(CounterEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting QE events = " << 
    double(CounterQEEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting MEC events = " << 
    double(CounterMECEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting RES events = " << 
    double(CounterRESEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting DIS events = " << 
    double(CounterDISEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

    std::cout << "Success percentage in selecting COH events = " << 
    double(CounterCOHEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;	

    if (fOutputFile.Contains("GiBUU")) {
        std::cout << "Mode percentage for RES events" << std::endl;
        std::cout << "10: " << double(RESMode[10]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "11: " << double(RESMode[11]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "12: " << double(RESMode[12]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "13: " << double(RESMode[13]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "17: " << double(RESMode[17]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "22: " << double(RESMode[22]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
        std::cout << "23: " << double(RESMode[23]) / double(NoFSICounterRESEventsPassedSelection)*100. << " %" << std::endl;
    }

    //----------------------------------------//	
    //----------------------------------------//	

    // Division by bin width to get the cross sections	
    // Loop over the interaction processes

    for (int inte = 0; inte < NInte; inte++) {

        //----------------------------------------//
    
        Reweight(TrueVertexXPlot[inte]);
        Reweight(TrueVertexYPlot[inte]);
        Reweight(TrueVertexZPlot[inte]);
        Reweight(TrueMuonCosThetaPlot[inte]);
        Reweight(TrueLeadingProtonCosThetaPlot[inte]);
        Reweight(TrueRecoilProtonCosThetaPlot[inte]);
        Reweight(TrueLeadingProtonMomentumPlot[inte]);
        Reweight(TrueRecoilProtonMomentumPlot[inte]);
        Reweight(TrueMuonMomentumPlot[inte]);
        Reweight(TrueCosOpeningAngleProtonsPlot[inte]);
        Reweight(TrueCosOpeningAngleMuonTotalProtonPlot[inte]);
        Reweight(TrueTransverseMomentumPlot[inte]);
        Reweight(TrueDeltaAlphaTPlot[inte]);
	Reweight(TrueInvariantMassPlot[inte]);
	Reweight(TrueCosOpeningAngleLProtonMuonPlot[inte]);
	Reweight(TrueCosOpeningAngleRProtonMuonPlot[inte]);

        Reweight(TrueNoFSILeadingProtonCosThetaPlot[inte]);
        Reweight(TrueNoFSIRecoilProtonCosThetaPlot[inte]);
        Reweight(TrueNoFSILeadingProtonMomentumPlot[inte]);
        Reweight(TrueNoFSIRecoilProtonMomentumPlot[inte]);
        Reweight(TrueNoFSIMuonMomentumPlot[inte]);
        Reweight(TrueNoFSICosOpeningAngleProtonsPlot[inte]);
        Reweight(TrueNoFSICosOpeningAngleMuonTotalProtonPlot[inte]);
        Reweight(TrueNoFSITransverseMomentumPlot[inte]);
        Reweight(TrueNoFSIDeltaAlphaTPlot[inte]);
	Reweight(TrueNoFSIInvariantMassPlot[inte]);
        Reweight(TrueNoFSICosOpeningAngleLProtonMuonPlot[inte]);
        Reweight(TrueNoFSICosOpeningAngleRProtonMuonPlot[inte]);
        
	// GKI
        Reweight(TrueCosOpeningAngleMomentumTransferTotalProtonPlot[inte]);
        Reweight(TrueMissingMomentumPlot[inte]);
        Reweight(TrueAlphaThreeDPlot[inte]);

        Reweight(TrueNoFSICosOpeningAngleMomentumTransferTotalProtonPlot[inte]);
        Reweight(TrueNoFSIMissingMomentumPlot[inte]);
        Reweight(TrueNoFSIAlphaThreeDPlot[inte]);


        //----------------------------------------//

    } // End of the loop over the interaction processes		

    //----------------------------------------//		
        
    file->cd();
    file->Write();
    fFile->Close();

    std::cout << std::endl;
    std::cout << "File " << FileNameAndPath +" has been created " << std::endl; 
    std::cout << std::endl;

    std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;

    //----------------------------------------//		

} // End of the program

//----------------------------------------//		

void Reweight(TH1D* h) {
    int NBins = h->GetXaxis()->GetNbins();

    for (int i = 0; i < NBins; i++) {

        double CurrentEntry = h->GetBinContent(i+1);
        double NewEntry = CurrentEntry / h->GetBinWidth(i+1);

        double CurrentError = h->GetBinError(i+1);
        double NewError = CurrentError / h->GetBinWidth(i+1);

        h->SetBinContent(i+1,NewEntry); 
        h->SetBinError(i+1,NewError); 
        //h->SetBinError(i+1,0.000001); 

    }
}
//----------------------------------------//		
