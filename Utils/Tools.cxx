// Class created by Afroditi Papadopoulou (apapadop@mit.edu)

//----------------------------------------//

#ifndef Tools_cxx
#define Tools_cxx

// STD includes
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <sstream>

// Root Includes
#include "TRandom3.h"

#include "Tools.h"

using namespace std;

//----------------------------------------//

std::vector< std::vector<double> > Tools::CollapseMatrixIntoArray(std::vector< std::vector< std::vector<double> > > Matrix) {

	// First make sure that the new vector has the correct size

	int NRows = Matrix.size();
	std::vector< std::vector<double> > SerialArray; SerialArray.resize(NRows);
	int NColumns = Matrix.at(0).size();

	for (int icolumn = 0; icolumn < NColumns; icolumn++) {

		for (int irow = 0; irow < NRows; irow++) {		

			for (int ielement = 0; ielement < (int)(Matrix.at(irow).at(icolumn).size()); ielement++) {

				SerialArray[irow].push_back(Matrix.at(irow).at(icolumn).at(ielement));

			}

		}

	}

	return SerialArray;

}

//----------------------------------------//

TH2D* Tools::Get2DHistoBins(TH2D* h,int LowBin,int HighBin,double ScaleFactor,std::vector<double> Binning, bool Scale = true) {

	TString PlotName = TString(h->GetName()) + "_"+TString(std::to_string(LowBin)) + "_"+TString(std::to_string(HighBin));
	TString XaxisTitle = h->GetXaxis()->GetTitle();
	TString YaxisTitle = h->GetYaxis()->GetTitle();	
	int NBins = HighBin - LowBin + 1;
 
	TH2D* clone = new TH2D(PlotName,";" + XaxisTitle + ";" + YaxisTitle, NBins, &Binning[0], NBins, &Binning[0] );

	// Loop over the bin indices, aka starting from 1
	for (int iBin = 1; iBin <= NBins; iBin++) {

		for (int jBin = 1; jBin <= NBins; jBin++) {

			double BinWidth = ( Binning.at(iBin) - Binning.at(iBin - 1) ) * ( Binning.at(jBin) - Binning.at(jBin - 1) );
			double TotalScaleFactor = BinWidth * TMath::Power(ScaleFactor,2.);	
			if (!Scale) { TotalScaleFactor = 1.; }		

			double CurrentHistoEntry = h->GetBinContent(LowBin + iBin - 1,LowBin + jBin - 1);		
			double CurrentHistoError = h->GetBinError(LowBin + iBin - 1,LowBin + jBin - 1);	

			double NewHistoEntry = CurrentHistoEntry / TotalScaleFactor;
			double NewHistoError = CurrentHistoError / TotalScaleFactor;

			clone->SetBinContent(iBin,jBin,NewHistoEntry);
			clone->SetBinError(iBin,jBin,NewHistoError);

		}		

	} // End of the loop over the bin indices

	return clone;

}

//----------------------------------------//

TH1D* Tools::GetHistoBins(TH1D* h,int LowBin,int HighBin,double ScaleFactor,std::vector<double> Binning, TString Name) {

	TString PlotName = Name + TString(h->GetName()) + "_"+TString(std::to_string(LowBin)) + "_"+TString(std::to_string(HighBin));
	TString XaxisTitle = h->GetXaxis()->GetTitle();
	TString YaxisTitle = h->GetYaxis()->GetTitle();	
	int NBins = HighBin - LowBin + 1;

	TH1D* clone = new TH1D(PlotName,";" + XaxisTitle + ";" + YaxisTitle, NBins, &Binning[0] );

	// Loop over the bin indices, aka starting from 1
	for (int iBin = 1; iBin <= NBins; iBin++) {

		double BinWidth = Binning.at(iBin) - Binning.at(iBin - 1);
		double TotalScaleFactor = BinWidth * ScaleFactor;

		double CurrentHistoEntry = h->GetBinContent(LowBin + iBin - 1);		
		double CurrentHistoError = h->GetBinError(LowBin + iBin - 1);	

		double NewHistoEntry = CurrentHistoEntry / TotalScaleFactor;
		double NewHistoError = CurrentHistoError / TotalScaleFactor;		

		clone->SetBinContent(iBin,NewHistoEntry);
		clone->SetBinError(iBin,NewHistoError);		

	} // End of the loop over the bin indices

	return clone;

}

TH1D* Tools::GetHistoBinsNoScale(TH1D* h,int LowBin,int HighBin,std::vector<double> Binning, TString Name) {
	TString PlotName = Name + TString(h->GetName()) + "_"+TString(std::to_string(LowBin)) + "_"+TString(std::to_string(HighBin));
	TString XaxisTitle = h->GetXaxis()->GetTitle();
	TString YaxisTitle = h->GetYaxis()->GetTitle();
	int NBins = HighBin - LowBin + 1;

	TH1D* clone = new TH1D(PlotName,";" + XaxisTitle + ";" + YaxisTitle, NBins, &Binning[0] );

	// Loop over the bin indices, aka starting from 1
	for (int iBin = 1; iBin <= NBins; iBin++) {
		double CurrentHistoEntry = h->GetBinContent(LowBin + iBin - 1);		
		double CurrentHistoError = h->GetBinError(LowBin + iBin - 1);		
		clone->SetBinContent(iBin,CurrentHistoEntry);
		clone->SetBinError(iBin,CurrentHistoError);		
	} // End of the loop over the bin indices

	return clone;
}

//----------------------------------------//

std::vector<TMatrixD> Tools::MatrixDecomp(int nbins,TVectorD matrix_pred,TMatrixD matrix_syst) {

	// MiniBooNE note from Mike Schaevitz
	// https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=5926&filename=tn253.pdf&version=1
	
	TMatrixD matrix_shape(nbins, nbins);
	TMatrixD matrix_mixed(nbins, nbins);
	TMatrixD matrix_norm(nbins, nbins);

	///
	double N_T = 0;
	for (int idx = 0; idx < nbins; idx++) { N_T += matrix_pred(idx); }

	///
	double M_kl = 0;

	for (int i = 0; i < nbins; i++) {
		
		for (int j = 0; j < nbins; j++) {
			
			M_kl += matrix_syst(i,j);
	
		}

	}

	///
	for (int i = 0; i < nbins; i++) {

		for (int j = 0; j < nbins; j++) {	
  
			double N_i = matrix_pred(i);
			double N_j = matrix_pred(j);
			double M_ij = matrix_syst(i,j);	  
			double M_ik = 0; for(int k=0; k<nbins; k++) M_ik += matrix_syst(i,k);
			double M_kj = 0; for(int k=0; k<nbins; k++) M_kj += matrix_syst(k,j);
			matrix_shape(i,j) = M_ij - N_j*M_ik/N_T - N_i*M_kj/N_T + N_i*N_j*M_kl/N_T/N_T;
			matrix_mixed(i,j) = N_j*M_ik/N_T + N_i*M_kj/N_T - 2*N_i*N_j*M_kl/N_T/N_T;	
			matrix_norm(i,j) = N_i*N_j*M_kl/N_T/N_T;

			// debug
			//if (i == j) {
	
			//	cout << "matrix_syst(" << i << "," << j <<") = " << matrix_syst(i,j) << " "; 
			//	cout << "matrix_norm(" << i << "," << j <<") = " << matrix_norm(i,j) << " "; 
			//	cout << "matrix_shape(" << i << "," << j <<") = " << matrix_shape(i,j) << " "; 
			//	cout << "matrix_mixed(" << i << "," << j <<") = " << matrix_mixed(i,j) << endl; 

			//} // end of debugging

		}

	}

	//cout << endl;

	std::vector<TMatrixD> NormShapeVector = {matrix_norm,matrix_shape+matrix_mixed};
	return NormShapeVector;

}

//----------------------------------------//

int Tools::ReturnIndexIn3DList(std::vector< std::vector< std::vector<double> > > BinEdgeVector, int FirstSliceIndex, int SecondSliceIndex, double ValueInSlice) { 

	int BinIndex = 1; // TH1D bin index, thus starting from 1
	int VectorRowSize = BinEdgeVector.size();

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int VectorColumnSize = BinEdgeVector.at(irow).size();

		for (int icolumn = 0; icolumn < VectorColumnSize; icolumn++){

			if (irow != FirstSliceIndex || icolumn != SecondSliceIndex) {

				BinIndex += BinEdgeVector.at(irow).at(icolumn).size()-1;

			} else {

				int LocalBins = BinEdgeVector.at(irow).at(icolumn).size();
				BinIndex += ReturnIndex(ValueInSlice, BinEdgeVector.at(irow).at(icolumn));
				return BinIndex;

			}

		}	

	}

	return BinIndex+1; // Offset to account for bin number vs array index

}

//----------------------------------------//

std::vector<double> Tools::Return3DBinIndices(std::vector< std::vector< std::vector<double> > > BinEdgeVector) { 

	int BinCounter = 0;
	int VectorRowSize = BinEdgeVector.size();
	std::vector<double> BinIndices;

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int NElements = BinEdgeVector.at(irow).size();

		for (int ielement = 0; ielement < NElements; ielement++) {

			int NElementsColumn = BinEdgeVector.at(irow).at(ielement).size();	

			for (int icolumn = 0; icolumn < NElementsColumn-1; icolumn++) {

				// Lower bin edges in the form of indices
				// + 0.5 so that the bins are centered at an integer (e.g. Bin 1, 2, 3 et al)
				BinIndices.push_back(BinCounter+0.5);
				BinCounter++;

			}	

		}

	}
	// Upper bin edge
	BinIndices.push_back(BinCounter+0.5);
	return BinIndices;

}	

//----------------------------------------//

int Tools::Return3DNBins(std::vector< std::vector< std::vector<double> > > BinEdgeVector) { 

	int NBins = 0;
	int VectorRowSize = BinEdgeVector.size();

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int NElements = BinEdgeVector.at(irow).size();

		for (int icolumn = 0; icolumn < NElements; icolumn++) {

			int NElementsColumn = BinEdgeVector.at(irow).at(icolumn).size();

			// Number of bins for each subvector
			NBins += NElementsColumn-1;

		}

	}

	return NBins;

}

//----------------------------------------//

int Tools::ReturnIndexIn2DList(std::vector< std::vector<double> > BinEdgeVector, int SliceIndex, double ValueInSlice) { 

	int BinIndex = 1; // TH1D bin index, thus starting from 1
	int VectorRowSize = BinEdgeVector.size();

	for (int irow = 0; irow < VectorRowSize; irow++) {

		if (irow != SliceIndex) {

			BinIndex += BinEdgeVector.at(irow).size()-1;

		} else {

			int LocalBins = BinEdgeVector.at(irow).size();
			BinIndex += ReturnIndex(ValueInSlice, BinEdgeVector.at(irow));
			return BinIndex;

		}


	}

	return BinIndex+1; // Offset to account for bin number vs array index

}

//----------------------------------------//

std::vector<double> Tools::Return2DBinIndices(std::vector< std::vector<double> > BinEdgeVector) { 

	int BinCounter = 0;
	int VectorRowSize = BinEdgeVector.size();
	std::vector<double> BinIndices;

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int NElements = BinEdgeVector.at(irow).size();

		for (int ielement = 0; ielement < NElements-1; ielement++) {

			// Lower bin edges in the form of indices
			// + 0.5 so that the bins are centered at an integer (e.g. Bin 1, 2, 3 et al)
			BinIndices.push_back(BinCounter+0.5);
			BinCounter++;

		}

	}

	// Upper bin edge
	BinIndices.push_back(BinCounter+0.5);
	return BinIndices;

}	

//----------------------------------------//

int Tools::Return2DNBins(std::vector< std::vector<double> > BinEdgeVector) { 

	int NBins = 0;
	int VectorRowSize = BinEdgeVector.size();

	for (int irow = 0; irow < VectorRowSize; irow++) {

		int NElements = BinEdgeVector.at(irow).size();

		// Number of bins for each subvector
		NBins += NElements-1;

	}

	return NBins;

}	

//----------------------------------------//

int Tools::ConcatRunSubRunEvent(int run, int subrun, int event, int universe) const { 

	// Convert the integers to string 
	std::string srun    = std::to_string(run); 
	std::string ssubrun = std::to_string(subrun);
	std::string sevent  = std::to_string(event);
	std::string suniv  = std::to_string(universe);	

	// Concatenate the subrun and universe. Dont add the run/subrun because it makes the number too long for storing as an int
	std::string s =  ssubrun + suniv; 

	// std::cout << srun << "  " << ssubrun << "  " << sevent<< std::endl;
  
	// Convert the concatenated string to integer 
	int c = stoi(s); 
  
	// return the formed integer 
	return c; 

}

//----------------------------------------//

double Tools::PoissonRandomNumber(int seed) const {

	// Set the seed of the TRandom 3 based on the run,subrun,event
	TRandom3* rand = new TRandom3();
	rand->SetSeed(seed); 

	// Generate the weight, using a poisson dist with mean 1
	double weight_poisson = rand->Poisson(1);

	delete rand;
	return weight_poisson;

}

//----------------------------------------//

bool Tools::is_meson_or_antimeson( int pdg_code ) {

	// Ignore differences between mesons and antimesons for this test. Mesons
	// will have positive PDG codes, while antimesons will have negative ones.
	int abs_pdg = std::abs( pdg_code );

	// Meson PDG codes have no more than seven digits. Seven-digit
	// codes beginning with "99" are reserved for generator-specific
	// particles
	if ( abs_pdg >= 9900000 ) return false;

	// Mesons have a value of zero for $n_{q1}$, the thousands digit
	int thousands_digit = ( abs_pdg / 1000 ) % 10;
	if ( thousands_digit != 0 ) return false;

	// They also have a nonzero value for $n_{q2}$, the hundreds digit
	int hundreds_digit = ( abs_pdg / 100 ) % 10;
	if ( hundreds_digit == 0 ) return false;

	// Reserved codes for Standard Model parton distribution functions
	if ( abs_pdg >= 901 && abs_pdg <= 930 ) return false;

	// Reggeon and pomeron
	if ( abs_pdg == 110 || abs_pdg == 990 ) return false;

	// Reserved codes for GEANT tracking purposes
	if ( abs_pdg == 998 || abs_pdg == 999 ) return false;

	// Reserved code for generator-specific pseudoparticles
	if ( abs_pdg == 100 ) return false;

	// If we've passed all of the tests above, then the particle is a meson
	return true;
}

//----------------------------------------//

bool Tools::inFVVector(TVector3 vector) {

	if(vector.X() < (FVx - borderx) && (vector.X() > borderx) && (vector.Y() < (FVy/2. - bordery)) && (vector.Y() > (-FVy/2. + bordery)) && 
	(vector.Z() < (FVz - borderz)) && (vector.Z() > borderz)) return true;
	else return false;
}

//----------------------------------------//

bool Tools::inFV(double x, double y, double z) {

	if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
	else return false;
}

//----------------------------------------//

bool Tools::IsContained(TVector3 TrackStart, TVector3 TrackEnd) {

	bool fContainedTrack = false;

	if ( inFV(TrackStart.X(),TrackStart.Y(),TrackStart.Z()) && inFV(TrackEnd.X(),TrackEnd.Y(),TrackEnd.Z()) ) { fContainedTrack = true; }

	return fContainedTrack;

}

//----------------------------------------//

double Tools::PToKE(int pdg, double momentum) {

	double TrackKEConvert = -99;

	if ( fabs(pdg) == MuonPdg) { TrackKEConvert = sqrt( TMath::Power(momentum,2) + TMath::Power(MuonMass,2) ) - MuonMass ; }
	if ( fabs(pdg) == ProtonPdg) { TrackKEConvert = sqrt( TMath::Power(momentum,2) + TMath::Power(ProtonMass,2) ) - ProtonMass; }

	return TrackKEConvert;

}

//----------------------------------------//

double Tools::KEToP(int pdg, double ke) {

	double TrackPConvert = -99;

	if ( fabs(pdg) == MuonPdg) { TrackPConvert = sqrt( TMath::Power(ke+MuonMass,2) - TMath::Power(MuonMass,2) ) ; }
	if ( fabs(pdg) == ProtonPdg) { TrackPConvert = sqrt( TMath::Power(ke+ProtonMass,2) - TMath::Power(ProtonMass,2) ); }

	return TrackPConvert;

}

//----------------------------------------//

TString Tools::to_string_with_precision(double a_value, const int n = 2) {

    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());

}

//----------------------------------------//

TString Tools::ConvertToString(double value) {

	TString StringValue = Tools::to_string_with_precision(value, 2);
	StringValue.ReplaceAll(".","_");
	StringValue.ReplaceAll("-","Minus");	

	return StringValue;

}

//----------------------------------------//

int Tools::ReturnIndex(double value, std::vector<double> vec) {

	int length = vec.size();
	int index = -1;

	for (int i = 0; i < length-1; i ++) {

		if (i == 0 && value == vec.at(0)) { return 0; }
		if (value > vec.at(i) && value <= vec.at(i+1)) { return i; }

	}	

	return index;

}

//----------------------------------------//

void Tools::Reweight(TH1D* h, double SF = 1.) {

	int NBins = h->GetXaxis()->GetNbins();

	for (int i = 0; i < NBins; i++) {

		double CurrentEntry = h->GetBinContent(i+1);
		double NewEntry = CurrentEntry * SF / h->GetBinWidth(i+1);

		double CurrentError = h->GetBinError(i+1);
		double NewError = CurrentError * SF / h->GetBinWidth(i+1);

		h->SetBinContent(i+1,NewEntry); 
//		h->SetBinError(i+1,NewError); 
		h->SetBinError(i+1,0.000001); 

	}
}

//----------------------------------------//

void Tools::Unweight(TH1D* h, double SF = 1.) {
	int NBins = h->GetXaxis()->GetNbins();
	for (int i = 0; i < NBins; i++) {
		double CurrentEntry = h->GetBinContent(i+1);
		double NewEntry = CurrentEntry * SF * h->GetBinWidth(i+1);

		double CurrentError = h->GetBinError(i+1);
		double NewError = CurrentError * SF * h->GetBinWidth(i+1);

		h->SetBinContent(i+1,NewEntry); 
		h->SetBinError(i+1,0.000001); 
	}
}

//----------------------------------------//

void Tools::Reweight2D(TH2D* h, double SF = 1.) {

	int NBinsX = h->GetXaxis()->GetNbins();
	int NBinsY = h->GetYaxis()->GetNbins();

	for (int i = 0; i < NBinsX; i++) {

		for (int j = 0; j < NBinsX; j++) {

			double CurrentEntry = h->GetBinContent(i+1,j+1);
			double NewEntry = CurrentEntry * SF / ( h->GetXaxis()->GetBinWidth(i+1) * h->GetYaxis()->GetBinWidth(j+1) );

			double CurrentError = h->GetBinError(i+1,j+1);
			double NewError = CurrentError * SF / ( h->GetXaxis()->GetBinWidth(i+1) * h->GetYaxis()->GetBinWidth(j+1) );

			h->SetBinContent(i+1,j+1,NewEntry); 
//			h->SetBinError(i+1,j+1,NewError); 
			h->SetBinError(i+1,j+1,0.000001); 

		}

	}

}

//----------------------------------------//

std::tuple<int, vector<double>, vector<int>, vector<int>, vector<int>> Tools::FlattenNDBins(vector<double> SliceDiscriminators, vector<vector<double>> SliceBinning) {
	int NSlices = 1;
	vector<double> SerialVectorRanges;
	vector<int> SerialVectorBins;
	vector<int> SerialVectorLowBin;
	vector<int> SerialVectorHighBin;

	int BinCounter = 1;

	// For the discriminator, how many slices do we have?
	int SliceDiscrimSize = SliceDiscriminators.size() - 1;
	NSlices *= SliceDiscrimSize; 

	for (int iSliceDiscrimSize = 0; iSliceDiscrimSize < SliceDiscrimSize; iSliceDiscrimSize++) {
		// Accessing the vector<double> with the bin ranges
		int SliceDiscrimValue = SliceBinning.at(iSliceDiscrimSize).size();

		// Storing the number of bins for a specific slice					
		SerialVectorBins.push_back(SliceDiscrimValue-1);
		for (int iBin = 0; iBin < SliceDiscrimValue; iBin++) {
			double BinValue = SliceBinning.at(iSliceDiscrimSize).at(iBin);
			// First bin number for a given slice
			if (iBin == 0) { SerialVectorLowBin.push_back(BinCounter); }
			// Last bin number for a given slice
			if (iBin == SliceDiscrimValue-2) { SerialVectorHighBin.push_back(BinCounter); }	
			// Storing the binning for a specific slice
			SerialVectorRanges.push_back(BinValue);
			// Increase the global bin counter
			// But not for the last bin
			if (iBin != SliceDiscrimValue-1) { BinCounter++; }

		} // End of the loop over the bins of a given slice
	} // End of the loop over the slices of the discriminator

	return {NSlices, SerialVectorRanges, SerialVectorBins, SerialVectorLowBin, SerialVectorHighBin};
}

void Tools::CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval, double &sigma) {
	// Clone them so we can scale them 
	TH1D* h_model_clone = (TH1D*)h_model->Clone();
	TH1D* h_data_clone  = (TH1D*)h_data->Clone();
	TH2D* h_cov_clone   = (TH2D*)cov->Clone();
	int NBins = h_cov_clone->GetNbinsX();

	// Getting covariance matrix in TMatrix form
	TMatrixD cov_m;
	cov_m.Clear();
	cov_m.ResizeTo(NBins,NBins);

	// loop over rows
	for (int i = 0; i < NBins; i++) {			
		// loop over columns
		for (int j = 0; j < NBins; j++) {
			cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
		}
	}
	TMatrixD copy_cov_m = cov_m;

	// Inverting the covariance matrix
	TMatrixD inverse_cov_m = cov_m.Invert();

	// Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
	// x = data, mu = model, E^(-1) = inverted covariance matrix 
	chi = 0.;
	
	for (int i = 0; i < NBins; i++) {
		//double XWidth = h_data_clone->GetBinWidth(i+1);
		for (int j = 0; j < NBins; j++) {
			//double YWidth = h_data_clone->GetBinWidth(i+1);
			double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
			double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
			double LocalChi = diffi * inverse_cov_m[i][j] * diffj; 
			chi += LocalChi;
		}
	}
	ndof = h_data_clone->GetNbinsX();
	pval = TMath::Prob(chi, ndof);
	sigma = TMath::Sqrt( TMath::ChisquareQuantile( 1-pval, 1 ) ); 

	delete h_model_clone;
	delete h_data_clone;
	delete h_cov_clone;
}

//---------------------//

const std::vector<std::string> Tools::GetInputFiles(const std::string TargetPath, bool print = false, bool detvar = false) {

	std::vector<std::string> Input;

	for (const auto & entry : std::filesystem::directory_iterator(TargetPath)) {

		if (detvar) {

			// MC2023A production or detvars
			if ((entry.path().string().find("flat.caf.root") != std::string::npos) && (entry.path().extension() == ".root")) {

				std::string XRootPath = "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/" + entry.path().string().substr(6);
				if (print == true) std::cout << entry.path().string().substr(6) << std::endl;
				Input.push_back(XRootPath);

			}

		} else {

			//MC2024B production
			//2nd iterator bc we have subdirectories
			for (const auto & subentry : std::filesystem::directory_iterator(entry)) {

				//3rd iterator bc we have sub-subdirectories
				for (const auto & subsubentry : std::filesystem::directory_iterator(subentry)) {

					if ((subsubentry.path().string().find("caf.flat.caf") != std::string::npos) && (subsubentry.path().extension() == ".root")) {

						std::string XRootPath = "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/" + subsubentry.path().string().substr(6);
						if (print == true) std::cout << subsubentry.path().string().substr(6) << std::endl;
						Input.push_back(XRootPath);

					} // end of requirement that the naming scheme is satisfied

				} // end of loop over 3rd directory iterator

			}// end of loop over 2nd directory iterator

		} // detvar or regular file

	} // end of loop over 1st directory iterator

	return Input;

}

//---------------------//

#endif
