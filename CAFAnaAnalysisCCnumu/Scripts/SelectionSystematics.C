// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/ISyst.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"

// std includes.
#include <filesystem>
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionSystematics(std::string SystName, int SystNUniv, bool ModifiedResponse) {
    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Systematic with name " << SystName << ", and number of universes " << SystNUniv <<  std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;

    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir_figs+"/Figs/CAFAna/Uncertainties/"+SystName);

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematics"+TString(SystName)+".root";
    TFile* SaveFile = new TFile(RootFilePath, "UPDATE");
    
    // Create shift depending on number of universes
    ISyst* syst = new SBNWeightSyst(SystName);
    std::vector<SystShifts> Shifts; std::vector<Var> Weis; std::vector<TruthVar> TruthWeis;

    if (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) {
    	// Add +1 sigma shift
        SystShifts SigP1Shift(syst, +1);
	    Shifts.push_back(SigP1Shift);
    } else {
        Weis.reserve(SystNUniv);
        for (int i = 0; i < SystNUniv; i++) {
            Weis.push_back(GetUniverseWeight(SystName, i));
            TruthWeis.push_back(GetTruthUniverseWeight(SystName, i));
        }
    }

    // We now have the option to either load all the spectra from a previous run or 
    // run the spectra in this run
    const bool ConstructSpectra = true;

    // Where we store spectra if we are going to construct them    
    std::vector<std::tuple<
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>
    >> Spectra;

    // Where we store spectra needed for computing response matrices for universes
    // This will only be populated for cross section systematics
    std::vector<std::tuple<
        std::vector<std::unique_ptr<EnsembleSpectrum>>,
        std::unique_ptr<EnsembleSpectrum>
    >> ResponseMatrixSpectra;

    // Where we load histograms if we do not construct spectra
    std::vector<std::vector<std::tuple<TH1D*, TH1D*, TH1D*>>> LoadedHistos;

    // Again will only be loaded for cross section systematics
    std::vector<std::vector<std::tuple<
        std::vector<TH1D*>, // vector with RS for each bin
        TH1D*               // S_univ
    >>> ResponseMatrixHistos;

    // Finally also only loaded for cross section systematics
    std::vector<TH1D*> CVTrueSignalHistos;

    if (ConstructSpectra) {
        // Construct all spectra
        for (std::size_t iVar = 0; iVar < Vars.size(); iVar++) {
            std::unique_ptr<EnsembleSpectrum> RecoSpectra;
            std::unique_ptr<EnsembleSpectrum> RecoTrueSpectra;
            std::unique_ptr<EnsembleSpectrum> RecoBkgSpectra;

            if (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) {
                RecoSpectra = std::make_unique<EnsembleSpectrum>(
                    NuLoader, HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))),
                    kNoSpillCut, kRecoIsSignal, Shifts
                );
                RecoTrueSpectra = std::make_unique<EnsembleSpectrum>(
                    NuLoader, HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))),
                    kNoSpillCut, kRecoIsTrueReco, Shifts
                );
                RecoBkgSpectra = std::make_unique<EnsembleSpectrum>(
                    NuLoader, HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))),
                    kNoSpillCut, kRecoIsBackground, Shifts
                );
            } else {
                RecoSpectra = std::make_unique<EnsembleSpectrum>(
                    NuLoader, HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))),
                    kNoSpillCut, kRecoIsSignal, Weis
                );
                RecoTrueSpectra = std::make_unique<EnsembleSpectrum>(
                    NuLoader, HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))),
                    kNoSpillCut, kRecoIsTrueReco, Weis
                );
                RecoBkgSpectra = std::make_unique<EnsembleSpectrum>(
                    NuLoader, HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))),
                    kNoSpillCut, kRecoIsBackground, Weis
                );
            }
            Spectra.push_back({std::move(RecoSpectra), std::move(RecoTrueSpectra), std::move(RecoBkgSpectra)});

            // If we are dealing with cross section systematics, we have to use the modified response
            //  matrices in which we use the true signal spectrum for the universe in the denominator
            // instead of the central value true signal spectrum
            if (ModifiedResponse) {
                const std::vector<double>& BinEdges = VarBins.at(iVar).Edges();
                std::vector<std::unique_ptr<EnsembleSpectrum>> InnerSpectra;
                Var kCurrentVar = std::get<1>(Vars.at(iVar));

                // This will get us the histograms neccesary to compute the response matrix using 
                // the signal spectrum from the universe
                for (int j = 0; j < VarBins.at(iVar).NBins(); j++) {
                    double BinMin = BinEdges.at(j);
                    double BinMax = (j == VarBins.at(iVar).NBins() - 1) ?  VarBins.at(iVar).Max() : BinEdges.at(j + 1);

                    const Cut TempCut([=](const caf::SRSliceProxy* slc) {
                        return (
                            kRecoIsTrueReco(slc) && 
                            kCurrentVar(slc) >= BinMin &&
                            kCurrentVar(slc) < BinMax
                        );
                    });

                    std::unique_ptr<EnsembleSpectrum> RecoBinValues = nullptr;
                    std::unique_ptr<EnsembleSpectrum> TruthValues = nullptr;
                    if (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) {
                        RecoBinValues = std::make_unique<EnsembleSpectrum>(
                            NuLoader, HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))),
                            kNoSpillCut, TempCut, Shifts
                        );
                        
                    } else {
                        RecoBinValues = std::make_unique<EnsembleSpectrum>(
                            NuLoader, HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))),
                            kNoSpillCut, TempCut, Weis
                        );
                    }
                    InnerSpectra.push_back(std::move(RecoBinValues));
                }

                // Get true signal for universes
                std::unique_ptr<EnsembleSpectrum> TruthValues = nullptr;
                if (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) {
                    TruthValues = std::make_unique<EnsembleSpectrum>(
                        VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<2>(Vars.at(iVar)),
                        kTruthIsSignal, kNoSpillCut, Shifts
                    );
                } else {
                    TruthValues = std::make_unique<EnsembleSpectrum>(
                        VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<2>(Vars.at(iVar)),
                        kTruthIsSignal, kNoSpillCut, TruthWeis
                    );
                }
                ResponseMatrixSpectra.push_back({std::move(InnerSpectra), std::move(TruthValues)});
            }
        }
        // Load spectra
        NuLoader.Go();

    } else {
        // Load previously constructed histograms from file
        for (std::size_t i = 0; i < Vars.size(); i++) {
            std::vector<std::tuple<TH1D*, TH1D*, TH1D*>> VarHistos;
            std::vector<std::tuple<std::vector<TH1D*>, TH1D*>> ResponseHistos;

            // Nominal plots
            TH1D* RecoHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_reco"));
            TH1D* RecoTrueHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_reco_true"));
            TH1D* RecoBkgHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_reco_bkg"));

            RecoHisto->Scale(Units / (IntegratedFlux * NTargets));
            RecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
            RecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

            VarHistos.push_back({std::move(RecoHisto), std::move(RecoTrueHisto), std::move(RecoBkgHisto)});

            // Var univ plots
            int NUniv = (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) ? 1 : SystNUniv;
            for (int iUniv = 0; iUniv < NUniv; iUniv++) {
                TString UnivString = TString(std::to_string(iUniv));
                TH1D* UnivRecoHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco"));
                TH1D* UnivRecoTrueHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_true"));
                TH1D* UnivRecoBkgHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_bkg"));

                UnivRecoHisto->Scale(Units / (IntegratedFlux * NTargets));
                UnivRecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
                UnivRecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

                VarHistos.push_back({std::move(UnivRecoHisto), std::move(UnivRecoTrueHisto), std::move(UnivRecoBkgHisto)});

                std::vector<TH1D*> ResponseVarBinHistos; TH1D* UnivTruthHisto;
                if (ModifiedResponse) {
                    for (int j = 0; j < VarBins.at(i).NBins(); j++) {
                        TString BinString = TString(std::to_string(j));
                        TH1D* TruthValuesBinHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_"+BinString));
                        TruthValuesBinHisto->Scale(Units / (IntegratedFlux * NTargets));
                        ResponseVarBinHistos.push_back(std::move(TruthValuesBinHisto));
                    }
                    UnivTruthHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_true"));
                    UnivTruthHisto->Scale(Units / (IntegratedFlux * NTargets));
                }
                ResponseHistos.push_back({std::move(ResponseVarBinHistos), std::move(UnivTruthHisto)});
            }
            
            TH1D* CVTrueSignalHisto;
            if (ModifiedResponse) {
                CVTrueSignalHisto = (TH1D*)(SaveFile->Get<TH1D>(PlotNames[i]+"_true"));
                CVTrueSignalHisto->Scale(Units / (IntegratedFlux * NTargets));
            }
            
            // Push all histos for given variable
            LoadedHistos.push_back(std::move(VarHistos));
            ResponseMatrixHistos.push_back(std::move(ResponseHistos));
            CVTrueSignalHistos.push_back(std::move(CVTrueSignalHisto));

            // Resize Spectra so everything compiles correctly
            Spectra.resize(Vars.size());
        }
    }

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {
        // Get histograms
        TH1D* RecoHisto; TH1D* RecoTrueHisto; TH1D* RecoBkgHisto; TH1D* CVTrueSignalHisto;
        auto& [RecoSpectra, RecoTrueSpectra, RecoBkgSpectra] = Spectra.at(i);
        if (ConstructSpectra) {
            RecoHisto = RecoSpectra->Nominal().ToTH1(TargetPOT);
            RecoTrueHisto = RecoTrueSpectra->Nominal().ToTH1(TargetPOT);
            RecoBkgHisto = RecoBkgSpectra->Nominal().ToTH1(TargetPOT);

            // Manage under/overflow bins
            RecoHisto->SetBinContent(RecoHisto->GetNbinsX(), RecoHisto->GetBinContent(RecoHisto->GetNbinsX()) + RecoHisto->GetBinContent(RecoHisto->GetNbinsX() + 1));
            RecoTrueHisto->SetBinContent(RecoTrueHisto->GetNbinsX(), RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX()) + RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX() + 1));
            RecoBkgHisto->SetBinContent(RecoBkgHisto->GetNbinsX(), RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX()) + RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX() + 1));

            RecoHisto->SetBinContent(1, RecoHisto->GetBinContent(0) + RecoHisto->GetBinContent(1));
            RecoTrueHisto->SetBinContent(1, RecoTrueHisto->GetBinContent(0) + RecoTrueHisto->GetBinContent(1));
            RecoBkgHisto->SetBinContent(1, RecoBkgHisto->GetBinContent(0) + RecoBkgHisto->GetBinContent(1));

            // Save histos to root file
            SaveFile->WriteObject(RecoHisto, PlotNames[i]+"_reco");
            SaveFile->WriteObject(RecoTrueHisto, PlotNames[i]+"_reco_true");
            SaveFile->WriteObject(RecoBkgHisto, PlotNames[i]+"_reco_bkg");

            // Scale histograms
            RecoHisto->Scale(Units / (IntegratedFlux * NTargets));
            RecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
            RecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

            // Get central value true signal spectrum
            if (ModifiedResponse) {
                CVTrueSignalHisto = std::get<1>(ResponseMatrixSpectra[i])->Nominal().ToTH1(TargetPOT);
                CVTrueSignalHisto->SetBinContent(CVTrueSignalHisto->GetNbinsX(), CVTrueSignalHisto->GetBinContent(CVTrueSignalHisto->GetNbinsX()) + CVTrueSignalHisto->GetBinContent(CVTrueSignalHisto->GetNbinsX() + 1));
                SaveFile->WriteObject(CVTrueSignalHisto, PlotNames[i]+"_true");
                CVTrueSignalHisto->Scale(Units / (IntegratedFlux * NTargets));
            }
        } else {
            RecoHisto = std::get<0>(LoadedHistos.at(i)[0]);
            RecoTrueHisto = std::get<1>(LoadedHistos.at(i)[0]);
            RecoBkgHisto = std::get<2>(LoadedHistos.at(i)[0]);

            if (ModifiedResponse) {
                CVTrueSignalHisto = CVTrueSignalHistos.at(i);
            }
        }

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);

        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(RecoHisto,"Reconstructed","l");
        RecoHisto->SetLineColor(kBlue+2);
        RecoHisto->SetLineWidth(4);

        // Style histograms
        RecoHisto->GetXaxis()->SetTitleFont(FontStyle);
        RecoHisto->GetXaxis()->SetLabelFont(FontStyle);
        RecoHisto->GetXaxis()->SetNdivisions(8);
        RecoHisto->GetXaxis()->SetLabelSize(TextSize);
        RecoHisto->GetXaxis()->SetTitleSize(TextSize);
        RecoHisto->GetXaxis()->SetTitleOffset(1.1);
        RecoHisto->GetXaxis()->CenterTitle();
        RecoHisto->GetXaxis()->SetTitle(("Reco " + VarLabels.at(i)).c_str());

        RecoHisto->GetYaxis()->SetTitleFont(FontStyle);
        RecoHisto->GetYaxis()->SetLabelFont(FontStyle);
        RecoHisto->GetYaxis()->SetNdivisions(6);
        RecoHisto->GetYaxis()->SetLabelSize(TextSize);
        RecoHisto->GetYaxis()->SetTitleSize(TextSize);
        RecoHisto->GetYaxis()->SetTitleOffset(1.3);
        RecoHisto->GetYaxis()->SetTickSize(0);
        RecoHisto->GetYaxis()->CenterTitle();

        double imax = RecoHisto->GetMaximum();
        double YAxisRange = 1.3*imax;
        RecoHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoTrueHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoBkgHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);

        TLegendEntry* legRecoTrue = leg->AddEntry(RecoTrueHisto,"True","l");
        RecoTrueHisto->SetLineColor(kRed+1);
        RecoTrueHisto->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(RecoBkgHisto,"Background","l");
        RecoBkgHisto->SetLineColor(kOrange+7);
        RecoBkgHisto->SetLineWidth(4);

        // Get bins for matrices
        const int NBins = VarBins.at(i).NBins();
        const std::vector<double>& BinEdges = VarBins.at(i).Edges();

        // Create covariance matrix
        std::string CovName = "Cov" + SystName;
        TH2* CovMatrix = new TH2D(
            (CovName + (std::string)PlotNames[i]).c_str(),
            CovName.c_str(),
            NBins, BinEdges.data(),
            NBins, BinEdges.data()
        );

        // Create fractional covariance matrix
        std::string FracCovName = "FracCov" + SystName;
        TH2* FracCovMatrix = new TH2D(
            (FracCovName + (std::string)PlotNames[i]).c_str(),
            FracCovName.c_str(),
            NBins, BinEdges.data(),
            NBins, BinEdges.data()
        );

        // Create correlation matrix
        std::string CorrName = "Corr" + SystName;
        TH2* CorrMatrix = new TH2D(
            (CorrName + (std::string)PlotNames[i]).c_str(),
            CorrName.c_str(),
            NBins, BinEdges.data(),
            NBins, BinEdges.data()
        );

        // Loop over all universes
        int NUniv = (SystNUniv == 6 || SystNUniv == 10 || SystNUniv == 4 || SystNUniv == 2 || SystNUniv == 7) ? 1 : SystNUniv;
        for (int iUniv = 0; iUniv < NUniv; iUniv++) {
            // Get histograms
            TH1D* UnivRecoHisto; TH1D* UnivRecoTrueHisto; TH1D* UnivRecoBkgHisto;
            std::vector<TH1D*> UnivTruthValuesHistos; TH1D* UnivTrueSignalHisto;
            if (ConstructSpectra) {
                UnivRecoHisto = RecoSpectra->Universe(iUniv).ToTH1(TargetPOT);
                UnivRecoTrueHisto = RecoTrueSpectra->Universe(iUniv).ToTH1(TargetPOT);
                UnivRecoBkgHisto = RecoBkgSpectra->Universe(iUniv).ToTH1(TargetPOT);

                // Manage under/overflow bins
                UnivRecoHisto->SetBinContent(UnivRecoHisto->GetNbinsX(), UnivRecoHisto->GetBinContent(UnivRecoHisto->GetNbinsX()) + UnivRecoHisto->GetBinContent(UnivRecoHisto->GetNbinsX() + 1));
                UnivRecoTrueHisto->SetBinContent(UnivRecoTrueHisto->GetNbinsX(), UnivRecoTrueHisto->GetBinContent(UnivRecoTrueHisto->GetNbinsX()) + UnivRecoTrueHisto->GetBinContent(UnivRecoTrueHisto->GetNbinsX() + 1));
                UnivRecoBkgHisto->SetBinContent(UnivRecoBkgHisto->GetNbinsX(), UnivRecoBkgHisto->GetBinContent(UnivRecoBkgHisto->GetNbinsX()) + UnivRecoBkgHisto->GetBinContent(UnivRecoBkgHisto->GetNbinsX() + 1));

                UnivRecoHisto->SetBinContent(1, UnivRecoHisto->GetBinContent(0) + UnivRecoHisto->GetBinContent(1));
                UnivRecoTrueHisto->SetBinContent(1, UnivRecoTrueHisto->GetBinContent(0) + UnivRecoTrueHisto->GetBinContent(1));
                UnivRecoBkgHisto->SetBinContent(1, UnivRecoBkgHisto->GetBinContent(0) + UnivRecoBkgHisto->GetBinContent(1));

                // Save to root file
                TString UnivString = TString(std::to_string(iUniv));
                SaveFile->WriteObject(UnivRecoHisto, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco");
                SaveFile->WriteObject(UnivRecoTrueHisto, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_true");
                SaveFile->WriteObject(UnivRecoBkgHisto, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_bkg");

                // Scale histograms
                UnivRecoHisto->Scale(Units / (IntegratedFlux * NTargets));
                UnivRecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
                UnivRecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

                // Load histograms needed for response matrix if systematic is cross section
                if (ModifiedResponse) {
                    for (int iBin = 0; iBin < VarBins.at(i).NBins(); ++iBin) {
                        TH1D* TruthValuesBinHisto = std::get<0>(ResponseMatrixSpectra[i])[iBin]->Universe(iUniv).ToTH1(TargetPOT);
                        TruthValuesBinHisto->SetBinContent(TruthValuesBinHisto->GetNbinsX(), TruthValuesBinHisto->GetBinContent(TruthValuesBinHisto->GetNbinsX()) + TruthValuesBinHisto->GetBinContent(TruthValuesBinHisto->GetNbinsX() + 1));

                        // Save to root file
                        TString BinString = TString(std::to_string(iBin));
                        SaveFile->WriteObject(TruthValuesBinHisto, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_reco_"+BinString);
                        UnivTruthValuesHistos.push_back(TruthValuesBinHisto);

                        // Scale
                        TruthValuesBinHisto->Scale(Units / (IntegratedFlux * NTargets));
                    }
                    UnivTrueSignalHisto = std::get<1>(ResponseMatrixSpectra[i])->Universe(iUniv).ToTH1(TargetPOT);
                    UnivTrueSignalHisto->SetBinContent(UnivTrueSignalHisto->GetNbinsX(), UnivTrueSignalHisto->GetBinContent(UnivTrueSignalHisto->GetNbinsX()) + UnivTrueSignalHisto->GetBinContent(UnivTrueSignalHisto->GetNbinsX() + 1));
                    SaveFile->WriteObject(UnivTrueSignalHisto, PlotNames[i]+"_"+(TString)SystName+"_"+UnivString+"_true");
                    UnivTrueSignalHisto->Scale(Units / (IntegratedFlux * NTargets));
                } 
            } else {
                UnivRecoHisto = std::get<0>(LoadedHistos.at(i)[iUniv + 1]);
                UnivRecoTrueHisto = std::get<1>(LoadedHistos.at(i)[iUniv + 1]);
                UnivRecoBkgHisto = std::get<2>(LoadedHistos.at(i)[iUniv + 1]);

                if (ModifiedResponse) {
                    for (int iBin = 0; iBin < VarBins.at(i).NBins(); ++iBin) {
                        UnivTruthValuesHistos.push_back(std::get<0>(ResponseMatrixHistos[i][iUniv])[iBin]);
                    }
                    UnivTrueSignalHisto = std::get<1>(ResponseMatrixHistos[i][iUniv]);
                }
            }

            for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
                double XEventRateCV = RecoHisto->GetBinContent(x);
                double XEventRateVar = UnivRecoHisto->GetBinContent(x);
                for (int y = 1; y < VarBins.at(i).NBins() + 1; y++) {
                    double YEventRateCV = RecoHisto->GetBinContent(y);
                    double YEventRateVar = UnivRecoHisto->GetBinContent(y);

                    if (ModifiedResponse) {
                        XEventRateVar = UnivRecoBkgHisto->GetBinContent(x);
                        YEventRateVar = UnivRecoBkgHisto->GetBinContent(y);

                        for (int iBin = 0; iBin < VarBins.at(i).NBins(); ++iBin) {
                            double factor = CVTrueSignalHisto->GetBinContent(iBin + 1) / UnivTrueSignalHisto->GetBinContent(iBin + 1);
                            XEventRateVar += UnivTruthValuesHistos[iBin]->GetBinContent(x) * factor;
                            YEventRateVar += UnivTruthValuesHistos[iBin]->GetBinContent(y) * factor;
                        }
                    }
                    double Value = ((XEventRateVar - XEventRateCV) * (YEventRateVar - YEventRateCV)) / NUniv;

                    // Fill covariance matrix
                    if (TMath::Abs(Value) < 1e-14) Value = 1e-14;
                    CovMatrix->Fill(
                        RecoHisto->GetXaxis()->GetBinCenter(x),
                        RecoHisto->GetXaxis()->GetBinCenter(y),
                        Value
                    );
                }
	        }
        }

        // Create fractional covariance and correlation matrices
        for (int x = 1; x < VarBins.at(i).NBins() + 1; x++) {
            double XEventRateCV = RecoHisto->GetBinContent(x);
            for (int y = 1; y < VarBins.at(i).NBins() + 1; y++) {
                double YEventRateCV = RecoHisto->GetBinContent(y);
                double CovBinValue = CovMatrix->GetBinContent(x,y);
                double XBinValue = CovMatrix->GetBinContent(x,x);
                double YBinValue = CovMatrix->GetBinContent(y,y);

                // Fill frac cov matrix
                double FracValue = (XEventRateCV == 0. || YEventRateCV == 0.) ? 0. : CovBinValue / (XEventRateCV * YEventRateCV);
                if (TMath::Abs(FracValue) < 1e-14) FracValue = 1e-14;
                FracCovMatrix->SetBinContent(x, y, FracValue);

                // Fill corr matrix
                double CorrValue = (XBinValue == 0. || YBinValue == 0.) ? 0. : CovBinValue / (TMath::Sqrt(XBinValue) * TMath::Sqrt(YBinValue));
                if (TMath::Abs(CorrValue) < 1e-14) CorrValue = 1e-14;
                CorrMatrix->SetBinContent(x, y, CorrValue);
            }
        }
            
        // Plot cov matrix
        double CovMin = CovMatrix->GetMinimum();
        double CovMax = CovMatrix->GetMaximum();
        CovMatrix->GetZaxis()->SetRangeUser(CovMin,CovMax);
        CovMatrix->GetZaxis()->CenterTitle();
        CovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        CovMatrix->GetZaxis()->SetTitleSize(TextSize);
        CovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        CovMatrix->GetZaxis()->SetLabelSize(TextSize);
        CovMatrix->GetZaxis()->SetNdivisions(5);

        CovMatrix->GetXaxis()->SetTitle(("bin i " + VarLabels.at(i)).c_str());
        CovMatrix->GetXaxis()->CenterTitle();
        CovMatrix->GetXaxis()->SetTitleOffset(1.1);
        CovMatrix->GetXaxis()->SetTitleFont(FontStyle);
        CovMatrix->GetXaxis()->SetTitleSize(TextSize);
        CovMatrix->GetXaxis()->SetLabelFont(FontStyle);
        CovMatrix->GetXaxis()->SetLabelSize(TextSize);
        CovMatrix->GetXaxis()->SetNdivisions(5);

        CovMatrix->GetYaxis()->SetTitle(("bin j " + VarLabels.at(i)).c_str());
        CovMatrix->GetYaxis()->CenterTitle();
        CovMatrix->GetYaxis()->SetTitleOffset(1.1);
        CovMatrix->GetYaxis()->SetTitleFont(FontStyle);
        CovMatrix->GetYaxis()->SetTitleSize(TextSize);
        CovMatrix->GetYaxis()->SetLabelFont(FontStyle);
        CovMatrix->GetYaxis()->SetLabelSize(TextSize);
        CovMatrix->GetYaxis()->SetNdivisions(5);

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        CovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/Cov"+PlotNames[i]+".pdf");

        // Plot frac cov matrix
        double FracCovMin = FracCovMatrix->GetMinimum();
        double FracCovMax = FracCovMatrix->GetMaximum();
        FracCovMatrix->GetZaxis()->SetRangeUser(FracCovMin,FracCovMax); // set the ranges accordingly, for frac cov should be [0,100], for corr matrices [-1,1]
        FracCovMatrix->GetZaxis()->CenterTitle();
        FracCovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        FracCovMatrix->GetZaxis()->SetTitleSize(TextSize);
        FracCovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        FracCovMatrix->GetZaxis()->SetLabelSize(TextSize);
        FracCovMatrix->GetZaxis()->SetNdivisions(5);

        FracCovMatrix->GetXaxis()->SetTitle(("bin i " + VarLabels.at(i)).c_str());
        FracCovMatrix->GetXaxis()->CenterTitle();
        FracCovMatrix->GetXaxis()->SetTitleOffset(1.1);
        FracCovMatrix->GetXaxis()->SetTitleFont(FontStyle);
        FracCovMatrix->GetXaxis()->SetTitleSize(TextSize);
        FracCovMatrix->GetXaxis()->SetLabelFont(FontStyle);
        FracCovMatrix->GetXaxis()->SetLabelSize(TextSize);
        FracCovMatrix->GetXaxis()->SetNdivisions(5);

        FracCovMatrix->GetYaxis()->SetTitle(("bin j " + VarLabels.at(i)).c_str());
        FracCovMatrix->GetYaxis()->CenterTitle();
        FracCovMatrix->GetYaxis()->SetTitleOffset(1.1);
        FracCovMatrix->GetYaxis()->SetTitleFont(FontStyle);
        FracCovMatrix->GetYaxis()->SetTitleSize(TextSize);
        FracCovMatrix->GetYaxis()->SetLabelFont(FontStyle);
        FracCovMatrix->GetYaxis()->SetLabelSize(TextSize);
        FracCovMatrix->GetYaxis()->SetNdivisions(5);

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        FracCovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/FracCov"+PlotNames[i]+".pdf");

        // Plot correlation matrix
        CorrMatrix->GetZaxis()->SetRangeUser(-1,1);
        CorrMatrix->GetZaxis()->CenterTitle();
        CorrMatrix->GetZaxis()->SetTitleFont(FontStyle);
        CorrMatrix->GetZaxis()->SetTitleSize(TextSize);
        CorrMatrix->GetZaxis()->SetLabelFont(FontStyle);
        CorrMatrix->GetZaxis()->SetLabelSize(TextSize);
        CorrMatrix->GetZaxis()->SetNdivisions(5);

        CorrMatrix->GetXaxis()->SetTitle(("bin i " + VarLabels.at(i)).c_str());
        CorrMatrix->GetXaxis()->CenterTitle();
        CorrMatrix->GetXaxis()->SetTitleOffset(1.1);
        CorrMatrix->GetXaxis()->SetTitleFont(FontStyle);
        CorrMatrix->GetXaxis()->SetTitleSize(TextSize);
        CorrMatrix->GetXaxis()->SetLabelFont(FontStyle);
        CorrMatrix->GetXaxis()->SetLabelSize(TextSize);
        CorrMatrix->GetXaxis()->SetNdivisions(5);

        CorrMatrix->GetYaxis()->SetTitle(("bin j " + VarLabels.at(i)).c_str());
        CorrMatrix->GetYaxis()->CenterTitle();
        CorrMatrix->GetYaxis()->SetTitleOffset(1.1);
        CorrMatrix->GetYaxis()->SetTitleFont(FontStyle);
        CorrMatrix->GetYaxis()->SetTitleSize(TextSize);
        CorrMatrix->GetYaxis()->SetLabelFont(FontStyle);
        CorrMatrix->GetYaxis()->SetLabelSize(TextSize);
        CorrMatrix->GetYaxis()->SetNdivisions(5);

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        CorrMatrix->Draw("colz text");
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/Corr"+PlotNames[i]+".pdf");

        // Save objects
        SaveFile->WriteObject(CovMatrix, PlotNames[i]+"_cov");
        SaveFile->WriteObject(FracCovMatrix, PlotNames[i]+"_fraccov");
        SaveFile->WriteObject(CorrMatrix, PlotNames[i]+"_corr");

        if (ConstructSpectra) {
            // Plot histograms with error bands, only when constructing spectra
            TGraphAsymmErrors* RecoErrorBand = RecoSpectra->ErrorBand(TargetPOT);
            TGraphAsymmErrors* RecoTrueErrorBand = RecoTrueSpectra->ErrorBand(TargetPOT);
            TGraphAsymmErrors* RecoBkgErrorBand = RecoBkgSpectra->ErrorBand(TargetPOT);

            PlotCanvas->cd();

            PlotCanvas->SetTopMargin(0.13);
            PlotCanvas->SetLeftMargin(0.17);
            PlotCanvas->SetRightMargin(0.05);
            PlotCanvas->SetBottomMargin(0.16);

            // Undo scale
            RecoHisto->Scale((IntegratedFlux * NTargets) / Units);
            RecoTrueHisto->Scale((IntegratedFlux * NTargets) / Units);
            RecoBkgHisto->Scale((IntegratedFlux * NTargets) / Units);

            double imax = RecoHisto->GetMaximum();
            double YAxisRange = 1.3*imax;
            RecoHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
            RecoTrueHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
            RecoBkgHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);

            RecoHisto->Draw("hist");
            ana::DrawErrorBand(RecoHisto, RecoErrorBand);
            RecoTrueHisto->Draw("hist same");
            ana::DrawErrorBand(RecoTrueHisto, RecoTrueErrorBand);
            RecoBkgHisto->Draw("hist same");
            ana::DrawErrorBand(RecoBkgHisto, RecoBkgErrorBand);
            leg->Draw();
            PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/Uncertainties/"+(TString)SystName+"/"+PlotNames[i]+".pdf");
        }
        delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();
}
