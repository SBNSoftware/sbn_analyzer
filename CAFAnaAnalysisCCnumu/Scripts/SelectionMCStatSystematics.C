// SBNAna includes.
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/HistAxis.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"

// std includes.
#include <filesystem>
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"
#include "Helpers.cpp"

// Utils includes.
#include "../../Utils/Constants.h"
#include "../../Utils/Tools.cxx"

using namespace std;
using namespace ana;
using namespace Constants;

void SelectionMCStatSystematics() {
    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    int FontStyle = 132;
    double TextSize = 0.06;

    // Initialize tools
    Tools tools;

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);

    // Create directory for this sytematic if it does not exist yet
    std::filesystem::create_directory((std::string)dir_figs+"/Figs/CAFAna/Uncertainties/MCStat");

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/SelectionSystematicsMCStat.root";
    TFile* SaveFile = new TFile(RootFilePath, "recreate");

    const int NVars = Vars.size();

    // Create ISyst that computes weight for each event
    class MCStatSyst : public ISyst {
        public:
            MCStatSyst(int UnivIndex) : ISyst("MCStat"+std::to_string(UnivIndex), "MCStat"+std::to_string(UnivIndex)), UnivIndex(UnivIndex) {}

            // Use fmatch time to create random seed
            // This might need to be revisited in the future
            void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const {
                int seed = 1;
                if (slc->fmatch.time > 0) {
                    seed = slc->fmatch.time * 10000 + UnivIndex;
                } 
                // in the case time is negative, the event is
                // not going to be used so we don't care about the seed
                weight *= tools.PoissonRandomNumber(seed);
            }

            // does not get called in our systematics, so we don't really care about the seed
            void Shift(double sigma, caf::SRTrueInteractionProxy* nu, double& weight) const {
                weight *= tools.PoissonRandomNumber(1);
            }

            const Tools tools;
            int UnivIndex;
    };

    std::vector<SystShifts> Shifts; 
    for (int i = 0; i < 100; ++i) {
        ISyst* Syst = new MCStatSyst(i);
        SystShifts Shift(Syst, +1);
        Shifts.push_back(Shift);
    }

    // Construct all spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>,
        std::unique_ptr<EnsembleSpectrum>
    >> Spectra;
    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsSignal); 
        auto RecoTrueSignals = std::make_unique<Spectrum> (VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsTrueReco); 
        auto RecoBkgSignals = std::make_unique<Spectrum>(VarLabels.at(iVar), VarBins.at(iVar), NuLoader, std::get<0>(Vars.at(iVar)), kNoSpillCut, kRecoIsBackground); 

        auto UnivRecoSignals = std::make_unique<EnsembleSpectrum>(
            NuLoader,
            HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))), 
            kNoSpillCut,
            kRecoIsSignal, 
            Shifts
        );
        auto UnivRecoTrueSignals = std::make_unique<EnsembleSpectrum>(
            NuLoader,
            HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))), 
            kNoSpillCut,
            kRecoIsTrueReco, 
            Shifts
        ); 
        auto UnivRecoBkgSignals = std::make_unique<EnsembleSpectrum>(
            NuLoader,
            HistAxis(VarLabels.at(iVar), VarBins.at(iVar), std::get<0>(Vars.at(iVar))), 
            kNoSpillCut,
            kRecoIsBackground, 
            Shifts
        ); 
        Spectra.push_back({
            std::move(RecoSignals), 
            std::move(RecoTrueSignals), 
            std::move(RecoBkgSignals),
            std::move(UnivRecoSignals), 
            std::move(UnivRecoTrueSignals), 
            std::move(UnivRecoBkgSignals)
        });
    }
    NuLoader.Go();

    for (int iVar = 0; iVar < NVars; ++iVar) {
        auto& [
            RecoSpectra, 
            RecoTrueSpectra, 
            RecoBkgSpectra,
            UnivRecoSpectra, 
            UnivRecoTrueSpectra, 
            UnivRecoBkgSpectra
        ] = Spectra.at(iVar);

        // Get nominal plots
        TH1D* RecoHisto = RecoSpectra->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSpectra->ToTH1(TargetPOT);
        TH1D* RecoBkgHisto = RecoBkgSpectra->ToTH1(TargetPOT);

        // Get universe plots
        std::vector<TH1D*> UnivRecoHistos;
        std::vector<TH1D*> UnivRecoTrueHistos;
        std::vector<TH1D*> UnivRecoBkgHistos;
        for (int i = 0; i < 100; ++i) {
            TH1D* ModifiedRecoHisto = UnivRecoSpectra->Universe(i).ToTH1(TargetPOT);
            UnivRecoHistos.push_back(std::move(ModifiedRecoHisto));
            TH1D* ModifiedRecoTrueHisto = UnivRecoTrueSpectra->Universe(i).ToTH1(TargetPOT);
            UnivRecoTrueHistos.push_back(std::move(ModifiedRecoTrueHisto));
            TH1D* ModifiedRecoBkgHisto = UnivRecoBkgSpectra->Universe(i).ToTH1(TargetPOT);
            UnivRecoBkgHistos.push_back(std::move(ModifiedRecoBkgHisto));
        }

        // Plot histograms with error band
        SelectionHelpers::DrawHistosWithErrorBands(
            RecoHisto,
            RecoTrueHisto,
            RecoBkgHisto,
            UnivRecoHistos,
            UnivRecoTrueHistos,
            UnivRecoBkgHistos,
            dir_figs, 
            "MCStat",
            PlotNames[iVar]
        );

        // Scale histograms for cov matrices
        RecoHisto->Scale(Units / (IntegratedFlux * NTargets));
        RecoTrueHisto->Scale(Units / (IntegratedFlux * NTargets));
        RecoBkgHisto->Scale(Units / (IntegratedFlux * NTargets));

        // Scale modified histograms for cov matrices
        for (int i = 0; i < 100; ++i) {
            UnivRecoHistos[i]->Scale(Units / (IntegratedFlux * NTargets));
            UnivRecoTrueHistos[i]->Scale(Units / (IntegratedFlux * NTargets));
            UnivRecoBkgHistos[i]->Scale(Units / (IntegratedFlux * NTargets));
        }

        // Plot cov, frac cov, and corr matrices
        SelectionHelpers::DrawMatrices(
            RecoHisto,
            RecoTrueHisto,
            RecoBkgHisto,
            UnivRecoHistos,
            UnivRecoTrueHistos,
            UnivRecoBkgHistos,
            dir_figs,
            "MCStat",
            VarLabels[iVar],
            PlotNames[iVar],
            SaveFile
        );
    }
    // Close file
    SaveFile->Close();
}
