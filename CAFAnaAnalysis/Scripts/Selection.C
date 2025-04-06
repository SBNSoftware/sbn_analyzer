// SBNAna includes.
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

// ROOT includes.
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TH1D.h"
#include "TPad.h"

// std includes.
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void Selection() {

    //----------------------------------------//

    // Set defaults and load tools
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();

    //----------------------------------------//
    
    // MC spectrum loader

    // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
    SpectrumLoader NuLoader(InputFiles);
    SpectrumLoader DataNuLoader(DataInputFiles);

    // We now create overlaid plots for several reconstructed variables and three lines:
    //     1. all selected reconstructed events
    //     2. reco signal events
    //     3. reco background events

    // Root file to store objects in
    TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Selection.root";
    TFile* SaveFile = new TFile(RootFilePath, "RECREATE");

    // Construct all spectra
    std::vector<std::tuple<
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>,
        std::unique_ptr<Spectrum>
    >> Spectra;

    std::vector< std::unique_ptr<Spectrum> > data_Spectra;    

    for (std::size_t i = 0; i < Vars.size(); i++) {

        auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsSignal); 
        auto RecoTrueSignals = std::make_unique<Spectrum> (VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsTrueReco); 
        auto RecoBkgSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsBackground); 
        Spectra.push_back({std::move(RecoSignals), std::move(RecoTrueSignals), std::move(RecoBkgSignals)});

        auto data_RecoSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), DataNuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsSignal);
        data_Spectra.push_back( std::move(data_RecoSignals) );        

    }

    // We now create the single-bin spectra that will help us get the efficiency and purity data for each of the cuts

    // MC

    // Spectrum with all events
    Spectrum sAllEvents("AllEvents", bEventCount, NuLoader, kTrueEventCount, kNoTruthCut, kNoSpillCut);
    // Spectrum with all events that were reconstructed
    Spectrum sAllRecoEvents("AllRecoEvents", bEventCount, NuLoader, kTrueEventCount, kNoTruthCut, kNoSpillCut, kNoCut);
    // Spectrum with all true signal events
    Spectrum sAllTrueEvents("AllTrueEvents", bEventCount, NuLoader, kTrueEventCount, kTruthIsSignal, kNoSpillCut);
    // Spectrum with all true signal events that were reconstructed
    Spectrum sAllTrueRecoEvents("AllTrueRecoEvents", bEventCount, NuLoader, kTrueEventCount, kTruthIsSignal, kNoSpillCut, kNoCut);
    // Spectrum with first cut (cosmic)
    Spectrum sFirstCut("FirstCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFirstCut);
    Spectrum sFirstCutTrue("FirstCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFirstCutTrue);
    // Spectrum with second cut (cosmic and vertex FV)
    Spectrum sSecondCut("SecondCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kSecondCut);
    Spectrum sSecondCutTrue("SecondCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kSecondCutTrue);
    // Spectrum with second cut (cosmic, vertex FV, and one muon)
    Spectrum sThirdCut("ThirdCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kThirdCut);
    Spectrum sThirdCutTrue("ThirdCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kThirdCutTrue);
    // Spectrum with second cut (cosmic, vertex FV, one muon, and two protons)
    Spectrum sFourthCut("FourthCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFourthCut);
    Spectrum sFourthCutTrue("FourthCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFourthCutTrue);
    // Spectrum with second cut (cosmic, vertex FV, one muon, two protons, and no charged pions)
    Spectrum sFifthCut("FifthCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFifthCut);
    Spectrum sFifthCutTrue("FifthCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kFifthCutTrue);
    // Spectrum with second cut (cosmic, vertex FV, one muon, two protons, no charged pions, and no neutral pions)
    Spectrum sSixthCut("SixthCut", bEventCount, NuLoader, kEventCount, kNoSpillCut, kSixthCut);
    Spectrum sSixthCutTrue("SixthCutTrue", bEventCount, NuLoader, kEventCount, kNoSpillCut, kSixthCutTrue);
    // Spectrum with overall signal definition to sanity check it matches
    Spectrum sRecoSignal("RecoSignal", bEventCount, NuLoader, kEventCount, kNoSpillCut, kRecoIsSignal); 
    Spectrum sRecoTrueSignal("RecoTrueSignal", bEventCount, NuLoader, kEventCount, kNoSpillCut, kRecoIsTrueReco); 

    // Data

    // Spectrum with all events that were reconstructed
    Spectrum sDataAllRecoEvents("DataAllRecoEvents", bEventCount, DataNuLoader, kEventCount, kNoSpillCut, kNoCut);
    // Spectrum with first cut (cosmic)
    Spectrum sDataFirstCut("DataFirstCut", bEventCount, DataNuLoader, kEventCount, kNoSpillCut, kFirstCut);
    // Spectrum with second cut (cosmic and vertex FV)
    Spectrum sDataSecondCut("DataSecondCut", bEventCount, DataNuLoader, kEventCount, kNoSpillCut, kSecondCut);
    // Spectrum with second cut (cosmic, vertex FV, and one muon)
    Spectrum sDataThirdCut("DataThirdCut", bEventCount, DataNuLoader, kEventCount, kNoSpillCut, kThirdCut);
    // Spectrum with second cut (cosmic, vertex FV, one muon, and two protons)
    Spectrum sDataFourthCut("DataFourthCut", bEventCount, DataNuLoader, kEventCount, kNoSpillCut, kFourthCut);
    // Spectrum with second cut (cosmic, vertex FV, one muon, two protons, and no charged pions)
    Spectrum sDataFifthCut("DataFifthCut", bEventCount, DataNuLoader, kEventCount, kNoSpillCut, kFifthCut);
    // Spectrum with second cut (cosmic, vertex FV, one muon, two protons, no charged pions, and no neutral pions)
    Spectrum sDataSixthCut("DataSixthCut", bEventCount, DataNuLoader, kEventCount, kNoSpillCut, kSixthCut);
    // Spectrum with overall signal definition
    Spectrum sDataRecoSignal("DataRecoSignal", bEventCount, DataNuLoader, kEventCount, kNoSpillCut, kRecoIsSignal); 

    DataNuLoader.Go();
    NuLoader.Go();

    //----------------------------------------//    

    // Loop over variables
    for (std::size_t i = 0; i < Vars.size(); i++) {

        auto& [RecoSignals, RecoTrueSignals, RecoBkgSignals] = Spectra.at(i);
        auto& data_RecoSignals = data_Spectra.at(i);        

        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.15);        

        TH1D* RecoHisto = RecoSignals->ToTH1(TargetPOT);
        TH1D* RecoTrueHisto = RecoTrueSignals->ToTH1(TargetPOT);
        TH1D* RecoBkgHisto = RecoBkgSignals->ToTH1(TargetPOT);   
        TH1D* data_RecoHisto = data_RecoSignals->ToTH1(TargetPOT);

        //overflow bins
        RecoHisto->SetBinContent(RecoHisto->GetNbinsX(), RecoHisto->GetBinContent(RecoHisto->GetNbinsX()) + RecoHisto->GetBinContent(RecoHisto->GetNbinsX() + 1));
        RecoTrueHisto->SetBinContent(RecoTrueHisto->GetNbinsX(), RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX()) + RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX() + 1));
        RecoBkgHisto->SetBinContent(RecoBkgHisto->GetNbinsX(), RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX()) + RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX() + 1));
        data_RecoHisto->SetBinContent(data_RecoHisto->GetNbinsX(), data_RecoHisto->GetBinContent(data_RecoHisto->GetNbinsX()) + data_RecoHisto->GetBinContent(data_RecoHisto->GetNbinsX() + 1));

        // undeflow bins
        RecoHisto->SetBinContent(1, RecoHisto->GetBinContent(0) + RecoHisto->GetBinContent(1));
        RecoTrueHisto->SetBinContent(1, RecoTrueHisto->GetBinContent(0) + RecoTrueHisto->GetBinContent(1));
        RecoBkgHisto->SetBinContent(1, RecoBkgHisto->GetBinContent(0) + RecoBkgHisto->GetBinContent(1));
        data_RecoHisto->SetBinContent(1, data_RecoHisto->GetBinContent(0) + data_RecoHisto->GetBinContent(1));

        TLegend* leg = new TLegend(0.15,0.9,0.7,0.98);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(RecoHisto,"MC","l");
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
        RecoHisto->GetXaxis()->SetTitle((VarLabels.at(i)).c_str());
        RecoHisto->GetXaxis()->SetNdivisions(6);        

        RecoHisto->GetYaxis()->SetTitleFont(FontStyle);
        RecoHisto->GetYaxis()->SetLabelFont(FontStyle);
        RecoHisto->GetYaxis()->SetNdivisions(6);
        RecoHisto->GetYaxis()->SetLabelSize(TextSize);
        RecoHisto->GetYaxis()->SetTitleSize(TextSize);
        RecoHisto->GetYaxis()->SetTitleOffset(1.3);
        RecoHisto->GetYaxis()->CenterTitle();
        RecoHisto->GetYaxis()->SetNdivisions(6);             

        TLegendEntry* legRecoData = leg->AddEntry(data_RecoHisto,"Data","ep");
        TLegendEntry* legRecoTrue = leg->AddEntry(RecoTrueHisto,"CC2p0#pi","l");
        RecoTrueHisto->SetLineColor(kAzure-4);
        RecoTrueHisto->SetFillColor(kAzure-4);        
        RecoTrueHisto->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(RecoBkgHisto,"non-CC2p0#pi","l");
        RecoBkgHisto->SetLineColor(kOrange+7);
        RecoBkgHisto->SetFillColor(kOrange+7);        
        RecoBkgHisto->SetLineWidth(4);

        data_RecoHisto->SetLineColor(kBlack);
        data_RecoHisto->SetFillColor(kBlack);        
        data_RecoHisto->SetLineWidth(2);   
        data_RecoHisto->SetMarkerStyle(8);
        data_RecoHisto->SetMarkerSize(2.);                        

        PlotCanvas->cd();
        TH1D* total_mc = (TH1D*)(RecoTrueHisto->Clone());
        total_mc->Add(RecoBkgHisto);

        //------------------------------------//

        // Area normalized
        double mc_sf = 1./total_mc->Integral("width");
        total_mc->Scale(mc_sf);
        RecoHisto->Scale(mc_sf);
        RecoTrueHisto->Scale(mc_sf);
        RecoBkgHisto->Scale(mc_sf);        

        double data_sf = 1./data_RecoHisto->Integral("width");
        data_RecoHisto->Scale(data_sf);

        //------------------------------------// 

        double imax = RecoHisto->GetMaximum();
        double YAxisRange = 1.3*imax;
        RecoHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoTrueHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoBkgHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        data_RecoHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);           

        //------------------------------------//        

        RecoHisto->Draw("hist same");        
        //RecoTrueHisto->Draw("hist same");       
        total_mc->Draw("hist same");
        RecoBkgHisto->Draw("hist same");
        //redraw so that it is on top
        RecoHisto->Draw("hist same");   
        // finally the data
        data_RecoHisto->Draw("same e1x0");
        
        leg->Draw();

        gPad->RedrawAxis();

        // Save as pdf
        PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/"+PlotNames[i]+".pdf");

        // Save to root file
        SaveFile->WriteObject(RecoHisto, PlotNames[i]+"_reco");
        SaveFile->WriteObject(RecoTrueHisto, PlotNames[i]+"_reco_true");
        SaveFile->WriteObject(RecoBkgHisto, PlotNames[i]+"_bkg");

        delete PlotCanvas;
    }
    // Close file
    SaveFile->Close();

    // Get histograms for all cuts

    // MC

    TH1D* AllEventsHisto = sAllEvents.ToTH1(TargetPOT);
    TH1D* AllRecoEventsHisto = sAllRecoEvents.ToTH1(TargetPOT);
    TH1D* AllTrueEventsHisto = sAllTrueEvents.ToTH1(TargetPOT);
    TH1D* AllTrueRecoEventsHisto = sAllTrueRecoEvents.ToTH1(TargetPOT);

    TH1D* FirstCutHisto = sFirstCut.ToTH1(TargetPOT);
    TH1D* FirstCutTrueHisto = sFirstCutTrue.ToTH1(TargetPOT);

    TH1D* SecondCutHisto = sSecondCut.ToTH1(TargetPOT);
    TH1D* SecondCutTrueHisto = sSecondCutTrue.ToTH1(TargetPOT);

    TH1D* ThirdCutHisto = sThirdCut.ToTH1(TargetPOT);
    TH1D* ThirdCutTrueHisto = sThirdCutTrue.ToTH1(TargetPOT);

    TH1D* FourthCutHisto = sFourthCut.ToTH1(TargetPOT);
    TH1D* FourthCutTrueHisto = sFourthCutTrue.ToTH1(TargetPOT);

    TH1D* FifthCutHisto = sFifthCut.ToTH1(TargetPOT);
    TH1D* FifthCutTrueHisto = sFifthCutTrue.ToTH1(TargetPOT);

    TH1D* SixthCutHisto = sSixthCut.ToTH1(TargetPOT);
    TH1D* SixthCutTrueHisto = sSixthCutTrue.ToTH1(TargetPOT);

    TH1D* RecoSignalHisto = sRecoSignal.ToTH1(TargetPOT);
    TH1D* RecoTrueSignalHisto = sRecoTrueSignal.ToTH1(TargetPOT);

    // Data

    TH1D* DataAllRecoEventsHisto = sDataAllRecoEvents.ToTH1(TargetPOT);
    TH1D* DataFirstCutHisto = sDataFirstCut.ToTH1(TargetPOT);
    TH1D* DataSecondCutHisto = sDataSecondCut.ToTH1(TargetPOT);
    TH1D* DataThirdCutHisto = sDataThirdCut.ToTH1(TargetPOT);
    TH1D* DataFourthCutHisto = sDataFourthCut.ToTH1(TargetPOT);
    TH1D* DataFifthCutHisto = sDataFifthCut.ToTH1(TargetPOT);
    TH1D* DataSixthCutHisto = sDataSixthCut.ToTH1(TargetPOT);
    TH1D* DataRecoSignalHisto = sDataRecoSignal.ToTH1(TargetPOT);

    // Get integrals for all cuts for efficiency / purity study

    // MC

    double AllEventsInt = AllEventsHisto->Integral("width");
    double AllRecoEventsInt = AllRecoEventsHisto->Integral("width");
    double AllTrueEventsInt = AllTrueEventsHisto->Integral("width");
    double AllTrueRecoEventsInt = AllTrueRecoEventsHisto->Integral("width");

    double FirstCutInt = FirstCutHisto->Integral("width");
    double FirstCutTrueInt = FirstCutTrueHisto->Integral("width");

    double SecondCutInt = SecondCutHisto->Integral("width");
    double SecondCutTrueInt = SecondCutTrueHisto->Integral("width");

    double ThirdCutInt = ThirdCutHisto->Integral("width");
    double ThirdCutTrueInt = ThirdCutTrueHisto->Integral("width");

    double FourthCutInt = FourthCutHisto->Integral("width");
    double FourthCutTrueInt = FourthCutTrueHisto->Integral("width");

    double FifthCutInt = FifthCutHisto->Integral("width");
    double FifthCutTrueInt = FifthCutTrueHisto->Integral("width");

    double SixthCutInt = SixthCutHisto->Integral("width");
    double SixthCutTrueInt = SixthCutTrueHisto->Integral("width");

    double RecoSignalInt = RecoSignalHisto->Integral("width");
    double RecoTrueSignalInt = RecoTrueSignalHisto->Integral("width");

    // Data

    double DataAllRecoEventsInt = DataAllRecoEventsHisto->Integral("width");
    double DataFirstCutInt = DataFirstCutHisto->Integral("width");
    double DataSecondCutInt = DataSecondCutHisto->Integral("width");
    double DataThirdCutInt = DataThirdCutHisto->Integral("width");
    double DataFourthCutInt = DataFourthCutHisto->Integral("width");
    double DataFifthCutInt = DataFifthCutHisto->Integral("width");
    double DataSixthCutInt = DataSixthCutHisto->Integral("width");
    double DataRecoSignalInt = DataRecoSignalHisto->Integral("width");    

    // Print results

    // MC

    std::cout << std::endl;
    std::cout << "============== MC ==================" << std::endl;
    std::cout << "All events: " << AllEventsInt << std::endl;
    std::cout << "Reconstructed events: " << AllRecoEventsInt << std::endl;
    std::cout << "True signal events: " << AllTrueEventsInt << std::endl;
    std::cout << "True signal events that were reconstructed: " << AllTrueRecoEventsInt << std::endl;
    std::cout << std::endl;
    std::cout << "Cuts: " << std::endl;
    std::cout << "    Cosmic cut: " << FirstCutInt << ". G.E: " <<  (FirstCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (FirstCutTrueInt / AllTrueEventsInt) * 100. << ". Purity: " << (FirstCutTrueInt / FirstCutInt) * 100. << std::endl;
    std::cout << "    Vertex in FV cut: " << SecondCutInt << ". G.E: " <<  (SecondCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (SecondCutTrueInt / AllTrueEventsInt) * 100. << ". Purity: " << (SecondCutTrueInt / SecondCutInt) * 100. << std::endl;
    std::cout << "    One muon cut: " << ThirdCutInt << ". G.E: " <<  (ThirdCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (ThirdCutTrueInt / AllTrueEventsInt) * 100. << ". Purity: " << (ThirdCutTrueInt / ThirdCutInt) * 100. << std::endl;
    std::cout << "    Two protons cut: " << FourthCutInt << ". G.E: " <<  (FourthCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (FourthCutTrueInt / AllTrueEventsInt) * 100. << ". Purity: " << (FourthCutTrueInt / FourthCutInt) * 100. << std::endl;
    std::cout << "    No charged pions cut: " << FifthCutInt << ". G.E: " <<  (FifthCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (FifthCutTrueInt / AllTrueEventsInt) * 100. << ". Purity: " << (FifthCutTrueInt / FifthCutInt) * 100. << std::endl;
    std::cout << "    No neutral pions cut: " << SixthCutInt << ". G.E: " <<  (SixthCutInt / AllRecoEventsInt) * 100. << ". S.E.: " << (SixthCutTrueInt / AllTrueEventsInt) * 100. << ". Purity: " << (SixthCutTrueInt / SixthCutInt) * 100. << std::endl;
    std::cout << std::endl;
    std::cout << "Reconstructed events satisfying signal definition: " << RecoSignalInt << ". Final signal efficiency: " << (RecoTrueSignalInt / AllTrueEventsInt) << ". Purity: " << (RecoTrueSignalInt / RecoSignalInt) * 100. << std::endl;
    std::cout << "Cross check!" << endl;
    std::cout << "Reconstructed true signal: " << RecoTrueSignalInt << ", divided by signal efficiency: " << RecoTrueSignalInt * (AllTrueEventsInt / RecoTrueSignalInt) << std::endl;
    std::cout << "True signal events: " << AllTrueEventsInt << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << std::endl;

    // Data

    std::cout << std::endl;
    std::cout << "============== Data ==================" << std::endl;
    std::cout << "Reconstructed events: " << DataAllRecoEventsInt << std::endl;
    std::cout << std::endl;
    std::cout << "Cuts: " << std::endl;
    std::cout << "    Cosmic cut: " << DataFirstCutInt << ". G.E: " <<  (DataFirstCutInt / DataAllRecoEventsInt) * 100. << std::endl;
    std::cout << "    Vertex in FV cut: " << DataSecondCutInt << ". G.E: " <<  (DataSecondCutInt / DataAllRecoEventsInt) * 100. << std::endl;
    std::cout << "    One muon cut: " << DataThirdCutInt << ". G.E: " <<  (DataThirdCutInt / DataAllRecoEventsInt) * 100. << std::endl;
    std::cout << "    Two protons cut: " << DataFourthCutInt << ". G.E: " <<  (DataFourthCutInt / DataAllRecoEventsInt) * 100.  << std::endl;
    std::cout << "    No charged pions cut: " << DataFifthCutInt << ". G.E: " <<  (DataFifthCutInt / DataAllRecoEventsInt) * 100. << std::endl;
    std::cout << "    No neutral pions cut: " << DataSixthCutInt << ". G.E: " <<  (DataSixthCutInt / DataAllRecoEventsInt) * 100. << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << std::endl;   

}
