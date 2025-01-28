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

// std includes.
#include <vector>
#include <memory>

// Definitions for Vars and Cuts.
#include "Definitions.h"
#include "Utils.h"

// Utils includes.
#include "../../Utils/Constants.h"

using namespace std;
using namespace ana;
using namespace Constants;

void Selection() {
  bool debug=true;
  // Set defaults and load tools
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  if(debug) cout<<"l34"<<endl;
  int FontStyle = 132;
  double TextSize = 0.06;	
  if(debug)cout<<"l37"<<endl;
  // The SpectrumLoader object handles the loading of CAFs and the creation of Spectrum.
  SpectrumLoader NuLoader(InputFiles);
  if(debug)cout<<"l40"<<endl;
  // We now create overlaid plots for several reconstructed variables and three lines:
  //     1. all selected reconstructed events
  //     2. reco signal events
  //     3. reco background events

  // Root file to store objects in
  TString RootFilePath = "/exp/sbnd/data/users/" + (TString)UserName + "/CAFAnaOutput/Selection.root";
  TFile* SaveFile = new TFile(RootFilePath, "RECREATE");
  if(debug)  cout<<"l49"<<endl;
  // Construct all spectra
  std::vector<std::tuple<
		std::unique_ptr<Spectrum>,
		std::unique_ptr<Spectrum>,
		std::unique_ptr<Spectrum>
		>> Spectra;
 
 static const std::vector<std::tuple<Var, Var, TruthVar>> Vars=GetSISVars();
 static const std::vector<Binning> VarBins = GetSISBins();
 static const std::vector<std::string>VarLabels= GetSISVarLabels();

  

  cout<<"vars:"<<VarsOld.size()<<"\t VarLabels"<<VarLabels.size()<<"\t VarBins"<<VarBins.size()<<endl;
 
  for (std::size_t i = 0; i < Vars.size(); i++) {
    cout<<"l57 i:"<<i<<"  VarBins.NBins(): "<<VarBins[i].NBins()<<" ["<<VarBins[i].Min()<<","<<VarBins[i].Max()<<"]"<<endl;
    auto RecoSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsSignal); 
    auto RecoTrueSignals = std::make_unique<Spectrum> (VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsTrueReco); 
    auto RecoBkgSignals = std::make_unique<Spectrum>(VarLabels.at(i), VarBins.at(i), NuLoader, std::get<0>(Vars.at(i)), kNoSpillCut, kRecoIsBackground); 
    Spectra.push_back({std::move(RecoSignals), std::move(RecoTrueSignals), std::move(RecoBkgSignals)});
    if(debug)cout<<"l63 i:"<<i<<endl;
  }
  if(debug)cout<<"l63"<<endl;
  // We now create spectra that will help us get the efficiency and purity data for each of the cuts
  if(debug) cout<<"l63"<<endl;  
  // Spectrum with all events
  /*    Spectrum sAllEvents("AllEvents", bEventCount, NuLoader, kTrueEventCount, kNoTruthCut, kNoSpillCut);
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
  Spectrum sRecoTrueSignal("RecoTrueSignal", bEventCount, NuLoader, kEventCount, kNoSpillCut, kRecoIsTrueReco); */
  if(debug) cout<<"l94"<<endl;  
  NuLoader.Go();
  if(debug) cout<<"l96"<<endl;  
  // Loop over variables
  double mcPOT;
  for (std::size_t i = 0; i < Vars.size(); i++) {
    auto& [RecoSignals, RecoTrueSignals, RecoBkgSignals] = Spectra.at(i);
    
    if(i==0) mcPOT=RecoSignals->POT();
    
    TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);
    TH1D* RecoHisto = RecoSignals->ToTH1(TargetPOT);
    TH1D* RecoTrueHisto = RecoTrueSignals->ToTH1(TargetPOT);
    TH1D* RecoBkgHisto = RecoBkgSignals->ToTH1(TargetPOT);

    // Manage under/overflow bins
    RecoHisto->SetBinContent(RecoHisto->GetNbinsX(), RecoHisto->GetBinContent(RecoHisto->GetNbinsX()) + RecoHisto->GetBinContent(RecoHisto->GetNbinsX() + 1));
    RecoTrueHisto->SetBinContent(RecoTrueHisto->GetNbinsX(), RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX()) + RecoTrueHisto->GetBinContent(RecoTrueHisto->GetNbinsX() + 1));
    RecoBkgHisto->SetBinContent(RecoBkgHisto->GetNbinsX(), RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX()) + RecoBkgHisto->GetBinContent(RecoBkgHisto->GetNbinsX() + 1));

    RecoHisto->SetBinContent(1, RecoHisto->GetBinContent(0) + RecoHisto->GetBinContent(1));
    RecoTrueHisto->SetBinContent(1, RecoTrueHisto->GetBinContent(0) + RecoTrueHisto->GetBinContent(1));
    RecoBkgHisto->SetBinContent(1, RecoBkgHisto->GetBinContent(0) + RecoBkgHisto->GetBinContent(1));

    PlotCanvas->SetTopMargin(0.13);
    PlotCanvas->SetLeftMargin(0.17);
    PlotCanvas->SetRightMargin(0.05);
    PlotCanvas->SetBottomMargin(0.16);

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

    PlotCanvas->cd();
    RecoHisto->Draw("hist same");
    RecoTrueHisto->Draw("hist same");
    RecoBkgHisto->Draw("hist same");
    leg->Draw();

    // Save as pdf
    PlotCanvas->SaveAs(dir_figs+"/Figs/CAFAna/"+PlotNames[i]+".pdf");

    // Save to root file
    SaveFile->WriteObject(RecoHisto, PlotNames[i]+"_reco");
    SaveFile->WriteObject(RecoTrueHisto, PlotNames[i]+"_reco_true");
    SaveFile->WriteObject(RecoBkgHisto, PlotNames[i]+"_bkg");

    delete PlotCanvas;
  }
  // Close file
  WritePOT(SaveFile, mcPOT);
  SaveFile->Close();

  // Get histograms for all cuts
  /*TH1D* AllEventsHisto = sAllEvents.ToTH1(TargetPOT);
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

    // Get integrals for all cuts
    double AllEventsInt = AllEventsHisto->Integral();
    double AllRecoEventsInt = AllRecoEventsHisto->Integral();
    double AllTrueEventsInt = AllTrueEventsHisto->Integral();
    double AllTrueRecoEventsInt = AllTrueRecoEventsHisto->Integral();

    double FirstCutInt = FirstCutHisto->Integral();
    double FirstCutTrueInt = FirstCutTrueHisto->Integral();

    double SecondCutInt = SecondCutHisto->Integral();
    double SecondCutTrueInt = SecondCutTrueHisto->Integral();

    double ThirdCutInt = ThirdCutHisto->Integral();
    double ThirdCutTrueInt = ThirdCutTrueHisto->Integral();

    double FourthCutInt = FourthCutHisto->Integral();
    double FourthCutTrueInt = FourthCutTrueHisto->Integral();

    double FifthCutInt = FifthCutHisto->Integral();
    double FifthCutTrueInt = FifthCutTrueHisto->Integral();

    double SixthCutInt = SixthCutHisto->Integral();
    double SixthCutTrueInt = SixthCutTrueHisto->Integral();

    double RecoSignalInt = RecoSignalHisto->Integral();
    double RecoTrueSignalInt = RecoTrueSignalHisto->Integral();

    // Print results
    std::cout << std::endl;
    std::cout << "================================" << std::endl;
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
    std::cout << "Cross check. Reconstructed true signal: " << RecoTrueSignalInt << ", divided by signal efficiency: " << RecoTrueSignalInt * (AllTrueEventsInt / RecoTrueSignalInt) << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << std::endl;
  */
}
