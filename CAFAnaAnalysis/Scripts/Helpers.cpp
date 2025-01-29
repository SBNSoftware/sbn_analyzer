// ROOT includes.
#include "TH1D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"

// std includes.
#include <vector>

namespace SelectionHelpers {
    // https://github.com/SBNSoftware/sbnana/blob/245aecf3f422e54820ffcdc53db326e5cd8e8fe0/sbnana/CAFAna/Core/EnsembleSpectrum.cxx#L196C3-L224C4
    void DrawErrorBand(TH1* nom, TGraphAsymmErrors* band, int bandCol, double alpha) {
        if(bandCol == -1) bandCol = nom->GetLineColor()-10; // hopefully a lighter version

        // Check if this pad has already been drawn in
        const bool same = gPad && !gPad->GetListOfPrimitives()->IsEmpty();

        nom->Draw(same ? "hist same" : "hist");

        band->SetFillColorAlpha(bandCol, alpha);
        band->Draw("e2 same");

        nom->Draw("hist same");

        // If we are the first spectrum to draw, scale the y-axis appropriately to
        // fit the error band as well as the nominal
        if(!same){
        double maxY = 0;
        // Don't consider underflow or overflow bins when determining maximum
        for(int i = 1; i < band->GetN()-1; ++i){
            maxY = std::max(maxY, band->GetY()[i] + band->GetErrorYhigh(i));
        }

        // Use non-zero lower value so that SetLogy() still works
        nom->GetYaxis()->SetRangeUser(1e-10, 1.1 * maxY);
        }
        gPad->RedrawAxis();
    }

    // https://github.com/SBNSoftware/sbnana/blob/245aecf3f422e54820ffcdc53db326e5cd8e8fe0/sbnana/CAFAna/Core/Utilities.cxx#L974-L1000
    double FindQuantile(double frac, std::vector<double>& xs) {
        // This turns out to be a much more fraught issue than you would naively
        // expect. This algorithm is equivalent to R-6 here:
        // https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample

        // In principle we could use std::nth_element(). Probably doesn't matter
        // much in practice since this is only for plotting.
        std::sort(xs.begin(), xs.end());

        const int N = xs.size();
        // The index we would ideally be sampling at
        const double h = frac*(N+1);
        // The indices on either side where we have to actually evaluate
        const unsigned int h0 = std::floor(h);
        const unsigned int h1 = std::ceil(h);
        if(h0 == 0) return xs[0]; // Don't underflow indexing
        if(h1 > xs.size()) return xs.back(); // Don't overflow indexing
        // The values at those indices
        const double x0 = xs[h0-1]; // wikipedia is using 1-based indexing
        const double x1 = xs[h1-1];

        if(h0 == h1) return x0;

        // Linear interpolation
        return (h1-h)*x0 + (h-h0)*x1;
    }

    void DrawHistosWithErrorBands(
        TH1D* RecoHisto, 
        TH1D* RecoTrueHisto, 
        TH1D* RecoBkgHisto,
        std::vector<TH1D*> UnivRecoHisto, 
        std::vector<TH1D*> UnivRecoTrueHisto, 
        std::vector<TH1D*> UnivRecoBkgHisto,
        TString dir, 
        TString SystName, 
        TString PlotName,
        double TextSize = 0.06,
        int FontStyle = 132
    ) {
        // Make sure univ histos are of the same length
        assert (UnivRecoHisto.size() == UnivRecoTrueHisto.size());
        assert (UnivRecoTrueHisto.size() == UnivRecoBkgHisto.size());
        int NUniv = UnivRecoHisto.size();

        // Create canvas for plots
        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);

        // Legend for plot with error bands
        TLegend* leg = new TLegend(0.2,0.73,0.75,0.83);
        leg->SetBorderSize(0);
        leg->SetNColumns(3);
        leg->SetTextSize(TextSize*0.8);
        leg->SetTextFont(FontStyle);

        TLegendEntry* legReco = leg->AddEntry(RecoHisto,"Reconstructed","l");
        RecoHisto->SetLineColor(kBlue+2);
        RecoHisto->SetLineWidth(4);

        TLegendEntry* legRecoTrue = leg->AddEntry(RecoTrueHisto,"True","l");
        RecoTrueHisto->SetLineColor(kRed+1);
        RecoTrueHisto->SetLineWidth(4);

        TLegendEntry* legRecoBkg = leg->AddEntry(RecoBkgHisto,"Background","l");
        RecoBkgHisto->SetLineColor(kOrange+7);
        RecoBkgHisto->SetLineWidth(4);

        int n = RecoHisto->GetNbinsX();
        double edges[n+1];
        for (int i = 0; i < n+1; i++) { edges[i] = RecoHisto->GetBinLowEdge(i+1); }

        // Get graph with error bands
        TGraphAsymmErrors* RecoErrorBand = new TGraphAsymmErrors;
        TGraphAsymmErrors* RecoTrueErrorBand = new TGraphAsymmErrors;
        TGraphAsymmErrors* BkgErrorBand = new TGraphAsymmErrors;
        for (int binIdx = 0; binIdx < n + 2; ++binIdx) {
            const double recoxnom = RecoHisto->GetXaxis()->GetBinCenter(binIdx);
            const double recoynom = RecoHisto->GetBinContent(binIdx);
            RecoErrorBand->SetPoint(binIdx, recoxnom, recoynom);

            const double bkgxnom = RecoBkgHisto->GetXaxis()->GetBinCenter(binIdx);
            const double bkgynom = RecoBkgHisto->GetBinContent(binIdx);
            BkgErrorBand->SetPoint(binIdx, bkgxnom, bkgynom);

            const double recotruexnom = RecoTrueHisto->GetXaxis()->GetBinCenter(binIdx);
            const double recotrueynom = RecoTrueHisto->GetBinContent(binIdx);
            RecoTrueErrorBand->SetPoint(binIdx, recotruexnom, recotrueynom);

            const double dx = RecoHisto->GetXaxis()->GetBinWidth(binIdx);

            std::vector<double> recoys;
            for(int iUniv = 0; iUniv < NUniv; ++iUniv){
                recoys.push_back(UnivRecoHisto.at(iUniv)->GetBinContent(binIdx));
            }
            const double recoy0 = SelectionHelpers::FindQuantile(.5-0.6827/2, recoys);
            const double recoy1 = SelectionHelpers::FindQuantile(.5+0.6827/2, recoys);

            std::vector<double> bkgys;
            for(int iUniv = 0; iUniv < NUniv; ++iUniv){
                bkgys.push_back(UnivRecoBkgHisto.at(iUniv)->GetBinContent(binIdx));
            }
            const double bkgy0 = SelectionHelpers::FindQuantile(.5-0.6827/2, bkgys);
            const double bkgy1 = SelectionHelpers::FindQuantile(.5+0.6827/2, bkgys);

            std::vector<double> recotrueys;
            for(int iUniv = 0; iUniv < NUniv; ++iUniv){
                recotrueys.push_back(UnivRecoTrueHisto.at(iUniv)->GetBinContent(binIdx));
            }
            const double recotruey0 = SelectionHelpers::FindQuantile(.5-0.6827/2, recotrueys);
            const double recotruey1 = SelectionHelpers::FindQuantile(.5+0.6827/2, recotrueys);

            RecoErrorBand->SetPointError(binIdx, dx/2, dx/2, std::max(recoynom-recoy0, 0.), std::max(recoy1-recoynom, 0.));
            BkgErrorBand->SetPointError(binIdx, dx/2, dx/2, std::max(bkgynom-bkgy0, 0.), std::max(bkgy1-bkgynom, 0.));
            RecoTrueErrorBand->SetPointError(binIdx, dx/2, dx/2, std::max(recotrueynom-recotruey0, 0.), std::max(recotruey1-recotrueynom, 0.));
        }
        PlotCanvas->cd();
        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.17);
        PlotCanvas->SetRightMargin(0.05);
        PlotCanvas->SetBottomMargin(0.16);

        double imax = RecoHisto->GetMaximum();
        double YAxisRange = 1.3*imax;
        RecoHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoTrueHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);
        RecoBkgHisto->GetYaxis()->SetRangeUser(0.,YAxisRange);

        RecoHisto->Draw("hist");
        SelectionHelpers::DrawErrorBand(RecoHisto, RecoErrorBand, -1, 1);
        RecoTrueHisto->Draw("hist same");
        SelectionHelpers::DrawErrorBand(RecoTrueHisto, RecoTrueErrorBand, -1, 1);
        RecoBkgHisto->Draw("hist same");
        SelectionHelpers::DrawErrorBand(RecoBkgHisto, BkgErrorBand, -1, 1);
        leg->Draw();
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+SystName+"/"+PlotName+".pdf");
        delete PlotCanvas; delete leg;
    }

    void GetFracCovAndCorrMatrix(
        TH1* RecoHisto,
        TH2* CovMatrix,
        TH2* FracCovMatrix,
        TH2* CorrMatrix,
        int n
    ) {
        for (int x = 1; x < n + 1; x++) {
            double XEventRateCV = RecoHisto->GetBinContent(x);
            for (int y = 1; y < n + 1; y++) {
                double YEventRateCV = RecoHisto->GetBinContent(y);
                double CovBinValue = CovMatrix->GetBinContent(x,y);
                double XBinValue = CovMatrix->GetBinContent(x,x);
                double YBinValue = CovMatrix->GetBinContent(y,y);

                // Fill frac cov matrix
                double FracValue = (XEventRateCV == 0. || YEventRateCV == 0.) ? 1e-14 : CovBinValue / (XEventRateCV * YEventRateCV);
                if (TMath::Abs(FracValue) < 1e-14) FracValue = 1e-14;
                FracCovMatrix->SetBinContent(x, y, FracValue);

                // Fill corr matrix
                double CorrValue = (XBinValue == 0. || YBinValue == 0.) ? 1e-14 : CovBinValue / (TMath::Sqrt(XBinValue) * TMath::Sqrt(YBinValue));
                if (TMath::Abs(CorrValue) < 1e-14) CorrValue = 1e-14;
                CorrMatrix->SetBinContent(x, y, CorrValue);
            }
        }
    }

    void GetCovMatrix(
        TH1* RecoHisto,
        std::vector<TH1D*> UnivRecoHisto,
        TH2* CovMatrix,
        int n
    ) {
        int NUniv = UnivRecoHisto.size();

        // Create covariance matrices 
        for (int iUniv = 0; iUniv < NUniv; ++iUniv) {
            for (int x = 1; x < n + 1; x++) {
                double XEventRateCV = RecoHisto->GetBinContent(x);
                double XEventRateVar = UnivRecoHisto[iUniv]->GetBinContent(x);
                for (int y = 1; y < n + 1; y++) {
                    double YEventRateCV = RecoHisto->GetBinContent(y);
                    double YEventRateVar = UnivRecoHisto[iUniv]->GetBinContent(y);

                    double Value = ((XEventRateVar - XEventRateCV) * (YEventRateVar - YEventRateCV)) / NUniv;
                    // std::cout << Value << std::endl;
                    // std::cout << XEventRateVar << "  " << XEventRateCV << "  " << XEventRateVar - XEventRateCV << std::endl;
                    // std::cout << YEventRateVar << "  " << YEventRateCV << "  " << YEventRateVar - YEventRateCV << std::endl;
                    if (TMath::Abs(Value) < 1e-14) Value = 1e-14;
                    // std::cout << Value << std::endl;
                    // std::cout << std::endl;

                    CovMatrix->Fill(
                        RecoHisto->GetXaxis()->GetBinCenter(x),
                        RecoHisto->GetXaxis()->GetBinCenter(y),
                        Value
                    );
                }
            }
        }
    }

    void DrawMatrices(
        TH1D* RecoHisto, 
        TH1D* RecoTrueHisto, 
        TH1D* RecoBkgHisto,
        std::vector<TH1D*> UnivRecoHisto, 
        std::vector<TH1D*> UnivRecoTrueHisto, 
        std::vector<TH1D*> UnivRecoBkgHisto,
        TString dir, 
        TString SystName, 
        TString VarLabel,
        TString PlotName,
        TFile* SaveFile,
        double TextSize = 0.06,
        int FontStyle = 132
    ) {
        // Make sure univ histos are of the same length
        assert (UnivRecoHisto.size() == UnivRecoTrueHisto.size());
        assert (UnivRecoTrueHisto.size() == UnivRecoBkgHisto.size());

        int n = RecoHisto->GetXaxis()->GetNbins();
        double edges[n+1];
        for (int i = 0; i < n+1; i++) { edges[i] = RecoHisto->GetBinLowEdge(i+1); }

        // Create canvas for plots
        TCanvas* PlotCanvas = new TCanvas("Selection","Selection",205,34,1124,768);

        // Create covariance matrix
        std::string CovName = "Cov" + (std::string) SystName;
        TH2* CovMatrix = new TH2D(
            (CovName + (std::string)PlotName).c_str(),
            CovName.c_str(),
            n, edges,
            n, edges
        );

        // Create fractional covariance matrix
        std::string FracCovName = "FracCov" + (std::string) SystName;
        TH2* FracCovMatrix = new TH2D(
            (FracCovName + (std::string)PlotName).c_str(),
            FracCovName.c_str(),
            n, edges,
            n, edges
        );

        // Create correlation matrix
        std::string CorrName = "Corr" + (std::string) SystName;
        TH2* CorrMatrix = new TH2D(
            (CorrName + (std::string)PlotName).c_str(),
            CorrName.c_str(),
            n, edges,
            n, edges
        );

        // Create covariance matrix
        SelectionHelpers::GetCovMatrix(
            RecoHisto,
            UnivRecoHisto,
            CovMatrix,
            n
        );

        // Fill in frac cov and corr matrix
        SelectionHelpers::GetFracCovAndCorrMatrix(
            RecoHisto,
            CovMatrix,
            FracCovMatrix,
            CorrMatrix,
            n
        );

        // Plot cov matrix
        double CovMin = CovMatrix->GetMinimum();
        double CovMax = CovMatrix->GetMaximum();
        CovMatrix->GetZaxis()->SetRangeUser(CovMin,CovMax); // set the ranges accordingly, for frac cov should be [0,100], for corr matrices [-1,1]
        CovMatrix->GetZaxis()->CenterTitle();
        CovMatrix->GetZaxis()->SetTitleFont(FontStyle);
        CovMatrix->GetZaxis()->SetTitleSize(TextSize);
        CovMatrix->GetZaxis()->SetLabelFont(FontStyle);
        CovMatrix->GetZaxis()->SetLabelSize(TextSize);
        CovMatrix->GetZaxis()->SetNdivisions(5);

        CovMatrix->GetXaxis()->SetTitle(("bin i " + (std::string)VarLabel).c_str());
        CovMatrix->GetYaxis()->SetTitle(("bin j " + (std::string)VarLabel).c_str());

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        CovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+SystName+"/Cov"+PlotName+".pdf");

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

        FracCovMatrix->GetXaxis()->SetTitle(("bin i " + (std::string)VarLabel).c_str());
        FracCovMatrix->GetYaxis()->SetTitle(("bin j " + (std::string)VarLabel).c_str());

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        FracCovMatrix->Draw("colz");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+SystName+"/FracCov"+PlotName+".pdf");

        // Plot correlation matrix
        CorrMatrix->GetZaxis()->SetRangeUser(-1,1);
        CorrMatrix->GetZaxis()->CenterTitle();
        CorrMatrix->GetZaxis()->SetTitleFont(FontStyle);
        CorrMatrix->GetZaxis()->SetTitleSize(TextSize);
        CorrMatrix->GetZaxis()->SetLabelFont(FontStyle);
        CorrMatrix->GetZaxis()->SetLabelSize(TextSize);
        CorrMatrix->GetZaxis()->SetNdivisions(5);

        CorrMatrix->GetXaxis()->SetTitle(("bin i " + (std::string)VarLabel).c_str());
        CorrMatrix->GetYaxis()->SetTitle(("bin j " + (std::string)VarLabel).c_str());

        PlotCanvas->cd();

        PlotCanvas->SetTopMargin(0.13);
        PlotCanvas->SetLeftMargin(0.15);
        PlotCanvas->SetRightMargin(0.15);
        PlotCanvas->SetBottomMargin(0.16);

        CorrMatrix->Draw("colz text");
        PlotCanvas->SaveAs(dir+"/Figs/CAFAna/Uncertainties/"+SystName+"/Corr"+PlotName+".pdf");

        // Save objects
        SaveFile->WriteObject(CovMatrix, PlotName+"_cov");
        SaveFile->WriteObject(FracCovMatrix, PlotName+"_fraccov");
        SaveFile->WriteObject(CorrMatrix, PlotName+"_corr");

        delete PlotCanvas;
    }
}
