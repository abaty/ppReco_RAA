#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  
  int centBounds[7] = {0,5,10,30,50,70,90};
  for(int i = 0 ; i<6; i++){
  const char * histName = "PbPbTrackSpectrum";
  TFile * fold = TFile::Open("Spectra_Jun23_FullStatsNoJ40NewCorr.root","read");
  TFile * fnew = TFile::Open("Spectra_Jun9_noChi2Cut.root","read");
  TH1D * num = (TH1D*)fnew->Get(Form("%s_%d_%d",histName,centBounds[i],centBounds[i+1]));
  TH1D * den = (TH1D*)fold->Get(Form("%s_%d_%d",histName,centBounds[i],centBounds[i+1]));
  //den->Add(num);
  num->Divide(den);
  //num->Scale(2);
  num->GetYaxis()->SetTitle("(New Table)/(Old Table)");
  num->GetYaxis()->SetRangeUser(0.7,1.3);
  num->GetXaxis()->SetRangeUser(0.7,350);
  num->Print("All");
  num->Draw();
  c1->SetLogx();
  c1->SaveAs(Form("plots/comparisonPlots/Jun27_noChi2_%d_%d.png",centBounds[i],centBounds[i+1]));
  c1->SaveAs(Form("plots/comparisonPlots/Jun27_noChi2_%d_%d.pdf",centBounds[i],centBounds[i+1]));
  fnew->Close();
  fold->Close();
  }  

  /*const char * histName = "PbPbTrackSpectrum_0_5";
  TFile * fnew = TFile::Open("Spectra_Jun9_noChi2Cut.root","read");
  TFile * fold = TFile::Open("Spectra_Jun9_withChi2Cut.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  //den->Add(num);
  num->Divide(den);
  //num->Scale(2);
  num->GetYaxis()->SetTitle("(No Chi2)/(With Chi2)");
  num->GetYaxis()->SetRangeUser(0.7,1.5);
  num->GetXaxis()->SetRangeUser(0.7,350);
  num->Print("All");
  num->Draw();
  c1->SetLogx();
  c1->SaveAs("plots/comparisonPlots/Chi2CutTest_PbPb_0_5.png");
  c1->SaveAs("plots/comparisonPlots/Chi2CutTest_PbPb_0_5.pdf");
  */

  //c1->SaveAs("plots/comparisonPlots/ppChargeFraction2.C");
}
