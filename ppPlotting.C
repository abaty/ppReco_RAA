#include "TBox.h"
#include "Settings.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAttAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "THStack.h"

void makePlotsPP(Settings s)
{
  for(int j = 0; j<s.nTriggers; j++)
  {
    s.ppByTrigger[j]->SetLineColor(j+1);
    s.ppByTrigger[j]->SetLineWidth(1);
    s.ppByTrigger[j]->SetMarkerColor(j+1);
    s.ppByTrigger[j]->SetFillColor(j+1);
    s.ppUsedByTrigger[j]->SetLineColor(kBlack);
    s.ppUsedByTrigger[j]->SetLineWidth(2);
    s.ppUsedByTrigger[j]->SetMarkerColor(j+1);
    s.ppUsedByTrigger[j]->SetMarkerSize(0);
    s.ppUsedByTrigger[j]->SetFillColor(j+1);
    s.ppJetsByTrigger[j]->SetLineColor(j+1);
    s.ppJetsByTrigger[j]->SetLineWidth(1);
    s.ppJetsByTrigger[j]->SetMarkerColor(j+1);
    s.ppJetsByTrigger[j]->SetMarkerSize(0.8);
    s.ppJetsByTrigger[j]->SetFillColor(j+1);
  }
  for(int j = 0; j<s.nTriggers_trk; j++)
  {
    s.ppByTrigger_trk[j]->SetLineColor(j+1);
    s.ppByTrigger_trk[j]->SetLineWidth(1);
    s.ppByTrigger_trk[j]->SetMarkerColor(j+1);
    s.ppByTrigger_trk[j]->SetFillColor(j+1);
    s.ppUsedByTrigger_trk[j]->SetLineColor(kBlack);
    s.ppUsedByTrigger_trk[j]->SetLineWidth(2);
    s.ppUsedByTrigger_trk[j]->SetMarkerColor(j+1);
    s.ppUsedByTrigger_trk[j]->SetMarkerSize(0);
    s.ppUsedByTrigger_trk[j]->SetFillColor(j+1);
    s.ppMaxtrkByTrigger[j]->SetLineColor(j+1);
    s.ppMaxtrkByTrigger[j]->SetLineWidth(1);
    s.ppMaxtrkByTrigger[j]->SetMarkerColor(j+1);
    s.ppMaxtrkByTrigger[j]->SetMarkerSize(0.8);
    s.ppMaxtrkByTrigger[j]->SetFillColor(j+1);
  }
//******************************************************************************************************************
//************************************************JET TRIGGER PLOTS**********************************************************
//******************************************************************************************************************
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  c1->SetLogy();
  s.ppJets->Scale(100);
  float Ymin = 0.00000000001;
  float Ymax = 10;
  s.ppJets->GetYaxis()->SetRangeUser(Ymin,Ymax);
  s.ppJets->Draw("h");
  for(int i = 0; i<s.nTriggers; i++) s.ppJetsByTrigger[i]->Draw("same");
  TLegend * leg = new TLegend(0.47,0.62,0.87,0.92);
  leg->AddEntry(s.ppJets,"Jet Spectrum (x100)","l");
  leg->AddEntry(s.ppJetsByTrigger[0],"MB (I)","p");
  leg->AddEntry(s.ppJetsByTrigger[1],"jet 40 (II)","p");
  leg->AddEntry(s.ppJetsByTrigger[2],"jet 60 (III)","p");
  leg->AddEntry(s.ppJetsByTrigger[3],"jet 80 (IV)","p");
  leg->AddEntry((TObject*)0,"ak4Calo Jets, |#eta|<2","");
  leg->Draw("same");
  c1->SaveAs("plots/png/ppJets_FullSpectrum.png"); 
  c1->SaveAs("plots/pdf/ppJets_FullSpectrum.pdf"); 

  s.ppJets->GetXaxis()->SetRangeUser(20,200);
  Ymin = 0.0000001;
  Ymax = 100;
  s.ppJets->GetYaxis()->SetRangeUser(0.0000001,100);
  s.ppJets->Draw("h");
  for(int i = 0; i<s.nTriggers; i++) s.ppJetsByTrigger[i]->Draw("same");
  leg->Draw("same"); 
  TLine * l[3];
  for(int i = 0; i<3; i++) 
  {
    l[i] = new TLine(s.triggerBins[i+1],Ymin,s.triggerBins[i+1],Ymax); 
    l[i]->SetLineWidth(2);
    l[i]->SetLineStyle(2);
    l[i]->SetLineColor(1);
    l[i]->Draw("same");
  }
  TLatex * lat = new TLatex(1,1,"test");
  lat->DrawLatex(45,Ymin*3,"I");
  lat->DrawLatex(65,Ymin*3,"II");
  lat->DrawLatex(85,Ymin*3,"III");
  lat->DrawLatex(105,Ymin*3,"IV");
  c1->SaveAs("plots/png/ppJets_FullSpectrum_XZoom.png"); 
  c1->SaveAs("plots/pdf/ppJets_FullSpectrum_XZoom.pdf"); 

  c1->SetLogy(0);
  for(int i = 0; i<s.nTriggers-1; i++) s.ppJetsByTrigger[s.nTriggers-1-i]->Divide(s.ppJetsByTrigger[s.nTriggers-2-i]);
  s.ppJetsByTrigger[1]->GetYaxis()->SetRangeUser(0,2);
  s.ppJetsByTrigger[1]->GetXaxis()->SetRangeUser(20,140);
  s.ppJetsByTrigger[1]->GetXaxis()->SetTitle("Leading jet p_{T}");
  s.ppJetsByTrigger[1]->Draw();
  s.ppJetsByTrigger[2]->Draw("same");
  s.ppJetsByTrigger[3]->Draw("same");
  TLegend * leg2 = new TLegend(0.2,0.6,0.5,0.9);
  leg2->AddEntry(s.ppJetsByTrigger[1],"Jet40/MB","p");
  leg2->AddEntry(s.ppJetsByTrigger[2],"Jet60/Jet40","p");
  leg2->AddEntry(s.ppJetsByTrigger[3],"Jet80/Jet60","p");
  leg2->Draw("same");
  c1->SaveAs("plots/png/JetRelativeTurnOnes.png");
  c1->SaveAs("plots/pdf/JetRelativeTurnOnes.pdf");
  
  c1->SetLogy();
  c1->SetLogx();
  s.pp->SetMarkerSize(0.8);
  s.pp->GetXaxis()->SetRangeUser(0.7,400);
  s.pp->Draw();
  s.ppUsedByTrigger[0]->SetFillColor(kGray);
  s.ppUsedByTrigger[3]->SetFillColor(kCyan+2);
  s.ppUsedByTrigger[2]->Add(s.ppUsedByTrigger[3]);
  s.ppUsedByTrigger[1]->Add(s.ppUsedByTrigger[2]);
  s.ppUsedByTrigger[0]->Add(s.ppUsedByTrigger[1]);
  for(int i = 0; i<s.nTriggers; i++) s.ppUsedByTrigger[i]->Draw("HIST same");
  s.pp->Draw("sameaxis");
  s.pp->Draw("same");
  leg->Clear();
  leg->AddEntry(s.pp,"pp track Spectrum","p");
  leg->AddEntry(s.ppUsedByTrigger[0],"MB trigger","f");
  leg->AddEntry(s.ppUsedByTrigger[1],"Jet40 trigger","f");
  leg->AddEntry(s.ppUsedByTrigger[2],"Jet60 trigger","f");
  leg->AddEntry(s.ppUsedByTrigger[3],"Jet80 trigger","f");
  leg->AddEntry((TObject*)0,"|#eta|<1","");
  leg->Draw("same");
  
  
 
  c1->SaveAs("plots/png/ppTrack_FullSpectrum.png");
  c1->SaveAs("plots/pdf/ppTrack_FullSpectrum.pdf");
  c1->SetLogy(0);
  for(int i = 0; i<s.nTriggers; i++) s.ppUsedByTrigger[i]->Divide(s.pp);
  s.ppUsedByTrigger[0]->GetYaxis()->SetRangeUser(0,2);
  s.ppUsedByTrigger[0]->GetYaxis()->SetTitle("Relative Contribution to Bin");
  s.ppUsedByTrigger[0]->Draw("HIST");
  for(int i = 1; i<s.nTriggers; i++) s.ppUsedByTrigger[i]->Draw("HIST same");
  s.ppUsedByTrigger[0]->Draw("sameaxis");
  leg->Draw("same");
  c1->SaveAs("plots/png/ppTrack_FullSpectrum_relativeContribution.png");
  c1->SaveAs("plots/pdf/ppTrack_FullSpectrum_relativeContribution.pdf");
  c1->SetLogy();


  TH1D * jtVsTrk = (TH1D*)s.pp->Clone("jtVsTrkRatio");
  
  //xt plot
  TFile * outF = TFile::Open("Spectra.root","UPDATE");
  TH1D * pp_perMBTrigger_xt = (TH1D*)s.pp_perMBTrigger->Clone("pp_perMBTrigger_xt");
  pp_perMBTrigger_xt->Scale(TMath::Power(5020,4.9));
  TH1D * pp_xt = new TH1D("ppTrackSpectrum_xt",";x_{T};#sqrt{s}^{4.9}E#frac{d^{3}#sigma}{d^{3}p} (mb/GeV^{2})",s.ntrkBins,s.xt_xtrkbins);
  for(int i = 1; i< pp_xt->GetSize()-1; i++){
    pp_xt->SetBinContent(i,pp_perMBTrigger_xt->GetBinContent(i));
    pp_xt->SetBinError(i,pp_perMBTrigger_xt->GetBinError(i));
  }
  pp_xt->Draw();
  pp_xt->Print("All");
  pp_xt->Write();
  c1->SaveAs("plots/png/ppTrack_xt_FullSpectrum.png");
  c1->SaveAs("plots/pdf/ppTrack_xt_FullSpectrum.pdf");
  outF->Close();


  //RpPb bin shift correction
  
  //interpolation points
  float ppintyvals[37+1] = {4.1202/0.5946, 2.8116/0.6211, 1.9778/0.6552, 1.4330/0.6984, 1.0405/0.7219, 0.7719/0.7515, 0.5816/0.7809, 0.3893/0.825, 0.23381/0.866, 0.14395/0.901, 0.09103/0.925, 0.05937/0.965, 0.03906/0.984, 0.014787/1.023, 0.003806/1.052, 0.001181/1.056, 0.0004290/1.048, 0.0001787/1.054, 0.00008152/1.031, 0.00002216/1.023, 0.000004653/1.036, 0.000001402/1.054, 0.0000003180/1.072, 0.00000006850/1.142, 0.00000001971/1.189, 0.000000006013/1.259, 0.000000001910/1.308, 0.0000000007181/1.342, 0.0000000002083/1.382,0.00000000005311/1.407,0.00000000001639/1.363,0.000000000005354/1.381,0.000000000001709/1.316,0,0 };//first 2 points were 4.1202/0.5946*0.999, 2.8116/0.6211*1.003
  //float ppintyvals_pPb[s.ntrkBins+1] = {4.1202*0.999, 2.8116*1.003, 1.9778*1.002, 1.4330*0.994 , 1.0405*1.002, 0.7719*1.002, 0.5816*1.002, 0.3893*1.001, 0.23381*1.000, 0.14395*1.004, 0.09103*1.005, 0.05937*0.997, 0.03906*0.998, 0.014787*1.062, 0.003806*1.042, 0.001181*1.033, 0.0004290*1.024, 0.0001787*1.018, 0.00008152*1.015, 0.00002216* 1.109, 0.000004653*1.061, 0.000001402*1.040, 0.0000003180*1.111, 0.00000006850*1.065, 0.00000001971*1.044, 0.000000006013*1.051, 0.000000001910*1.033, 0.0000000007181*1.024, 0.0000000002083*1.078,0.00000000005311*1.05,0.00000000001639*1.034386,0.000000000005354*1.047316,0.000000000001709*1.032760,0,0 };
  float ppintyvals_pPb[37+1] = { 4.1202, 2.8116,1.9778, 1.4330 , 1.0405, 0.7719, 0.5816, 0.3893, 0.23381, 0.14395, 0.09103, 0.05937, 0.03906, 0.014787, 0.003806, 0.001181, 0.0004290, 0.0001787, 0.00008152, 0.00002216, 0.000004653, 0.000001402, 0.0000003180, 0.00000006850, 0.00000001971, 0.000000006013, 0.000000001910, 0.0000000007181, 0.0000000002083,0.00000000005311,0.00000000001639,0.000000000005354,0.000000000001709,0,0 };//first 2 points were 4.1202, 2.8116
  float ppintyvals_RpPb[37+1] = { 0.5946, 0.6211, 0.6552, 0.6984,0.7219,0.7515, 0.7809, 0.825, 0.866, 0.901, 0.925, 0.965, 0.984, 1.023, 1.052,1.056, 1.048, 1.054, 1.031, 1.023, 1.036, 1.054, 1.072, 1.142, 1.189, 1.259, 1.308, 1.342,1.382,1.407,1.363,1.381,1.316,0,0 };//first 2 points were 0.5946, 0.6211
  float ppintyerrRpPb[37+1] = {0.0027/0.5946,0.0028/0.6211, 0.0030/0.6552,0.0032/0.6984,0.0033/0.7219,0.0034/0.7515,0.0036/0.7809,0.003/0.825,0.003/0.866,0.003/0.901,0.003/0.925,0.003/0.965,0.003/0.984,0.002/1.023,0.002/1.052,0.002/1.056,0.002/1.048,0.003/1.054,0.003/1.031,0.003/1.023,0.005/1.036,0.008/1.054,0.001/1.072,0.002/1.142,0.001/1.189,0.002/1.259,0.002/1.308,0.004/1.342,0.004/1.382,0.008/1.407,0.013/1.363,0.019/1.381,0.030/1.316,0,0 };//first 2 points were 0.0027/0.5946,0.0028/0.6211
  float ppintyerr_pPb[37+1] = {0.0188/4.1202,0.0128/2.8116, 0.0090/1.9778,0.0065/1.4330,0.0047/1.0405,0.0035/0.7719,0.0026/0.5816,0.0013/0.3893,0.00076/0.23381,0.00047/0.14395,0.00029/0.09103,0.00019/0.05937,0.00013/0.03906,0.000026/0.014787,0.000007/0.003806,0.000002/0.001181,0.0000009/0.0004290,0.0000004/0.0001787,0.00000025/0.00008152,0.00000025/0.00002216,0.00000006/0.000004653,0.000000023/0.000001402,0.000000011/0.0000003180,0.0000000004/0.00000006850,0.00000000011/0.00000001971,0.00000000002/0.000000006013,0.000000000009/0.000000001910,0.000000000004/0.0000000007181,0.0000000000020/0.0000000002083,0.0000000000007/0.00000000005311,0.00000000000031/0.00000000001639,0.00000000000016/0.000000000005354,0.000000000000073/0.000000000001709,0,0 };//first 2 points were 0.0188/4.1202,0.0128/2.8116
  TH1D * RpPb = new TH1D("pp_interpolation",";pt;E#frac{d^{3}N}{dPtdyd#phi}",s.ntrkBins,s.xtrkbins);
  TH1D * RpPbSpectrum = new TH1D("RpPb_interp",";pt;RpPb",s.ntrkBins,s.xtrkbins);
  TH1D * pPbSpectrum = new TH1D("pPb_data",";pt;E#frac{d^{3}N}{dPtdyd#phi}",s.ntrkBins,s.xtrkbins);
  for(int i=1; i<s.ntrkBins+1; i++) RpPb->SetBinContent(i,ppintyvals[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) RpPb->SetBinError(i,TMath::Power(TMath::Abs(TMath::Power(ppintyerrRpPb[i-1],2)-TMath::Power(ppintyerr_pPb[i-1],2)),0.5)*ppintyvals[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) RpPbSpectrum->SetBinContent(i,ppintyvals_RpPb[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) RpPbSpectrum->SetBinError(i,ppintyvals_RpPb[i-1]*ppintyerrRpPb[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) pPbSpectrum->SetBinContent(i,ppintyvals_pPb[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) pPbSpectrum->SetBinError(i,ppintyvals_pPb[i-1]*ppintyerr_pPb[i-1]);
  
  //float ppScale = RpPb->GetBinContent(17)/s.pp->GetBinContent(17);
  //s.pp->Scale(ppScale); 

  TFile * binShift = TFile::Open("BinningAndResolutionCorrection_TrackTrigger.root","read");
  TH1D *hBinningAndResol_extendedPt = (TH1D*)binShift->Get("hPt_pseudo2_copy1");
  hBinningAndResol_extendedPt->SetDirectory(0);
  binShift->Close();
  /*TH1D *pp_BinShift = (TH1D*)s.pp->Clone("pp_BinShift"); 
  for(int i = 2; i<hBinningAndResol_extendedPt->GetSize()+1; i++)
  {
    pp_BinShift->SetBinContent(i-1,pp_BinShift->GetBinContent(i-1)/hBinningAndResol_extendedPt->GetBinContent(i));
    pp_BinShift->SetBinError(i-1,pp_BinShift->GetBinError(i-1)/hBinningAndResol_extendedPt->GetBinError(i));
  }*/
  TFile * outBinF = TFile::Open("Spectra.root","update");
  hBinningAndResol_extendedPt->Write("RpPb_binShift_Correction");
  TH1D *pp_BinShift = (TH1D*)s.pp_perMBTrigger->Clone("pp_BinShift");
  TH1D *pp_BinShift_MB = (TH1D*)s.ppByTrigger[0]->Clone("pp_BinShift_MB");
  for(int i = 2; i<hBinningAndResol_extendedPt->GetSize()+1; i++)
  {
    pp_BinShift->SetBinContent(i-1,pp_BinShift->GetBinContent(i-1)/hBinningAndResol_extendedPt->GetBinContent(i));
    pp_BinShift->SetBinError(i-1,pp_BinShift->GetBinError(i-1)/hBinningAndResol_extendedPt->GetBinContent(i));
    pp_BinShift_MB->SetBinContent(i-1,pp_BinShift_MB->GetBinContent(i-1)/hBinningAndResol_extendedPt->GetBinContent(i));
    pp_BinShift_MB->SetBinError(i-1,pp_BinShift_MB->GetBinError(i-1)/hBinningAndResol_extendedPt->GetBinContent(i));
  }

  pp_BinShift->Write();
  pp_BinShift_MB->Write();
  RpPbSpectrum->Write("RpPb_Published");
  
  c1->SetLogy(0);
  //pPbSpectrum->Divide(s.pp);
  pPbSpectrum->Divide(pp_BinShift);
  pPbSpectrum->Scale(1.0/0.098);//<T_{pPb}> From HIN-012-017
  pPbSpectrum->Print("All");
  pPbSpectrum->GetYaxis()->SetTitle("RpPb");
  pPbSpectrum->GetYaxis()->SetRangeUser(0.5,1.5);
  pPbSpectrum->Draw();
  pPbSpectrum->Write("RpPb");

  c1->SaveAs("plots/png/pp_RpPb.png");
  c1->SaveAs("plots/png/pp_RpPb.pdf");
  c1->SaveAs("plots/png/pp_RpPb.C");
  
  /*s.pp->Divide(RpPb);
  s.pp->GetYaxis()->SetTitle("data/interp");
  s.pp->Draw(); */
  pp_BinShift->Divide(RpPb);
  pp_BinShift->GetYaxis()->SetTitle("data/interp");
  pp_BinShift->Scale(0.098);
  pp_BinShift->GetYaxis()->SetRangeUser(0.5,1.5);
  pp_BinShift->Draw();
  pp_BinShift->Write("RpPb_ppVsInterpolation");
  
  float TAAUncert = 0.177;
  TBox* bTAA = new TBox(0.15,0.1,0.2,0.2); 
  bTAA->SetFillColor(kBlue-9);
  bTAA->SetFillStyle(1001);
  bTAA->SetLineWidth(0);
  bTAA->DrawBox(200,1-TAAUncert,300,1+TAAUncert); 
 
  c1->SaveAs("plots/png/pp_dataVsInterp.png");
  c1->SaveAs("plots/pdf/pp_dataVsInterp.pdf");
  c1->SaveAs("plots/png/pp_dataVsInterp.C");
  
  pp_BinShift_MB->Divide(RpPb);
  pp_BinShift_MB->GetYaxis()->SetTitle("MB data/interp");
  //take lumi normalization from the previous made plot (should be same at low pt)
  float MBScale = pp_BinShift_MB->GetBinContent(5)/pp_BinShift->GetBinContent(5);
  pp_BinShift_MB->Scale(1/MBScale);
  pp_BinShift_MB->GetYaxis()->SetRangeUser(0.5,1.5);
  pp_BinShift_MB->SetMarkerColor(kRed);
  pp_BinShift_MB->SetLineColor(kRed);
  pp_BinShift_MB->Draw("same");
  TLegend * ppCompareLeg = new TLegend(0.2,0.2,0.5,0.5);
  ppCompareLeg->AddEntry(pp_BinShift,"Trigger Combined","p");
  ppCompareLeg->AddEntry(pp_BinShift_MB,"Only MB","p");
  ppCompareLeg->Draw("same");
  c1->SaveAs("plots/png/pp_dataVsInterp_MBData_PlusTrigger.png");
  c1->SaveAs("plots/pdf/pp_dataVsInterp_MBData_PlusTrigger.pdf");
  c1->SaveAs("plots/png/pp_dataVsInterp_MBData_PlusTrigger.C");
  delete ppCompareLeg;
  pp_BinShift_MB->SetMarkerColor(kBlack);
  pp_BinShift_MB->SetLineColor(kBlack);
  pp_BinShift_MB->Draw();
  pp_BinShift_MB->Write("RpPb_ppVsInterpolation_MBData");
  c1->SaveAs("plots/png/pp_dataVsInterp_MBData.png");
  c1->SaveAs("plots/pdf/pp_dataVsInterp_MBData.pdf");
  c1->SaveAs("plots/png/pp_dataVsInterp_MBData.C");
  outBinF->Close();
  
//************************************************************************************************************
//***********************************************TRACK TRIGGER PLOTS******************************************
//************************************************************************************************************
  TCanvas * c2 = new TCanvas("c2","c2",800,600);
  s.ppMaxtrkByTrigger[4]->SetMarkerColor(kOrange);
  s.ppMaxtrkByTrigger[4]->SetLineColor(kOrange);
  c2->SetLogy();
  s.ppMaxtrk->Scale(100);
  Ymin = 0.00000000001;
  Ymax = 10;
  s.ppMaxtrk->GetYaxis()->SetRangeUser(Ymin,Ymax);
  s.ppMaxtrk->Draw("h");
  for(int i = 0; i<s.nTriggers_trk; i++) s.ppMaxtrkByTrigger[i]->Draw("same");
  TLegend * leg_trk = new TLegend(0.6,0.6,0.9,0.9);
  leg_trk->AddEntry(s.ppMaxtrk,"Leading Trk p_{T}Spectrum (x100)","l");
  leg_trk->AddEntry(s.ppMaxtrkByTrigger[0],"MB (I)","p");
  leg_trk->AddEntry(s.ppMaxtrkByTrigger[1],"Track 18 (II)","p");
  leg_trk->AddEntry(s.ppMaxtrkByTrigger[2],"Track 24 (III)","p");
  leg_trk->AddEntry(s.ppMaxtrkByTrigger[3],"Track 34 (IV)","p");
  leg_trk->AddEntry(s.ppMaxtrkByTrigger[4],"Track 45 (V)","p");
  leg_trk->AddEntry(s.ppMaxtrkByTrigger[5],"Track 53 (VI)","p");
  leg_trk->Draw("same");
  c2->SaveAs("plots/png/ppMaxtrk_FullSpectrum.png"); 
  c2->SaveAs("plots/pdf/ppMaxtrk_FullSpectrum.pdf"); 

  s.ppMaxtrk->GetXaxis()->SetRangeUser(10,90);
  Ymin = 0.0000000001;
  Ymax = 1;
  s.ppMaxtrk->GetYaxis()->SetRangeUser(Ymin,Ymax);
  s.ppMaxtrk->Draw("h");
  for(int i = 0; i<s.nTriggers_trk; i++) s.ppMaxtrkByTrigger[i]->Draw("same");
  leg_trk->Draw("same"); 
  TLine * l_trk[5];
  for(int i = 0; i<5; i++) 
  {
    l_trk[i] = new TLine(s.triggerBins_trk[i+1],Ymin,s.triggerBins_trk[i+1],Ymax); 
    l_trk[i]->SetLineWidth(2);
    l_trk[i]->SetLineStyle(2);
    l_trk[i]->SetLineColor(1);
    l_trk[i]->Draw("same");
  }
  lat->DrawLatex(16,Ymin*3,"I");
  lat->DrawLatex(22,Ymin*3,"II");
  lat->DrawLatex(28,Ymin*3,"III");
  lat->DrawLatex(40,Ymin*3,"IV");
  lat->DrawLatex(49,Ymin*3,"V");
  lat->DrawLatex(57,Ymin*3,"VI");
  c2->SaveAs("plots/png/ppMaxtrk_FullSpectrum_XZoom.png"); 
  c2->SaveAs("plots/pdf/ppMaxtrk_FullSpectrum_XZoom.pdf"); 

  c2->SetLogy(0);
  for(int i = 0; i<s.nTriggers_trk; i++) s.ppMaxtrkByTrigger[i]->Rebin(2);
  for(int i = 0; i<s.nTriggers_trk-1; i++) s.ppMaxtrkByTrigger[s.nTriggers_trk-1-i]->Divide(s.ppMaxtrkByTrigger[s.nTriggers_trk-2-i]);
  s.ppMaxtrkByTrigger[1]->GetYaxis()->SetRangeUser(0,2);
  s.ppMaxtrkByTrigger[1]->GetXaxis()->SetRangeUser(10,90);
  s.ppMaxtrkByTrigger[1]->GetXaxis()->SetTitle("Leading Track p_{T}");
  s.ppMaxtrkByTrigger[1]->Draw();
  s.ppMaxtrkByTrigger[2]->Draw("same");
  s.ppMaxtrkByTrigger[3]->Draw("same");
  s.ppMaxtrkByTrigger[4]->Draw("same");
  s.ppMaxtrkByTrigger[5]->Draw("same");
  TLegend * leg2_trk = new TLegend(0.2,0.6,0.5,0.9);
  leg2_trk->AddEntry(s.ppMaxtrkByTrigger[1],"Trk18/MB","p");
  leg2_trk->AddEntry(s.ppMaxtrkByTrigger[2],"Trk24/Trk18","p");
  leg2_trk->AddEntry(s.ppMaxtrkByTrigger[3],"Trk34/Trk24","p");
  leg2_trk->AddEntry(s.ppMaxtrkByTrigger[4],"Trk45/Trk34","p");
  leg2_trk->AddEntry(s.ppMaxtrkByTrigger[5],"Trk53/Trk45","p");
  leg2_trk->Draw("same");
  c2->SaveAs("plots/png/TrackRelativeTurnOnes.png");
  c2->SaveAs("plots/pdf/TrackRelativeTurnOnes.pdf");
  
  c2->SetLogy();
  c2->SetLogx();
  s.pp_trk->SetMarkerSize(0.8);
  s.pp_trk->GetXaxis()->SetRangeUser(0.7,400);
  s.pp_trk->Draw();
  s.ppUsedByTrigger_trk[0]->SetFillColor(kGray);
  s.ppUsedByTrigger_trk[3]->SetFillColor(kBlue);
  s.ppUsedByTrigger_trk[4]->SetFillColor(kViolet);
  s.ppUsedByTrigger_trk[5]->SetFillColor(kCyan+2);
  s.ppUsedByTrigger_trk[4]->Add(s.ppUsedByTrigger_trk[5]);
  s.ppUsedByTrigger_trk[3]->Add(s.ppUsedByTrigger_trk[4]);
  s.ppUsedByTrigger_trk[2]->Add(s.ppUsedByTrigger_trk[3]);
  s.ppUsedByTrigger_trk[1]->Add(s.ppUsedByTrigger_trk[2]);
  s.ppUsedByTrigger_trk[0]->Add(s.ppUsedByTrigger_trk[1]);
  for(int i = 0; i<s.nTriggers_trk; i++) s.ppUsedByTrigger_trk[i]->Draw("HIST same");
  s.pp_trk->Draw("sameaxis");
  s.pp_trk->Draw("same");
  leg_trk->Clear();
  leg_trk->AddEntry(s.pp_trk,"pp track Spectrum","p");
  leg_trk->AddEntry(s.ppUsedByTrigger_trk[0],"MB trigger","f");
  leg_trk->AddEntry(s.ppUsedByTrigger_trk[1],"Track18 trigger","f");
  leg_trk->AddEntry(s.ppUsedByTrigger_trk[2],"Track24 trigger","f");
  leg_trk->AddEntry(s.ppUsedByTrigger_trk[3],"Track34 trigger","f");
  leg_trk->AddEntry(s.ppUsedByTrigger_trk[4],"Track45 trigger","f");
  leg_trk->AddEntry(s.ppUsedByTrigger_trk[5],"Track53 trigger","f");
  leg_trk->AddEntry((TObject*)0,"|#eta|<1","");
  leg_trk->Draw("same");
   
  c2->SaveAs("plots/png/ppTrack_FullSpectrum_trk.png");
  c2->SaveAs("plots/pdf/ppTrack_FullSpectrum_trk.pdf");
  c2->SetLogy(0);
  for(int i = 0; i<s.nTriggers_trk; i++) s.ppUsedByTrigger_trk[i]->Divide(s.pp_trk);
  s.ppUsedByTrigger_trk[0]->GetYaxis()->SetRangeUser(0,2);
  s.ppUsedByTrigger_trk[0]->GetYaxis()->SetTitle("Relative Contribution to Bin");
  s.ppUsedByTrigger_trk[0]->Draw("HIST");
  for(int i = 1; i<s.nTriggers_trk; i++) s.ppUsedByTrigger_trk[i]->Draw("HIST same");
  s.ppUsedByTrigger_trk[0]->Draw("sameaxis");
  leg_trk->Draw("same");
  c2->SaveAs("plots/png/ppTrack_FullSpectrum_trk_relativeContribution.png");
  c2->SaveAs("plots/pdf/ppTrack_FullSpectrum_trk_relativeContribution.pdf");
  c2->SetLogy();

  c2->SetLogy(0);
  jtVsTrk->Divide(s.pp_trk);
  jtVsTrk->GetYaxis()->SetTitle("pp Jet triggers/pp Track triggers");
  jtVsTrk->GetYaxis()->SetRangeUser(0.5,1.5);
  jtVsTrk->Draw();
  TLegend *leg3 = new TLegend(0.2,0.7,0.4,0.9);
  leg3->AddEntry((TObject*)0,"|#eta|<1 pp","");
  leg3->Draw("same");
  c2->SaveAs("plots/png/pp_jetVsTrkTriggers.png");
  c2->SaveAs("plots/pdf/pp_JetVsTrkTriggers.pdf");

  delete c1;
  delete c2;
  delete lat;
  delete leg;
  delete leg2;
  delete leg3;
  delete leg_trk; 
  delete leg2_trk; 
  delete RpPb;
  delete RpPbSpectrum;
  delete pPbSpectrum;
  return;
}
