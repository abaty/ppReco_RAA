#ifndef PRETTYPLOT
#define PRETTYPLOT
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TAttPad.h"
#include "TGraph.h"
#include "TAttMarker.h"
#include "TLine.h"
#include "TAttLine.h"
#include "TColor.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TAttAxis.h"
#include "TAttText.h"
#include "TCanvas.h"
#include "TAttCanvas.h"
#include "TBox.h"
#include "TAttFill.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "Settings.h"
#include "hyperon_check/hyperonCorrection.C"
#include "TFrame.h"
#include <vector>
#include <string>
#include <fstream>


inline int getHypInd(int c){
  gStyle->SetLegendBorderSize(0);

  int hyperonIndex = 7;
  if(c==0) hyperonIndex=6;
  if(c==1) hyperonIndex=6;
  if(c==20) hyperonIndex=7;
  if(c==23) hyperonIndex=2;
  if(c==24) hyperonIndex=3;
  if(c==25 || c==30) hyperonIndex=4;
  //if(c==30) hyperonIndex=5;
  if(c==31 || c==32) hyperonIndex=6;
  return hyperonIndex;
}
void getSPS(TGraphErrors * SPS);
void getNA49(TGraphErrors * NA49);
void getPHENIX(TGraphErrors * PHENIX);
void getSTAR(TGraphErrors * gSTAR);
void getAlice276(TGraphErrors * Alice276, int centBin = 0);
void getAtlas276(TGraphErrors * Atlas276, int centBin = 0);
void getCMS276(TGraphErrors * CMS276, TBox ** boxes, int centBin = 0);
void get276RAA(TCanvas * c276, Settings s, int centralityBin, bool doAddTheory=false, bool doWithSystBoxes=true);
void gettheoryRAA(TCanvas * c_th, Settings s, int centralityBin, std::string saveString, TGraph * vitev, TGraph * JiechenXu, TGraph * santiago, TGraph * BBMG, TGraph * JiechenXu_2, TGraph * hybrid, TGraph * jetscape);

double Quad(double a, double b)
{
  return TMath::Power(TMath::Power(a,2) + TMath::Power(b,2),0.5);
}

void prettyPlotting(Settings s){
  bool isForDMeson = false;

  TFile * inputPlots = TFile::Open("Spectra.root","Update");
  TH1D * h[s.nCentBins];
  TH1D * RCP[s.nCentBins];
  TH1D * pbpbSpec[s.nCentBins];
  TH1D * ppSpec;

  for(int c = 0; c<s.nCentBins; c++){
    if(s.lowCentBin[c]*5<50) continue; 
    h[c] = (TH1D*)inputPlots->Get(Form("RAA_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
  }
  for(int c = 0; c<s.nCentBins; c++){
    if(s.lowCentBin[c]*5<50) continue; 
    pbpbSpec[c] = (TH1D*)inputPlots->Get(Form("PbPbTrackSpectrum_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
  }
  ppSpec = (TH1D*)inputPlots->Get(Form("pp_NotperMBTrigger"));
  bool hasPP = false;
  for(int c = 0; c<s.nCentBins; c++){
     if(s.lowCentBin[c]*5<50) continue; 
     h[c]->SetDirectory(0);
     h[c]->GetYaxis()->CenterTitle(true);
     h[c]->GetYaxis()->SetTitle("R_{AA}");
     h[c]->GetXaxis()->SetTitle("p_{T} (GeV)");
     h[c]->GetXaxis()->CenterTitle(true);
     s.RAA_totSyst[c] = (TH1D*)h[c]->Clone(Form("RAA_totSyst_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
     s.RAA_totSyst[c]->Reset();
     s.RAA_totSyst[c]->SetDirectory(inputPlots);
     s.PbPb_totSyst[c] = (TH1D*)h[c]->Clone(Form("RAA_totSyst_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
     s.PbPb_totSyst[c]->Reset();
     s.PbPb_totSyst[c]->SetDirectory(inputPlots);
     if(hasPP==false){
       //s.pp_totSyst = (TH1D*)h[c]->Clone(Form("RAA_totSyst_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
       s.pp_totSyst = (TH1D*)h[c]->Clone("ppSpectrum_NormSyst");
       s.pp_totSyst->Reset();
       s.pp_totSyst->SetTitle("ppSpectrum_NormSyst");
       s.pp_totSyst->SetDirectory(inputPlots);
       hasPP=true;
     }
  }
  std::cout << "loaded stuff" << std::endl;
  //accounting for the 99% event selection efficiency by scaling by 1.01 in 0-100%
  h[20]->Scale(0.99);
  //which hyperon correction to pick up (cent dependent) (inclusive is the default)
  TH1D * hyperonPbPb[8];
  for(int i = 0; i<8; i++){
    std::cout << i << std::endl;
    hyperonPbPb[i] = (TH1D*)h[10]->Clone("hyperonPbPb");
    returnHyperonCorrection(0,hyperonPbPb[i],i,"hyperon_check/CurrentHyperonFractions_PASResult/");
    //hyperonPbPb[i]->Print("All");
  }
  TH1D * hyperonpp = (TH1D*)h[10]->Clone("hyperonpp");
  returnHyperonCorrection(1,hyperonpp,0,"hyperon_check/CurrentHyperonFractions_PASResult/");//need to change to pp
  //hyperonpp->Print("All");
  hasPP=false;
  for(int c = 0; c<s.nCentBins; c++){
    if(s.lowCentBin[c]*5<50) continue; 
    h[c]->Multiply(hyperonPbPb[getHypInd(c)]);
    h[c]->Divide(hyperonpp);

    //hyperon correction for spectra
    if(c==0){ ppSpec->Multiply(hyperonpp); hasPP=true;}
    pbpbSpec[c]->Multiply(hyperonPbPb[getHypInd(c)]);

  }
  std::cout << "loaded hyperon stuff" << std::endl;

  //hyperon correctoin plot
/*  TCanvas * c = new TCanvas("c","c",600,600);
  c->SetLogx();
  hyperonPbPb[7]->GetYaxis()->SetRangeUser(0.95,1.2);
  hyperonPbPb[7]->GetYaxis()->SetTitle("1+(Species Correction)");
  hyperonPbPb[7]->GetXaxis()->SetRangeUser(0.7,14);
  hyperonPbPb[7]->GetXaxis()->SetTitle("p_{T} (GeV)");
  hyperonPbPb[7]->Draw("p");
  for(int i = 2; i<7; i++){
    if(i==5) continue; 
    hyperonPbPb[i]->SetLineColor(i+1); 
    if(i==4) hyperonPbPb[i]->SetLineColor(6); 
    hyperonPbPb[i]->Draw("same h");
  }
  TLegend * hypleg = new TLegend(0.2,0.65,0.45,0.9);
  hypleg->AddEntry(hyperonPbPb[7],"0-100% (From Approval)","p"); 
  //hypleg->AddEntry(hyperonPbPb[0],"0-5%","l"); 
  //hypleg->AddEntry(hyperonPbPb[1],"5-10%","l"); 
  hypleg->AddEntry(hyperonPbPb[6],"0-10%","l");
  hypleg->AddEntry(hyperonPbPb[2],"10-30%","l"); 
  hypleg->AddEntry(hyperonPbPb[3],"30-50%","l"); 
  hypleg->AddEntry(hyperonPbPb[4],"50-100%","l"); 
  //hypleg->AddEntry(hyperonPbPb[6],"0-10%","l");
  hypleg->Draw("same"); 

  c->SaveAs("plots/png/HyperonCorrections.png");
  c->SaveAs("plots/pdf/HyperonCorrections.pdf");
  c->SaveAs("plots/png/HyperonCorrections.C");
*/

  std::cout << "opening fake files" << std::endl;
  TFile * fakeFile = TFile::Open("fakeRateUncertainties/Closure_pp.root","read");
  TH1D * fake_pp = (TH1D*)fakeFile->Get("Fake_0");
  fake_pp->SetDirectory(0);
  fakeFile->Close();
/*  fakeFile = TFile::Open("fakeRateUncertainties/Closure_PbPb.root","read");
  TH1D * fake_PbPb = (TH1D*)fakeFile->Get("Fake_0");
  fake_PbPb->SetDirectory(0);
  fakeFile->Close();*/

  setTDRStyle();
  TLine * line1;
  TLatex * tex = new TLatex(0.1,0.1,"cent");
  TLatex * tex2 = new TLatex(0.1,0.1,"cent");
  TBox* bLumi = new TBox(0.1,0.1,0.15,0.2); 
  TBox* bTAA = new TBox(0.15,0.1,0.2,0.2); 
  TBox* b[s.ntrkBins];
  for(int i = 0; i<s.ntrkBins; i++) b[i] = new TBox(0.1,0.1,0.2,0.2); 
  
 
  int W = 800;
  int H = 700;//700
  int H_ref = 700;//700
  int W_ref = 800;
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.15*W_ref;
  float R = 0.04*W_ref;
    
  TCanvas* canv = new TCanvas("RAA","RAA",50,50,W,H);
  canv->SetLogx();
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(1);
  canv->SetTicky(1);
 
  gStyle->SetErrorX(0);
  

  float TAAUncert;
  float TAAUncertlo;
  float RCP_TAAUncert;
  float lumiUncert;//12% for pp lumi
  for(int c = 0; c<s.nCentBins; c++){
    std::cout << c << std::endl;
    if(s.lowCentBin[c]*5<50) continue; 
    //adding up uncertainties
    for(int i = 1; i<s.RAA_totSyst[c]->GetSize()-1; i++){
      s.RAA_totSyst[c]->SetBinContent(i,0);
      s.PbPb_totSyst[c]->SetBinContent(i,0);

      //start of data/MC differences
      /*if(i<6){//low pt part
        s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));
      }else if(c==24){//30-50
        if(i<9) s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.045));
        else s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.05));//5% difference in data/MC (PbPb)
      }else if(c==25){//50-70
        if(i<9) s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.035));
        else s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.03));//5% difference in data/MC (PbPb)
      }else if(c==30){//70-90
        if(i<9) s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.03));
        else s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.02));//5% difference in data/MC (PbPb)
      }else{
        if(i<9) s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.052));
        else{
          s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.05));//5% difference in data/MC (PbPb)
          s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));//4% difference in data/MC (pp)
        }
      }*/

      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.05));//5% difference in data/MC (PbPb)
      //end data/MC differences     
 
      //nonclosure
      //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.025));//5% for nonclosure PbPb
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.025));//5% for nonclosure (PbPb)
      //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));//1% for nonclosure pp 
      
      //finite MC statistics
      if(c==24 || c==25 || c==30){//30-50,50-70,70-90
        if(i>23){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.03));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.03));
        }else if(i>21){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.04));
        }else if(i>19){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.02));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.02));
        }else if(i>14){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.01));
        }
      }else{
        if(!isForDMeson){
        if(i>27 && c!=23){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.04));
        }else if(i>21 && c!=23){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.045));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.045));
        }else if(i>16){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.01));
        }
        }
        else{//added for DMeson comparison, just changes 'i' values
        if(i>13 && c!=23){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.04));
        }else if(i>7 && c!=23){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.045));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.045));
        }else if(i>3){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.01));
        }
        }

        if(i>23 && c==23){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.035));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.035));
        }else if(i>21 && c==23){
          //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));
          s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.04));
        }
      }

      //!this sytematic is largely bullshit since we don't know the data fake rate!
      //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),fake_PbPb->GetBinContent(fake_PbPb->FindBin(s.RAA_totSyst[c]->GetBinCenter(i)))-1));//3% for MC-based fake rate PbPb
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),fake_pp->GetBinContent(fake_pp->FindBin(s.PbPb_totSyst[c]->GetBinCenter(i)))-1));//3% difference in data/MC (PbPb)
      //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),fake_pp->GetBinContent(fake_pp->FindBin(s.RAA_totSyst[c]->GetBinCenter(i)))-1));//for MC-based fake rate pp 
      //if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),fake_pp->GetBinContent(fake_pp->FindBin(s.pp_totSyst->GetBinCenter(i)))-1));//for MC-based faked rate pp
      //s.RCP_totSyst[c]->SetBinContent(i,Quad(s.RCP_totSyst[c]->GetBinContent(i),0.03));//for MC-based fake rate  */
      
      //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));//1% resolution for not unfolding
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.01));//1% resolution for not unfolding
      //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));//1% resolution for not unfolding
      //if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.01));//1% resolution for not unfolding
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),s.h_HInormSyst[c]->GetBinContent(i)));//add in PbPb normalization uncert
      std::cout << s.RAA_totSyst[c]->GetBinContent(i) << std::endl;
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),s.h_HInormSyst[c]->GetBinContent(i)));//add in PbPb normalization uncert
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),s.h_normSyst->GetBinContent(i)));//add in pp normalization uncert
      std::cout << s.RAA_totSyst[c]->GetBinContent(i) << std::endl;
      //if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),s.h_normSyst->GetBinContent(i)));//add in pp normalization uncert
      //s.RCP_totSyst[c]->SetBinContent(i,Quad(s.RCP_totSyst[c]->GetBinContent(i),s.h_HInormSyst[c]->GetBinContent(i)));//add in PbPb normalization uncert
      //s.RCP_totSyst[c]->SetBinContent(i,Quad(s.RCP_totSyst[c]->GetBinContent(i),s.h_HInormSyst[32]->GetBinContent(i)));//add in PbPb normalization uncert
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.03));//pp uncertainty for pileup
      std::cout << s.RAA_totSyst[c]->GetBinContent(i) << std::endl;
      //if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.03));//pp uncertainty for pileup
      
      //s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));//tight selection data/MC
      //s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.04));//tight selection data/MC
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),TMath::Max(hyperonPbPb[getHypInd(c)]->GetBinContent(i)-1,0.015)));//PbPb hyperon study
      std::cout << s.RAA_totSyst[c]->GetBinContent(i) <<std::endl;
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),TMath::Max(hyperonPbPb[getHypInd(c)]->GetBinContent(i)-1,0.015)));//PbPb hyperon study
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),TMath::Max(hyperonpp->GetBinContent(i)-1,0.015)));//pp hyperon study
      std::cout << s.RAA_totSyst[c]->GetBinContent(i)<< std::endl;
      //if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),TMath::Max(hyperonpp->GetBinContent(i)-1,0.015)));//pp hyperon study
      //s.RCP_totSyst[c]->SetBinContent(i,Quad(s.RCP_totSyst[c]->GetBinContent(i),TMath::Max(hyperonPbPb[getHypInd(c)]->GetBinContent(i)-1,0.015)));//PbPb hyperon study
      //s.RCP_totSyst[c]->SetBinContent(i,Quad(s.RCP_totSyst[c]->GetBinContent(i),TMath::Max(hyperonPbPb[getHypInd(32)]->GetBinContent(i)-1,0.015)));//PbPb hyperon study
      
      //if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.023));//pplumi uncertainty for spectrum
      //if(c==21)s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.02));//2% event selection uncertainty for 0-100% RAA
    }
    //s.RAA_totSyst[c]->Print("All");
    TFile * systFile;
    if(c!=10) systFile = TFile::Open("Spectra_syst.root","update");
    else     systFile = TFile::Open("Spectra_syst.root","recreate");
    s.RAA_totSyst[c]->Write();
    //s.PbPb_totSyst[c]->Write();
    systFile->Close();

    //plotting
    canv->Clear();
    h[c]->GetXaxis()->SetRangeUser(0.7,350);
    h[c]->GetXaxis()->SetLabelOffset(-0.005);
    h[c]->GetYaxis()->SetRangeUser(0,1.6);
    h[c]->SetMarkerSize(1.3);
    h[c]->Draw();
 
    TAAUncert = s.TAAuncert[c]/100.0;
    TAAUncertlo = s.TAAuncertlo[c]/100.0;
    lumiUncert = 0.023;//12% for pp lumi
    bLumi->SetFillColor(kGray);
    bTAA->SetFillColor(kBlue-9);
    bLumi->SetLineWidth(0);
    bTAA->SetLineWidth(0);
    bTAA->DrawBox(0.9,1-TAAUncertlo,TMath::Power(10,TMath::Log10(0.9)+(TMath::Log10(1.5)-TMath::Log10(0.9))/2.0),1+TAAUncert);
    bLumi->DrawBox(TMath::Power(10,TMath::Log10(0.9)+(TMath::Log10(1.5)-TMath::Log10(0.9))/2.0),1-lumiUncert,1.5,1+lumiUncert);
  
    line1 = new TLine(h[c]->GetXaxis()->GetBinLowEdge(3),1,h[c]->GetXaxis()->GetBinUpEdge(h[c]->GetSize()-2),1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw("same");
  
    tex2->DrawLatex(0.9,0.1,Form("%d-%d%s",5*s.lowCentBin[c],5*s.highCentBin[c],"%"));
    tex->SetTextFont(42);
    tex->SetTextSize(lumiTextSize*0.08);
    tex->DrawLatex(1.8,1.03,"T_{AA} and lumi. uncertainty");
    tex->DrawLatex(1.8,0.93,"|#eta|<1");
  
    for(int i = 1; i<h[c]->GetSize()-1; i++){
      if(i<3 && !isForDMeson) continue;
      b[i-1]->SetFillColor(kOrange);
      b[i-1]->SetX1(h[c]->GetXaxis()->GetBinLowEdge(i));
      b[i-1]->SetX2(h[c]->GetXaxis()->GetBinUpEdge(i));
      b[i-1]->SetY1((h[c]->GetBinContent(i))*(1-s.RAA_totSyst[c]->GetBinContent(i)));
      b[i-1]->SetY2(h[c]->GetBinContent(i)*(1+s.RAA_totSyst[c]->GetBinContent(i)));
      b[i-1]->Draw("same");
    }
    h[c]->SetMarkerSize(1.3);
    h[c]->Draw("same");
  
    int iPeriod = 0;
    lumi_sqrtS = "27.4 pb^{-1} (5.02 TeV pp) + 404 #mub^{-1} (5.02 TeV PbPb)";
    writeExtraText = false;  
    extraText  = "Preliminary";
    //extraText  = "Unpublished";
    CMS_lumi( canv, iPeriod, 11 );
 
    gStyle->SetPadTickY(1);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();    
    canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
    canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
    canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d.C",5*s.lowCentBin[c],5*s.highCentBin[c]));
    
    /*TLegend * legRaa276;
    TLegend * legRaa276_2;
    TCanvas * canv_th = (TCanvas*)canv->Clone("canv_th");
    canv->cd();
    if(c==0 || c==1 || c==23 || c==24 || c==25 || c==30){
      TGraphErrors * atlas276;
      if(c==0){
        atlas276 = new TGraphErrors("atlas276","atlas276");
        getAtlas276(atlas276,0);
        atlas276->Draw("same pZ");
      }
      TGraphErrors * alice276;
      if(c==0 || c==1){
        alice276 = new TGraphErrors("alice276","alice276");
        if(c==1) getAlice276(alice276,1);
        else getAlice276(alice276);
        alice276->Draw("same pZ");
      }
      TGraphErrors * cms276 = new TGraphErrors("cms276","cms276");
      TBox * bp[27];
      for(int i = 0; i<27; i++) bp[i] = new TBox(0.1,0.1,0.2,0.2);
      getCMS276(cms276,bp,c);   
      cms276->Draw("same pZ");
      for(int i = 0; i<27; i++) bp[i]->Draw("same");
      //TCanvas * canv_276 = (TCanvas*)canv->Clone("canv_276");
      //get276RAA(canv_276,s,c,false);
      //get276RAA(canv,s,c);
      h[c]->SetFillColor(kOrange);
      h[c]->Draw("same");
      
      //if(c==0)  legRaa276 = new TLegend(0.5,0.7,0.9,0.91);
      //else if(c==1) legRaa276 = new TLegend(0.5,0.75,0.9,0.91);
      //else legRaa276 = new TLegend(0.5,0.8,0.9,0.91);
      legRaa276 = new TLegend(0.18,0.7,0.55,0.81);
      legRaa276->SetFillStyle(0);
      legRaa276->AddEntry(h[c],"CMS 5.02 TeV","pf");
      legRaa276->SetTextFont(62);
      gStyle->SetPadTickY(1);
      legRaa276->Draw("same");
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.C",5*s.lowCentBin[c],5*s.highCentBin[c]));
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
      //canv_276->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.C",5*s.lowCentBin[c],5*s.highCentBin[c]));
      //canv_276->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
      //canv_276->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
      if(c==0){
        TCanvas * canv_276x2 = (TCanvas*)canv->Clone("canv_276x2");
        get276RAA(canv_276x2,s,c,true);
        //get276RAA(canv,s,c);
        h[c]->Draw("same");
        canv_276x2->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheoryWith276.C",5*s.lowCentBin[c],5*s.highCentBin[c]));
        canv_276x2->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheoryWith276.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
        canv_276x2->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheoryWith276.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
      }
    }*/
    /*if(c==0 || c==1  || c==24 || c==31){
      const int graphPts = 1952;
      TGraph * vitev = new TGraph(2*graphPts);
      const int graphPts2 = 370;
      TGraph * jx = new TGraph(graphPts2);
      TGraph * jx2 = new TGraph(graphPts2*2);
      const int graphPts3 = 38;
      TGraph * santiago = new TGraph(graphPts3);
      const int graphPts4 = 7;
      TGraph * BBMG = new TGraph(graphPts4);
      const int graphPts5 = 18*2;
      TGraph * hybrid = new TGraph(graphPts5);
      const int graphPts6 =  15*2;
      TGraph * jetscape = new TGraph(graphPts6); 

      gettheoryRAA(canv_th,s,c,"",vitev,jx,santiago,BBMG,jx2,hybrid,jetscape);
      vitev->SetFillStyle(3002);vitev->SetFillColor(kRed);vitev->SetLineWidth(0);
      vitev->Draw("same f");
      jx2->Draw("same f");
      jx->Draw("same");
      if(c==0 || c==31){BBMG->Draw("same"); santiago->Draw("same"); }
      if(c==0)  legRaa276->SetY1NDC(0.7);
      else if(c==1 ) legRaa276->SetY1NDC(0.75);
      else legRaa276->SetY1NDC(0.8);
      legRaa276->AddEntry(vitev,"SCET_{G} 0-10%","F");  
      legRaa276->AddEntry(jx,"CUJET 3.0 (h^{#pm}+#pi^{0}), 0-5%","L");
      if(c==0 || c==31) legRaa276->AddEntry(santiago,"Andr#acute{e}s et al. 0-5%","L");
      legRaa276->Draw("same");
      gStyle->SetPadTickY(1);
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheoryWith276.C",5*s.lowCentBin[c],5*s.highCentBin[c]));
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheoryWith276.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheoryWith276.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
     
      delete legRaa276;
      canv->Clear();
      tex2->DrawLatex(0.9,0.1,Form("%d-%d%s",5*s.lowCentBin[c],5*s.highCentBin[c],"%"));
      tex->SetTextFont(42);
      tex->SetTextSize(lumiTextSize*0.08);
      h[c]->Draw();
      tex->DrawLatex(1.8,1.03,"T_{AA} and lumi. uncertainty");
      tex->DrawLatex(1.8,0.93,"|#eta|<1");
      for(int i = 1; i<h[c]->GetSize()-1; i++){
        if(i<3) continue;
        b[i-1]->Draw("same");
      }
      h[c]->Draw("same");
      bTAA->DrawBox(0.9,1-TAAUncertlo,TMath::Power(10,TMath::Log10(0.9)+(TMath::Log10(1.5)-TMath::Log10(0.9))/2.0),1+TAAUncert);
      bLumi->DrawBox(TMath::Power(10,TMath::Log10(0.9)+(TMath::Log10(1.5)-TMath::Log10(0.9))/2.0),1-lumiUncert,1.5,1+lumiUncert);
      line1->Draw("same");
      if(c==0 || c==31) hybrid->Draw("same f");
      if(c==0 || c==31) jetscape->Draw("same f");
      vitev->Draw("same f");
      jx2->SetFillStyle(3001);jx2->SetFillColor(kCyan+2);jx2->SetLineWidth(3); jx2->SetLineColor(kBlue+1);
      jx2->Draw("same f");
      jx->Draw("same"); 
      if(c==0 || c==31){  BBMG->Draw("same"); santiago->Draw("same"); }
      double yLegbot = 0.7;
      if(c==0 || c==31) yLegbot = 0.67;
      TLegend  *legth = new TLegend(0.5,yLegbot,0.9,0.91);
      //TLegend  *legth = new TLegend(0.5,0.8,0.9,0.91);
      if(c==0 || c==1 || c==31) legth->AddEntry(vitev,"SCET_{G} (0-10%)","F"); 
      if(c==0 || c==1 || c==31) legth->AddEntry(hybrid,"Hybrid Model (0-10%)","F");
      if(c==0 || c==1 || c==31) legth->AddEntry(jetscape,"Bianchi et al. (0-10%)","F");
      if(c==0 || c==1 || c==31) legth->AddEntry(jx2,"CUJET 3.0 (h^{#pm}+#pi^{0}, 0-10%)","LF");
      if(c==24) legth->AddEntry(vitev,"SCET_{G} (30-50%)","F");  
      if(c==24) legth->AddEntry(jx2,"CUJET 3.0 (h^{#pm}+#pi^{0}, 30-50%)","LF");
      if(c==0 || c==31) legth->AddEntry(santiago,"Andr#acute{e}s et al. (0-5%)","L");
      if(c==0 || c==31) legth->AddEntry(BBMG,"v-USPhydro+BBMG (0-5%)","L");
      legth->Draw("same");
      tex2->DrawLatex(0.9,0.1,Form("%d-%d%s",5*s.lowCentBin[c],5*s.highCentBin[c],"%"));
      CMS_lumi( canv, iPeriod, 11 );
 
      gStyle->SetPadTickY(1);
      canv->Update();
      canv->RedrawAxis();
      canv->GetFrame()->Draw();    
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheory.C",5*s.lowCentBin[c],5*s.highCentBin[c]));
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheory.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
      canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheory.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
      delete legth; 
      delete vitev; delete jx; */
      //gettheoryRAA(canv,s,c,"With276");
    //}

    //RCP
    /*  if(c==0 || c==1 || c==23 || c==24 || c==31){
      canv->Clear();
      RCP[c]->GetXaxis()->SetRangeUser(0.7,350);
      RCP[c]->GetXaxis()->SetLabelOffset(-0.005);
      RCP[c]->GetXaxis()->CenterTitle();
      RCP[c]->GetYaxis()->CenterTitle();
      RCP[c]->GetYaxis()->SetRangeUser(0,1.6);
      RCP[c]->SetMarkerSize(1.3);
      RCP[c]->GetYaxis()->SetTitle("R_{CP}");
      RCP[c]->Draw();
   
      TAAUncert = TMath::Power(s.TAAuncert[c]*s.TAAuncert[c]+s.TAAuncert[32]*s.TAAuncert[32],0.5)/100.0;
      bTAA->SetFillColor(kBlue-9);
      bTAA->SetLineWidth(0);
      bTAA->DrawBox(0.9,1-TAAUncert,TMath::Power(10,TMath::Log10(0.9)+(TMath::Log10(1.5)-TMath::Log10(0.9))/2.0),1+TAAUncert);
    
      line1->Draw("same");
    
      tex2->DrawLatex(0.9,0.1,Form("(%d-%d%%)/(50-90%%)",5*s.lowCentBin[c],5*s.highCentBin[c]));
      tex->SetTextFont(42);
      tex->SetTextSize(lumiTextSize*0.08);
      tex->DrawLatex(1.8,1.03,"T_{AA} uncertainty     |#eta|<1");
    
      for(int i = 1; i<h[c]->GetSize()-1; i++){
        if(i<3) continue;
        b[i-1]->SetFillColor(kOrange);
        b[i-1]->SetX1(RCP[c]->GetXaxis()->GetBinLowEdge(i));
        b[i-1]->SetX2(RCP[c]->GetXaxis()->GetBinUpEdge(i));
        b[i-1]->SetY1((RCP[c]->GetBinContent(i))*(1-s.RCP_totSyst[c]->GetBinContent(i)));
        b[i-1]->SetY2(RCP[c]->GetBinContent(i)*(1+s.RCP_totSyst[c]->GetBinContent(i)));
        b[i-1]->Draw("same");
      }
      RCP[c]->SetMarkerSize(1.3);
      line1->Draw("same");
      RCP[c]->Draw("same");
    
      int iPeriod = 0;
      lumi_sqrtS = "404 #mub^{-1} (5.02 TeV PbPb)";
      writeExtraText = true;  
      extraText  = "Preliminary";
      //extraText  = "Unpublished";
      CMS_lumi( canv, iPeriod, 11 );
   
      gStyle->SetPadTickY(1);
      canv->Update();
      canv->RedrawAxis();
      canv->GetFrame()->Draw();    
      canv->SaveAs(Form("plots/prettyPlots/RCP_%d_%d.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
      canv->SaveAs(Form("plots/prettyPlots/RCP_%d_%d.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
      canv->SaveAs(Form("plots/prettyPlots/RCP_%d_%d.C",5*s.lowCentBin[c],5*s.highCentBin[c]));
      }*/
    delete line1;
  }

  /*TCanvas * canv2 = new TCanvas("canv2","canv2",700,800);
  canv2->SetBorderSize(0);
  TPad * pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0,0);
  TPad * pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3,0);
  canv2->SetLineWidth(0);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.15);
  pad1->SetTopMargin(0.08);
  pad1->SetBorderSize(0);
  pad1->Draw();
  pad2->SetTopMargin(0);
  pad2->SetLeftMargin(0.15);
  pad2->SetBottomMargin(0.3);
	  pad2->SetBorderSize(0);
  pad2->Draw();
  pad1->cd();
  pad1->SetLogx();
  pad1->SetLogy();
  //pbpbSpec[0]->GetXaxis()->SetRangeUser(0.7,390);
  TH1D * ppSpecD = new TH1D("specDummy1","",3,0.4,450);
  ppSpecD->GetXaxis()->SetRangeUser(0.4,390);
  ppSpecD->GetYaxis()->SetTitle("#frac{1}{N_{evt}} E #frac{d^{3}N}{dp^{3}} (GeV^{-2})");
  ppSpecD->GetYaxis()->SetTitleOffset(1.4);
  ppSpecD->GetYaxis()->SetTitleSize(0.045);
  ppSpecD->GetYaxis()->SetLabelSize(0.04);
  ppSpecD->GetYaxis()->CenterTitle();
  ppSpecD->GetYaxis()->SetLabelOffset(0.002);
  ppSpecD->GetYaxis()->SetRangeUser(1.1e-17,1e4);
  ppSpecD->GetXaxis()->SetRangeUser(0.4,390);
  ppSpecD->Draw();
  ppSpec->SetMarkerStyle(5);
  ppSpec->Scale(1/70.0);//scaled by inelastic xsection of 70 mb
  //ppSpec->Print("All");
  ppSpec->Draw("same");
  pbpbSpec[0]->GetXaxis()->SetRangeUser(0.7,390);
  pbpbSpec[0]->SetMarkerStyle(24);
  pbpbSpec[0]->Scale(10);
  pbpbSpec[0]->Draw("same");
  pbpbSpec[1]->SetMarkerColor(kBlue);
  pbpbSpec[1]->SetLineColor(kBlue);
  pbpbSpec[1]->SetMarkerStyle(25);
  pbpbSpec[1]->Scale(3);
  pbpbSpec[1]->GetXaxis()->SetRangeUser(0.7,390);
  pbpbSpec[1]->Draw("same");
  pbpbSpec[23]->SetMarkerColor(kRed);
  pbpbSpec[23]->SetLineColor(kRed);
  pbpbSpec[23]->SetMarkerStyle(28);
  pbpbSpec[23]->GetXaxis()->SetRangeUser(0.7,390);
  pbpbSpec[23]->Draw("same");
  pbpbSpec[24]->SetMarkerStyle(20);
  pbpbSpec[24]->GetXaxis()->SetRangeUser(0.7,390);
  pbpbSpec[25]->SetMarkerColor(kBlue);
  pbpbSpec[25]->SetLineColor(kBlue);
  pbpbSpec[24]->Draw("same");
  pbpbSpec[25]->SetMarkerStyle(21);
  pbpbSpec[25]->GetXaxis()->SetRangeUser(0.7,390);
  pbpbSpec[25]->Draw("same");
  pbpbSpec[30]->SetMarkerColor(kRed);
  pbpbSpec[30]->SetLineColor(kRed);
  pbpbSpec[30]->SetMarkerStyle(34);
  pbpbSpec[30]->GetXaxis()->SetRangeUser(0.7,390);
  pbpbSpec[30]->Draw("same");
  TLegend * specLeg = new TLegend(0.25,0.1,0.45,0.5);
  //specLeg->SetFillStyle(0);
  specLeg->AddEntry((TObject*)0,"|#eta|<1",""); 
  specLeg->AddEntry(pbpbSpec[0],Form("0-5%s (x10)","%"),"p");  
  specLeg->AddEntry(pbpbSpec[1],Form("5-10%s (x3)","%"),"p");  
  specLeg->AddEntry(pbpbSpec[23],Form("10-30%s","%"),"p");  
  specLeg->AddEntry(pbpbSpec[24],Form("30-50%s","%"),"p");  
  specLeg->AddEntry(pbpbSpec[25],Form("50-70%s","%"),"p");  
  specLeg->AddEntry(pbpbSpec[30],Form("70-90%s","%"),"p");  
  specLeg->AddEntry(ppSpec,"pp","p"); 
  specLeg->Draw("same"); 
 
  pad2->cd();
  pad2->SetLogx();
  TH1D * ppSpecD2 = new TH1D("specDummy2","",3,0.4,450);
  ppSpecD2->GetYaxis()->SetRangeUser(0.0,19.99);
  ppSpecD2->GetYaxis()->SetNdivisions(4,4,0,kTRUE);
  ppSpecD2->GetYaxis()->SetTitleOffset(0.6);
  ppSpecD2->GetYaxis()->SetTitleFont(42);
  ppSpecD2->GetYaxis()->SetTitleSize(0.095);
  ppSpecD2->GetYaxis()->SetLabelSize(0.095);
  ppSpecD2->GetXaxis()->SetTitleFont(42);
  ppSpecD2->GetYaxis()->SetTitle(Form("Syst. uncert. (%s)","%"));
  ppSpecD2->GetXaxis()->SetRangeUser(0.4,390);
  ppSpecD2->GetXaxis()->SetTitle("p_{T} (GeV)");
  ppSpecD2->GetXaxis()->SetTitleSize(0.1);
  ppSpecD2->GetXaxis()->SetLabelSize(0.1);
  ppSpecD2->GetXaxis()->SetTitleOffset(1.2);
  ppSpecD2->GetXaxis()->CenterTitle();
  ppSpecD2->GetXaxis()->SetTickLength(0.06);
  ppSpecD2->Draw();
  s.PbPb_totSyst[0]->SetFillColor(kOrange);
  s.PbPb_totSyst[0]->SetBinContent(1,0);
  s.PbPb_totSyst[0]->SetBinError(1,0);
  s.PbPb_totSyst[0]->SetBinContent(2,0);
  s.PbPb_totSyst[0]->SetBinError(2,0);
  s.PbPb_totSyst[0]->Scale(100);
  s.PbPb_totSyst[0]->GetXaxis()->SetRangeUser(0.7,400);
  s.PbPb_totSyst[0]->Draw("same");
  s.PbPb_totSyst[30]->SetFillColor(kRed);
  s.PbPb_totSyst[30]->SetFillStyle(3004);
  s.PbPb_totSyst[30]->Scale(100);
  s.PbPb_totSyst[30]->SetBinContent(1,0);
  s.PbPb_totSyst[30]->SetBinError(1,0);
  s.PbPb_totSyst[30]->SetBinContent(2,0);
  s.PbPb_totSyst[30]->SetBinError(2,0);
  s.PbPb_totSyst[30]->GetXaxis()->SetRangeUser(0.7,400);
  s.PbPb_totSyst[30]->Draw("same");
  s.pp_totSyst->SetFillColor(kBlack);
  s.pp_totSyst->SetFillStyle(3003);
  s.pp_totSyst->GetXaxis()->SetRangeUser(0.5,400);
  s.pp_totSyst->Scale(100);
  s.pp_totSyst->Draw("same");
  TLegend * systLeg = new TLegend(0.6,0.6,0.9,0.98);
  //systLeg->SetFillStyle(0);
  systLeg->AddEntry(s.PbPb_totSyst[0],Form("0-5%s","%"),"f");
  systLeg->AddEntry(s.PbPb_totSyst[30],Form("70-90%s","%"),"f");
  systLeg->AddEntry(s.pp_totSyst,"pp","f");
  gStyle->SetPadTickY(1);
  systLeg->Draw("same");
  ppSpecD2->Draw("sameaxis");
  ppSpecD2->GetXaxis()->Draw("same");

  CMS_lumi( canv2, 0,33);
  //canv2->Update();
  //canv2->RedrawAxis();
  //canv2->GetFrame()->Draw();    
  canv2->SaveAs("plots/prettyPlots/Spectra_perEventYield.png");
  canv2->SaveAs("plots/prettyPlots/Spectra_perEventYield.pdf");
  canv2->SaveAs("plots/prettyPlots/Spectra_perEventYield.C");*/

  //Huge results compilation plot
  /*TCanvas * canv3 = new TCanvas("canv3","canv3",1300,1300);
  canv3->SetBorderSize(0);
  canv3->SetLineWidth(0);
  canv3->SetLogx();
  TH1D * axisDummy = new TH1D("axisD","axisD",1,0.16,500);
  axisDummy->GetYaxis()->CenterTitle(true);
  axisDummy->GetYaxis()->SetTitle("R_{AA}");
  axisDummy->GetXaxis()->SetTitle("p_{T} (GeV)");
  axisDummy->GetXaxis()->CenterTitle(true);
  axisDummy->GetXaxis()->SetRangeUser(0.16,490);
  axisDummy->GetXaxis()->SetLabelOffset(99);
  axisDummy->GetYaxis()->SetRangeUser(0,2);
  axisDummy->Draw();
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
 //tex->DrawLatex(0.9082802,-0.09099604,"0.1");
  tex->DrawLatex(0.9082802,-0.09099604,"1");
  tex->DrawLatex(7.914521,-0.09969131,"10");
  tex->DrawLatex(72.06083,-0.09969131,"100");
  h[0]->SetBinContent(1,0); h[0]->SetBinContent(2,0);
  h[0]->SetBinError(1,0); h[0]->SetBinError(2,0);
  h[0]->SetMarkerSize(1.3);
  h[0]->Draw("same");
 
  
  for(int i = 1; i<h[0]->GetSize()-1; i++){
    if(i<3) continue;
    b[i-1]->SetFillColor(kOrange);
    b[i-1]->SetX1(h[0]->GetXaxis()->GetBinLowEdge(i));
    b[i-1]->SetX2(h[0]->GetXaxis()->GetBinUpEdge(i));
    //TAA Uncert is not included for now
    b[i-1]->SetY1(h[0]->GetBinContent(i)*(1-TMath::Power(TMath::Power(s.RAA_totSyst[0]->GetBinContent(i),2)+lumiUncert*lumiUncert,0.5)));
    b[i-1]->SetY2(h[0]->GetBinContent(i)*(1+TMath::Power(TMath::Power(s.RAA_totSyst[0]->GetBinContent(i),2)+lumiUncert*lumiUncert,0.5)));
    b[i-1]->Draw("same");
  }
  
  gStyle->SetErrorX(0);
  TGraphErrors * SPS = new TGraphErrors("SPS","SPS");
  getSPS(SPS);
  SPS->Draw("same p");
  TGraphErrors * NA49 = new TGraphErrors("SPS","SPS");
  getNA49(NA49);
  NA49->Draw("same p");
  TGraphErrors * PHENIX = new TGraphErrors("PHENIX","PHENIX");
  getPHENIX(PHENIX);
  PHENIX->Draw("same p");
  TGraphErrors * STAR = new TGraphErrors("STAR","STAR");
  getSTAR(STAR);
  STAR->Draw("same p");

  gStyle->SetErrorX(0);
  get276RAA(canv3,s,0,false,false);
  h[0]->SetMarkerSize(1.3);
  h[0]->Draw("same");

  //legend
  gStyle->SetLegendBorderSize(0); 
  TLegend * bigLegc1 = new TLegend(0.48,0.55,0.68,0.93);
  TLegend * bigLegc1x = new TLegend(0.44,0.55,0.64,0.93);
 //divide by number of entries on right/left
  TLegend * bigLegc2 = new TLegend(0.69,0.93-(0.93-0.55)/(10.0/2.0),0.89,0.93);
  TLegend * bigLegc2x = new TLegend(0.63,0.93-(0.93-0.55)/(10.0/2.0),0.85,0.93);
  bigLegc1->SetFillStyle(0); bigLegc1x->SetFillStyle(0);
  bigLegc2->SetFillStyle(0); bigLegc2x->SetFillStyle(0);
  bigLegc1->SetTextSize(0.4/20.0); bigLegc1x->SetTextSize(0.4/20.0);
  bigLegc2->SetTextSize(0.4/20.0); bigLegc2x->SetTextSize(0.4/20.0);
  //col1
  TH1D * dummyATLAS = new TH1D("dummyATLAS","dummyATLAS",10,0,10);
  dummyATLAS->SetMarkerColor(kBlue);
  dummyATLAS->SetLineColor(kBlue); 
  dummyATLAS->SetMarkerSize(1.4);
  dummyATLAS->SetMarkerStyle(32);
  TH1D * dummyCMS = new TH1D("dummyCMS","dummyCMS",10,0,10);
  dummyCMS->SetMarkerColor(kRed); 
  dummyCMS->SetLineColor(kRed); 
  dummyCMS->SetMarkerSize(1.4);
  dummyCMS->SetMarkerStyle(24);
  TH1D * dummyALICE = new TH1D("dummyALICE","dummyALICE",10,0,10);
  dummyALICE->SetMarkerStyle(27);
  dummyALICE->SetMarkerColor(kGreen+2);
  dummyALICE->SetLineColor(kGreen+2); 
  dummyALICE->SetMarkerSize(1.7); 
  bigLegc1->AddEntry((TObject*)0,"","");  bigLegc1x->AddEntry((TObject*)0,"SPS 17.3 GeV (PbPb)","");
  bigLegc1->AddEntry(SPS,"#pi^{0} WA98 (0-7%)","p");  bigLegc1x->AddEntry((TObject*)0,"","");
  bigLegc1->AddEntry(NA49,"#pi^{#pm} NA49 (0-5%)","p");  bigLegc1x->AddEntry((TObject*)0,"","");
  bigLegc1->AddEntry((TObject*)0,"","");  bigLegc1x->AddEntry((TObject*)0,"RHIC 200 GeV (AuAu)","");
  bigLegc1->AddEntry(PHENIX,"#pi^{0} PHENIX (0-5%)","p");  bigLegc1x->AddEntry((TObject*)0,"","");
  bigLegc1->AddEntry(STAR,"h^{#pm} STAR (0-5%)","p");  bigLegc1x->AddEntry((TObject*)0,"","");
  bigLegc1->AddEntry((TObject*)0,"","");  bigLegc1x->AddEntry((TObject*)0,"LHC 2.76 TeV (PbPb)","");
  bigLegc1->AddEntry(dummyALICE,"ALICE (0-5%)","p");  bigLegc1x->AddEntry((TObject*)0,"","");
  bigLegc1->AddEntry(dummyATLAS,"ATLAS (0-5%)","p");  bigLegc1x->AddEntry((TObject*)0,"","");
  bigLegc1->AddEntry(dummyCMS,"CMS (0-5%)","p");  bigLegc1x->AddEntry((TObject*)0,"","");
  bigLegc1->Draw("same");bigLegc1x->Draw("same");
  //col 2
  TH1D * dummyCMS5 = new TH1D("dummyCMS5","dummyCMS5",10,0,10);
  dummyCMS5->SetMarkerColor(kBlack); dummyCMS5->SetMarkerStyle(8); dummyCMS5->SetFillColor(kOrange);dummyCMS5->SetMarkerSize(1.3);
  TH1D * dummyVitev = new TH1D("dummyVitev","dummyVitev",10,0,10);
  dummyVitev->SetFillStyle(3002);dummyVitev->SetFillColor(kRed);dummyVitev->SetLineWidth(0);
  TH1D * dummyCUTEP = new TH1D("dummyCUTEP","dummyCUTEP",10,0,10);
  dummyCUTEP->SetLineColor(kBlue+1);dummyCUTEP->SetLineWidth(3);
  bigLegc2->AddEntry((TObject*)0,"","");  bigLegc2x->AddEntry((TObject*)0,"LHC 5.02 TeV (PbPb)","");
  bigLegc2->AddEntry(dummyCMS5,"CMS (0-5%)","pf");  bigLegc2x->AddEntry((TObject*)0,"","");
  bigLegc2->Draw("same");bigLegc2x->Draw("same");
  
  TLine * line3 = new TLine(0.16,1,h[0]->GetXaxis()->GetBinUpEdge(h[0]->GetSize()-2),1);
  line3->SetLineWidth(2);
  line3->SetLineStyle(2);
  line3->Draw("same");
  
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
  tex->DrawLatex(0.6071639,0.7285335,"SPS");
  tex->DrawLatex(4.099646,0.4698491,"RHIC");
  tex->DrawLatex(49.60057,0.1850789,"LHC");

  CMS_lumi( canv3,0, 11,false,true );
  gStyle->SetPadTickY(1);
  canv3->Update();
  canv3->RedrawAxis();
  canv3->SaveAs("plots/prettyPlots/RAA_Compilation_noTheory.png");
  canv3->SaveAs("plots/prettyPlots/RAA_Compilation_noTheory.pdf");
  canv3->SaveAs("plots/prettyPlots/RAA_Compilation_noTheory.C");

  //gettheoryRAA(canv3,s,0,"noSave");
  const int graphPts = 1952;
  TGraph * vitev = new TGraph(2*graphPts);
  const int graphPts2 = 370;
  TGraph * jx = new TGraph(graphPts2);
  TGraph * jx2 = new TGraph(graphPts2*2);
  const int graphPts3 = 38;
  TGraph * santiago = new TGraph(graphPts3);
  const int graphPts4 = 7;
  TGraph * BBMG = new TGraph(graphPts4);
  const int graphPts5 = 18*2;
  TGraph * hybrid = new TGraph(graphPts5);
  const int graphPts6 = 15*2;
  TGraph * jetscape= new TGraph(graphPts6);
  gettheoryRAA(canv3,s,0,"",vitev,jx,santiago,BBMG,jx2,hybrid,jetscape);
  vitev->SetFillStyle(3002);vitev->SetFillColor(kRed);vitev->SetLineWidth(0);
  vitev->Draw("same f");
  hybrid->Draw("same f");
  jetscape->Draw("same f");
  jx->Draw("same");
  santiago->Draw("same");
  BBMG->Draw("same");
  h[0]->Draw("same");
  bigLegc2->SetY1NDC(0.93-(0.93-0.55)/(10.0/9.0));
  bigLegc2x->SetY1NDC(0.93-(0.93-0.55)/(10.0/9.0));
  bigLegc2->AddEntry((TObject*)0,"","");  bigLegc2x->AddEntry((TObject*)0,"Models 5.02 TeV (PbPb)","");
  bigLegc2->AddEntry(dummyVitev,"SCET_{G} (0-10%)","f");  bigLegc2x->AddEntry((TObject*)0,"","");
  bigLegc2->AddEntry(hybrid,"Hybrid Model (0-10%)","f");  bigLegc2x->AddEntry((TObject*)0,"","");
  bigLegc2->AddEntry(jetscape,"Bianchi et al. (0-10%)","F");  bigLegc2x->AddEntry((TObject*)0,"","");
  bigLegc2->AddEntry(dummyCUTEP,"CUJET 3.0 (h^{#pm}+#pi^{0}, 0-5%)","l");  bigLegc2x->AddEntry((TObject*)0,"","");
  bigLegc2->AddEntry(santiago,"Andr#acute{e}s et al. (0-5%)","L"); bigLegc2x->AddEntry((TObject*)0,"","");
  BBMG->SetLineStyle(7);
  bigLegc2->AddEntry(BBMG,"v-USPhydro+BBMG (0-5%)","L"); bigLegc2x->AddEntry((TObject*)0,"","");
  bigLegc2->Draw("same");bigLegc2x->Draw("same");
  gStyle->SetPadTickY(1);
  canv3->Update();
  canv3->RedrawAxis();
  canv3->SaveAs("plots/prettyPlots/RAA_Compilation.png");
  canv3->SaveAs("plots/prettyPlots/RAA_Compilation.pdf");
  canv3->SaveAs("plots/prettyPlots/RAA_Compilation.C");*/
  inputPlots->Close();
  
  return;
}


void get276RAA(TCanvas * c276, Settings s, int centralityBin, bool doAddTheory,bool doWithSystBoxes){
  float raaaxis[28] = {1,1.1,1.2,1.4,1.6,1.8,2,2.2,2.4,3.2,4,4.8,5.6,6.4,7.2,9.6,12,14.4,19.2,24,28.8,35.2,41.6,48,60.8,73.6,86.4,103.6};
    int tempCentralityBin = centralityBin;
    if(centralityBin==23) tempCentralityBin=2;
    if(centralityBin==24) tempCentralityBin=3;
    if(centralityBin==25) tempCentralityBin=4;
    if(centralityBin==30) tempCentralityBin=5;
    float raaval[6][27] = {{3.591e-01, 3.701e-01, 3.887e-01, 4.059e-01, 4.181e-01, 4.204e-01, 4.202e-01, 4.104e-01, 3.622e-01, 2.640e-01, 1.882e-01, 1.487e-01, 1.3604e-01, 1.380e-01, 1.519e-01, 1.820e-01, 2.226e-01, 2.768e-01, 3.678e-01, 4.003e-01, 4.395e-01, 5.270e-01, 5.161e-01, 5.303e-01, 6.400e-01, 5.331e-01, 5.202e-01},
    {3.683e-01, 3.796e-01, 4.004e-01, 4.182e-01, 4.305e-01, 4.362e-01, 4.350e-01, 4.292e-01, 3.844e-01, 2.830e-01, 2.109e-01 , 1.717e-01, 1.616e-01, 1.628e-01 , 1.791e-01, 2.103e-01, 2.458e-01, 2.958e-01, 4.071e-01, 3.688e-01, 4.087e-01, 4.884e-01, 7.80e-01, 5.579e-01, 6.219e-01, 6.075e-01, 5.558e-01},
    { 3.999e-01 , 4.134e-01, 4.344e-01, 4.539e-01 , 4.673e-01, 4.731e-01,  4.748e-01, 4.693e-01, 4.290e-01 , 3.348e-01, 2.647e-01, 2.252e-01,  2.148e-01, 2.153e-01 , 2.358e-01, 2.682e-01, 3.176e-01, 3.950e-01, 4.970e-01, 4.893e-01, 4.702e-01, 6.024e-01, 8.08e-01, 7.215e-01, 7.242e-01, 5.498e-01 , 6.87e-01},
    {4.708e-01, 4.840e-01, 5.070e-01, 5.257e-01, 5.391e-01, 5.440e-01, 5.473e-01, 5.421e-01, 5.169e-01, 4.355e-01, 3.766e-01, 3.447e-01, 3.413e-01, 3.468e-01, 3.613e-01, 4.165e-01, 4.619e-01, 5.451e-01, 6.290e-01, 6.017e-01, 7.022e-01, 7.52e-01, 8.72e-01, 1.118e+00, 8.99e-01, 1.66e+00, 7.50e-01},
    {5.568e-01,  5.675e-01, 5.843e-01, 6.004e-01, 6.076e-01, 6.139e-01, 6.167e-01, 6.176e-01, 6.088e-01, 5.544e-01, 5.236e-01, 5.099e-01, 5.152e-01, 5.115e-01, 5.412e-01, 5.832e-01, 6.419e-01, 6.805e-01, 7.546e-01,  5.90e-01, 9.08e-01, 6.36e-01, 2.110e+00, 7.416e-01, 9.84e-01, 6.98e-01, 8.60e-01 },
    {6.053e-01, 6.102e-01, 6.154e-01, 6.250e-01, 6.297e-01, 6.341e-01, 6.406e-01, 6.456e-01, 6.492e-01, 6.390e-01 , 6.256e-01, 6.057e-01, 6.107e-01, 6.185e-01, 6.418e-01, 6.522e-01, 7.560e-01, 7.219e-01, 6.55e-01, 7.05e-01, 4.37e-01, 1.026e+00, 4.690e-01, 7.209e-01, 7.81e-01, 7.51e-01, 6.36e-01}};
    float raavalstat[6][27] = {{3e-04, 3e-04, 3e-04, 4e-04, 4e-04, 5e-04, 6e-04, 7e-04, 5e-04, 7e-04, 8e-04, 1.0e-03, 1.36e-03, 1.9e-03, 2.0e-03, 4.3e-03, 8.5e-03, 1.23e-02, 2.70e-02, 4.23e-02, 5.59e-02, 7.33e-02, 8.78e-02, 2.90e-02, 5.41e-02, 5.98e-02, 9.06e-02},
    {3e-04, 4e-04, 3e-04, 4e-04, 5e-04, 6e-04, 7e-04, 9e-04, 6e-04, 8e-04, 1.0e-03, 1.2e-03, 1.7e-03, 2.4e-03, 2.4e-03, 5.1e-03, 9.6e-03, 1.37e-02, 3.12e-02 ,  4.19e-02, 5.62e-02, 7.32e-02, 1.68e-01, 3.45e-02, 6.90e-02, 7.30e-02, 9.76e-02},
    {3e-04, 3e-04, 3e-04, 3e-04, 4e-04, 5e-04, 6e-04, 7e-04, 5e-04, 7e-04, 9e-04, 1.2e-03, 1.6e-03, 2.3e-03, 2.4e-03, 5.0e-03, 9.8e-03, 1.42e-02, 3.04e-02, 4.28e-02, 4.97e-02,  6.11e-02, 1.10e-01,  8.08e-02,  5.43e-02 , 5.70e-02, 1.09e-01},
    {4e-04, 5e-04, 4e-04, 5e-04,  7e-04, 8e-04, 1.0e-03, 1.1e-03, 8e-04, 1.1e-03, 1.6e-03, 2.2e-03 , 3.1e-03, 4.3e-03, 4.1e-03, 8.8e-03, 1.62e-02, 2.25e-02, 4.54e-02, 6.47e-02, 9.08e-02, 1.10e-01, 1.75e-01, 2.16e-01, 7.2e-02, 1.01e+00, 1.30e-01},
    {8e-04, 9e-04, 8e-04 , 1.0e-03, 1.2e-03, 1.5e-03, 1.7e-03 , 2.0e-03, 1.4e-03, 2.0e-03, 2.9e-03, 4.3e-03, 6.1e-03, 8.3e-03, 8.1e-03, 1.65e-02, 3.02e-02, 3.96e-02, 7.98e-02, 1.03e-01, 1.84e-01, 1.79e-01, 7.38e-01, 4.25e-02, 9.7e-02, 1.12e-01, 1.93e-01},
    {1.7e-03, 1.9e-03, 1.5e-03, 1.8e-03, 2.1e-03, 2.5e-03, 2.9e-03, 3.3e-03, 2.3e-03,  3.8e-03,  6.0e-03, 8.9e-03, 1.30e-02, 1.83e-02 , 1.81e-02, 3.63e-02,  6.86e-02, 8.55e-02, 1.56e-01, 2.66e-01, 2.50e-01, 7.14e-01, 5.06e-02, 7.99e-02, 1.49e-01, 2.08e-01, 2.77e-01}};
    float raavalsyst[6][27] = {{2.62e-02, 2.70e-02, 2.84e-02, 2.96e-02, 3.05e-02, 3.07e-02, 3.07e-02, 3.00e-02, 2.65e-02, 1.93e-02, 1.38e-02, 1.09e-02, 1.000e-02, 1.02e-02, 1.12e-02, 1.37e-02, 1.70e-02, 2.22e-02, 3.12e-02, 3.53e-02, 4.06e-02, 5.35e-02, 5.72e-02, 6.10e-02, 7.73e-02, 6.47e-02, 6.36e-02},
    {2.69e-02, 2.77e-02, 2.92e-02 , 3.05e-02,  3.14e-02, 3.19e-02, 3.18e-02,  3.14e-02  , 2.81e-02   , 2.07e-02, 1.55e-02, 1.26e-02, 1.19e-02, 1.20e-02, 1.33e-02, 1.58e-02, 1.88e-02, 2.37e-02,  3.45e-02,  3.26e-02,  3.77e-02, 4.96e-02, 8.6e-02, 6.42e-02, 7.51e-02 ,7.38e-02,  6.79e-02 },
    {2.84e-02,  2.93e-02,  3.08e-02, 3.22e-02, 3.32e-02, 3.36e-02, 3.37e-02, 3.33e-02, 3.05e-02,  2.38e-02 , 1.88e-02, 1.60e-02, 1.53e-02, 1.54e-02, 1.69e-02, 1.93e-02, 2.30e-02, 2.93e-02, 3.80e-02, 3.87e-02, 3.86e-02, 5.54e-02, 8.2e-02, 7.68e-02, 8.15e-02,  6.23e-02, 7.8e-02},
    {3.34e-02,  3.44e-02, 3.60e-02, 3.73e-02, 3.83e-02, 3.86e-02, 3.89e-02, 3.85e-02, 3.67e-02, 3.10e-02, 2.68e-02, 2.45e-02, 2.43e-02, 2.47e-02, 2.58e-02, 3.00e-02, 3.34e-02, 4.04e-02, 4.81e-02 , 4.76e-02, 5.76e-02 , 6.9e-02, 8.9e-02, 1.19e-01, 1.01e-01 , 1.9e-01, 8.6e-02},
    {3.95e-02, 4.03e-02, 4.15e-02, 4.26e-02, 4.31e-02 , 4.36e-02 , 4.38e-02, 4.39e-02,  4.33e-02, 3.94e-02, 3.72e-02, 3.63e-02, 3.67e-02, 3.65e-02, 3.87e-02 , 4.20e-02, 4.65e-02, 5.04e-02, 5.77e-02, 4.7e-02,  7.5e-02, 5.8e-02, 2.15e-01, 7.89e-02, 1.11e-01, 7.9e-02, 9.8e-02},
    {4.30e-02, 4.33e-02, 4.37e-02, 4.44e-02, 4.47e-02, 4.50e-02, 4.55e-02, 4.59e-02, 4.61e-02, 4.54e-02, 4.45e-02, 4.31e-02, 4.35e-02, 4.41e-02, 4.59e-02, 4.69e-02, 5.47e-02 , 5.35e-02, 5.0e-02, 5.6e-02, 3.6e-02, 9.4e-02, 4.78e-02, 7.67e-02, 8.8e-02, 8.5e-02,  7.3e-02 }};
    
    //turning off some points that are outliars
    raaval[4][22]=0; raavalstat[4][22]=0; raavalsyst[4][22]=0; 
    raaval[3][25]=0; raavalstat[3][25]=0; raavalsyst[3][25]=0; 
//ATLAS Data below
double p8800_d40x1y1_xval[] = { 0.5365, 0.615, 0.7050000000000001, 0.808, 0.9259999999999999, 1.0594999999999999, 1.21, 1.385, 1.5899999999999999, 
    1.825, 2.095, 2.4050000000000002, 2.755, 3.1550000000000002, 3.62, 4.15, 4.755, 5.455, 6.255, 
    7.17, 8.219999999999999, 9.39, 10.75, 12.35, 14.149999999999999, 16.2, 18.6, 21.35, 24.450000000000003, 
    28.05, 33.85, 42.6, 53.65, 67.55, 85.05, 106.9, 134.5 };
  double p8800_d40x1y1_xerrminus[] = { 0.03649999999999998, 0.04200000000000004, 0.04800000000000004, 0.05500000000000005, 0.06299999999999994, 0.0704999999999999, 0.08000000000000007, 0.09499999999999997, 0.10999999999999988, 
    0.125, 0.14500000000000024, 0.16500000000000004, 0.18500000000000005, 0.2150000000000003, 0.25, 0.28000000000000025, 0.3250000000000002, 0.375, 0.4249999999999998, 
    0.4900000000000002, 0.5599999999999987, 0.6100000000000012, 0.75, 0.8499999999999996, 0.9499999999999993, 1.0999999999999996, 1.3000000000000007, 1.4500000000000028, 1.6500000000000021, 
    1.9499999999999993, 3.8500000000000014, 4.899999999999999, 6.149999999999999, 7.75, 9.75, 12.100000000000009, 15.5 };
  double p8800_d40x1y1_xerrplus[] = { 0.03649999999999998, 0.04200000000000004, 0.04799999999999993, 0.05499999999999994, 0.06300000000000006, 0.07050000000000001, 0.08000000000000007, 0.09499999999999997, 0.1100000000000001, 
    0.125, 0.14500000000000002, 0.1649999999999996, 0.18500000000000005, 0.21499999999999986, 0.25, 0.27999999999999936, 0.3250000000000002, 0.375, 0.4249999999999998, 
    0.4900000000000002, 0.5600000000000005, 0.6099999999999994, 0.75, 0.8499999999999996, 0.9500000000000011, 1.1000000000000014, 1.2999999999999972, 1.4499999999999993, 1.6499999999999986, 
    1.9499999999999993, 3.8500000000000014, 4.899999999999999, 6.149999999999999, 7.75, 9.75, 12.099999999999994, 15.5 };
  double p8800_d40x1y1_yval[] = { 0.271, 0.289, 0.311, 0.336, 0.361, 0.388, 0.413, 0.434, 0.452, 
    0.452, 0.451, 0.427, 0.388, 0.336, 0.277, 0.221, 0.18, 0.154, 0.142, 
    0.141, 0.15, 0.164, 0.179, 0.201, 0.225, 0.245, 0.276, 0.322, 0.348, 
    0.351, 0.43, 0.51, 0.545, 0.562, 0.609, 0.614, 0.645 };
  double p8800_d40x1y1_yerrminus[] = { 0.013550227909083301, 0.013872270934854176, 0.014306324856094243, 0.01545642218553828, 0.01696752574610983, 0.01823669384637468, 0.019411929667099048, 0.019965226955774882, 0.02124573099330781, 
    0.02079447651772941, 0.02074962533106562, 0.019647615175639, 0.017467284258292703, 0.015130785486550263, 0.012479887553980606, 0.00996562982505371, 0.008309705409940835, 0.006824980471766933, 0.006475822777068562, 
    0.0063590999999999995, 0.006738323233564861, 0.007493920736170086, 0.00786989536143906, 0.00896423945463306, 0.009902556488099425, 0.011100960544025009, 0.012590562179664576, 0.016180299750004633, 0.017744587907302888, 
    0.017897558492710672, 0.022757524030526697, 0.029549208111216793, 0.039004641390480696, 0.045396939324143876, 0.057652714784301354, 0.0751115513087035, 0.11254279685968356 };
  double p8800_d40x1y1_yerrplus[] = { 0.013550227909083301, 0.013872270934854176, 0.014306324856094243, 0.01545642218553828, 0.01696752574610983, 0.01823669384637468, 0.019411929667099048, 0.019965226955774882, 0.02124573099330781, 
    0.02079447651772941, 0.02074962533106562, 0.019647615175639, 0.017467284258292703, 0.015130785486550263, 0.012479887553980606, 0.00996562982505371, 0.008309705409940835, 0.006824980471766933, 0.006475822777068562, 
    0.0063590999999999995, 0.006738323233564861, 0.007493920736170086, 0.00786989536143906, 0.00896423945463306, 0.009902556488099425, 0.011100960544025009, 0.012590562179664576, 0.016180299750004633, 0.017744587907302888, 
    0.017897558492710672, 0.022757524030526697, 0.029549208111216793, 0.039004641390480696, 0.045396939324143876, 0.057652714784301354, 0.0751115513087035, 0.11254279685968356 };
  double p8800_d40x1y1_ystatminus[] = { 7.859E-5, 8.669999999999999E-5, 9.641E-5, 1.1424000000000002E-4, 1.3356999999999998E-4, 1.5908E-4, 1.8998E-4, 2.2134E-4, 2.712E-4, 
    3.2091999999999993E-4, 3.8785999999999996E-4, 4.696999999999999E-4, 5.044E-4, 5.712000000000001E-4, 6.094000000000001E-4, 6.409E-4, 7.02E-4, 8.162E-4, 0.0010508, 
    0.0013958999999999998, 0.0019500000000000001, 0.002952, 0.002327, 0.0030150000000000003, 0.0036000000000000008, 0.004164999999999999, 0.005520000000000001, 0.008372000000000001, 0.007656, 
    0.007722, 0.0086, 0.01071, 0.010355, 0.016860000000000003, 0.029841000000000003, 0.052804, 0.094815 };
  double p8800_d40x1y1_ystatplus[] = { 7.859E-5, 8.669999999999999E-5, 9.641E-5, 1.1424000000000002E-4, 1.3356999999999998E-4, 1.5908E-4, 1.8998E-4, 2.2134E-4, 2.712E-4, 
    3.2091999999999993E-4, 3.8785999999999996E-4, 4.696999999999999E-4, 5.044E-4, 5.712000000000001E-4, 6.094000000000001E-4, 6.409E-4, 7.02E-4, 8.162E-4, 0.0010508, 
    0.0013958999999999998, 0.0019500000000000001, 0.002952, 0.002327, 0.0030150000000000003, 0.0036000000000000008, 0.004164999999999999, 0.005520000000000001, 0.008372000000000001, 0.007656, 
    0.007722, 0.0086, 0.01071, 0.010355, 0.016860000000000003, 0.029841000000000003, 0.052804, 0.094815 };
  int p8800_d40x1y1_numpoints = 37;
  //TGraphAsymmErrors p8800_d40x1y1 = TGraphAsymmErrors(p8800_d40x1y1_numpoints, p8800_d40x1y1_xval, p8800_d40x1y1_yval, p8800_d40x1y1_xerrminus, p8800_d40x1y1_xerrplus, p8800_d40x1y1_yerrminus, p8800_d40x1y1_yerrplus);
  TGraphAsymmErrors * p8800_d40x1y1 =new  TGraphAsymmErrors(p8800_d40x1y1_numpoints, p8800_d40x1y1_xval, p8800_d40x1y1_yval,0,0, p8800_d40x1y1_ystatminus, p8800_d40x1y1_ystatplus);
  p8800_d40x1y1->SetName("ATLAS_0_5");
  p8800_d40x1y1->SetTitle("ATLAS_0_5");
  p8800_d40x1y1->SetMarkerColor(kBlue);
  p8800_d40x1y1->SetMarkerStyle(32);
  p8800_d40x1y1->SetLineColor(kBlue);
//end ATLAS data

  TH1D * p;
  p = new TH1D("raa276",";p_{T};R_{AA}",27,raaaxis);
  for(int i = 1; i<p->GetSize()-1; i++){
    p->SetBinContent(i,raaval[tempCentralityBin][i-1]);
    p->SetBinError(i,raavalstat[tempCentralityBin][i-1]); 
  }
  p->SetMarkerColor(kRed);
  p->SetMarkerStyle(24);
  p->SetMarkerSize(1.3);
  p->SetLineColor(kRed);
  p->Draw("same");
 
  TBox *bp[27];
  if(doWithSystBoxes){
    for(int i = 0; i<27; i++) bp[i] = new TBox(0.1,0.1,0.2,0.2);
    for(int i = 1; i<p->GetSize()-1; i++){
      if(tempCentralityBin==3 && i-1==25) continue;
      if(tempCentralityBin==4 && i-1==22) continue;
      bp[i-1]->SetFillStyle(0);
      bp[i-1]->SetLineColor(kRed);
      bp[i-1]->SetLineWidth(1);
      bp[i-1]->SetX1(p->GetXaxis()->GetBinLowEdge(i));
      bp[i-1]->SetX2(p->GetXaxis()->GetBinUpEdge(i));
      bp[i-1]->SetY1(p->GetBinContent(i)-raavalsyst[tempCentralityBin][i-1]);
      bp[i-1]->SetY2(p->GetBinContent(i)+raavalsyst[tempCentralityBin][i-1]);
      bp[i-1]->Draw("same");
    }
  }else{
    //adds CMS systematics to the stat error bars, commented out for now
    //for(int i = 1; i<p->GetSize()-1; i++){
      //p->SetBinError(i,TMath::Power(raavalstat[tempCentralityBin][i-1]*raavalstat[tempCentralityBin][i-1]+raavalsyst[tempCentralityBin][i-1]*raavalsyst[tempCentralityBin][i-1],0.5));
    //}
  }

  TGraphErrors * alice276 = new TGraphErrors("alice276","alice276");
  if(centralityBin==1) getAlice276(alice276,1);
  else getAlice276(alice276);

  TLegend * legRaa276;
  /*if(centralityBin==0)  legRaa276 = new TLegend(0.5,0.7,0.9,0.91);
  else if(centralityBin==1) legRaa276 = new TLegend(0.5,0.75,0.9,0.91);
  else legRaa276 = new TLegend(0.5,0.8,0.9,0.91);*/
  if(centralityBin==0)  legRaa276 = new TLegend(0.35,0.7,0.75,0.81);
  else if(centralityBin==1) legRaa276 = new TLegend(0.35,0.7,0.75,0.81);
  else legRaa276 = new TLegend(0.35,0.7,0.75,0.81);
  legRaa276->SetFillStyle(0);
  TH1D * dummy = new TH1D("dummy","dummy",10,0,10);
  dummy->SetMarkerColor(kBlack); dummy->SetMarkerStyle(8); dummy->SetFillColor(kOrange); dummy->SetMarkerSize(1.3);
  legRaa276->AddEntry(dummy,"CMS 5.02 TeV","pf");
  legRaa276->AddEntry(p,"CMS 2.76 TeV","p");

  gStyle->SetErrorX(0.);
  if(centralityBin==0 || centralityBin==1){
    if(centralityBin==0){
      legRaa276->AddEntry(p8800_d40x1y1,"ATLAS 2.76 TeV","p");
      p8800_d40x1y1->SetMarkerSize(1.4);
      p8800_d40x1y1->Draw("P same");
    }

    legRaa276->AddEntry(alice276,"ALICE 2.76 TeV","p");
    alice276->Draw("P same");
    legRaa276->SetTextFont(62);
    legRaa276->Draw("same");
  }
  else legRaa276->Draw("same");
  
  std::cout << doAddTheory << std::endl;
  if(centralityBin==0 && doAddTheory==1){
    TH1D * dummy2 = new TH1D("dummy2","dummy2",10,0,10);
    dummy2->SetFillStyle(3002);dummy2->SetFillColor(kRed);dummy2->SetLineWidth(0);
    TH1D * dummy3 = new TH1D("dummy3","dummy3",10,0,10);
    dummy3->SetLineColor(kBlue+1);dummy3->SetLineWidth(3);
    legRaa276->AddEntry(dummy2,"SCET_{G} 0-10%","f");
    legRaa276->AddEntry(dummy3,"CUJET 3.0 0-10% (h^{#pm}+#pi^{0})","l");
    legRaa276->Draw("same");
    //gettheoryRAA(c276,s,centralityBin,"With276");
  }
  if(!doWithSystBoxes) delete legRaa276;
  if(!doWithSystBoxes) delete dummy;
  //if(doWithSystBoxes){for(int i = 0; i<27; i++) delete bp[i];}
  return;
}


void getCMS276(TGraphErrors * CMS276, TBox ** boxes, int centBin ){
        CMS276->SetMarkerColor(kRed);
        CMS276->SetLineColor(kRed);
        CMS276->SetMarkerStyle(24);
        CMS276->SetMarkerSize(1.4);
  float raaaxis[28] = {1,1.1,1.2,1.4,1.6,1.8,2,2.2,2.4,3.2,4,4.8,5.6,6.4,7.2,9.6,12,14.4,19.2,24,28.8,35.2,41.6,48,60.8,73.6,86.4,103.6};
    int tempCentralityBin = centBin;
    if(centBin==23) tempCentralityBin=2;
    if(centBin==24) tempCentralityBin=3;
    if(centBin==25) tempCentralityBin=4;
    if(centBin==30) tempCentralityBin=5;
    float raaval[6][27] = {{3.591e-01, 3.701e-01, 3.887e-01, 4.059e-01, 4.181e-01, 4.204e-01, 4.202e-01, 4.104e-01, 3.622e-01, 2.640e-01, 1.882e-01, 1.487e-01, 1.3604e-01, 1.380e-01, 1.519e-01, 1.820e-01, 2.226e-01, 2.768e-01, 3.678e-01, 4.003e-01, 4.395e-01, 5.270e-01, 5.161e-01, 5.303e-01, 6.400e-01, 5.331e-01, 5.202e-01},
    {3.683e-01, 3.796e-01, 4.004e-01, 4.182e-01, 4.305e-01, 4.362e-01, 4.350e-01, 4.292e-01, 3.844e-01, 2.830e-01, 2.109e-01 , 1.717e-01, 1.616e-01, 1.628e-01 , 1.791e-01, 2.103e-01, 2.458e-01, 2.958e-01, 4.071e-01, 3.688e-01, 4.087e-01, 4.884e-01, 7.80e-01, 5.579e-01, 6.219e-01, 6.075e-01, 5.558e-01},
    { 3.999e-01 , 4.134e-01, 4.344e-01, 4.539e-01 , 4.673e-01, 4.731e-01,  4.748e-01, 4.693e-01, 4.290e-01 , 3.348e-01, 2.647e-01, 2.252e-01,  2.148e-01, 2.153e-01 , 2.358e-01, 2.682e-01, 3.176e-01, 3.950e-01, 4.970e-01, 4.893e-01, 4.702e-01, 6.024e-01, 8.08e-01, 7.215e-01, 7.242e-01, 5.498e-01 , 6.87e-01},
    {4.708e-01, 4.840e-01, 5.070e-01, 5.257e-01, 5.391e-01, 5.440e-01, 5.473e-01, 5.421e-01, 5.169e-01, 4.355e-01, 3.766e-01, 3.447e-01, 3.413e-01, 3.468e-01, 3.613e-01, 4.165e-01, 4.619e-01, 5.451e-01, 6.290e-01, 6.017e-01, 7.022e-01, 7.52e-01, 8.72e-01, 1.118e+00, 8.99e-01, 1.66e+00, 7.50e-01},
    {5.568e-01,  5.675e-01, 5.843e-01, 6.004e-01, 6.076e-01, 6.139e-01, 6.167e-01, 6.176e-01, 6.088e-01, 5.544e-01, 5.236e-01, 5.099e-01, 5.152e-01, 5.115e-01, 5.412e-01, 5.832e-01, 6.419e-01, 6.805e-01, 7.546e-01,  5.90e-01, 9.08e-01, 6.36e-01, 2.110e+00, 7.416e-01, 9.84e-01, 6.98e-01, 8.60e-01 },
    {6.053e-01, 6.102e-01, 6.154e-01, 6.250e-01, 6.297e-01, 6.341e-01, 6.406e-01, 6.456e-01, 6.492e-01, 6.390e-01 , 6.256e-01, 6.057e-01, 6.107e-01, 6.185e-01, 6.418e-01, 6.522e-01, 7.560e-01, 7.219e-01, 6.55e-01, 7.05e-01, 4.37e-01, 1.026e+00, 4.690e-01, 7.209e-01, 7.81e-01, 7.51e-01, 6.36e-01}};
    float raavalstat[6][27] = {{3e-04, 3e-04, 3e-04, 4e-04, 4e-04, 5e-04, 6e-04, 7e-04, 5e-04, 7e-04, 8e-04, 1.0e-03, 1.36e-03, 1.9e-03, 2.0e-03, 4.3e-03, 8.5e-03, 1.23e-02, 2.70e-02, 4.23e-02, 5.59e-02, 7.33e-02, 8.78e-02, 2.90e-02, 5.41e-02, 5.98e-02, 9.06e-02},
    {3e-04, 4e-04, 3e-04, 4e-04, 5e-04, 6e-04, 7e-04, 9e-04, 6e-04, 8e-04, 1.0e-03, 1.2e-03, 1.7e-03, 2.4e-03, 2.4e-03, 5.1e-03, 9.6e-03, 1.37e-02, 3.12e-02 ,  4.19e-02, 5.62e-02, 7.32e-02, 1.68e-01, 3.45e-02, 6.90e-02, 7.30e-02, 9.76e-02},
    {3e-04, 3e-04, 3e-04, 3e-04, 4e-04, 5e-04, 6e-04, 7e-04, 5e-04, 7e-04, 9e-04, 1.2e-03, 1.6e-03, 2.3e-03, 2.4e-03, 5.0e-03, 9.8e-03, 1.42e-02, 3.04e-02, 4.28e-02, 4.97e-02,  6.11e-02, 1.10e-01,  8.08e-02,  5.43e-02 , 5.70e-02, 1.09e-01},
    {4e-04, 5e-04, 4e-04, 5e-04,  7e-04, 8e-04, 1.0e-03, 1.1e-03, 8e-04, 1.1e-03, 1.6e-03, 2.2e-03 , 3.1e-03, 4.3e-03, 4.1e-03, 8.8e-03, 1.62e-02, 2.25e-02, 4.54e-02, 6.47e-02, 9.08e-02, 1.10e-01, 1.75e-01, 2.16e-01, 7.2e-02, 1.01e+00, 1.30e-01},
    {8e-04, 9e-04, 8e-04 , 1.0e-03, 1.2e-03, 1.5e-03, 1.7e-03 , 2.0e-03, 1.4e-03, 2.0e-03, 2.9e-03, 4.3e-03, 6.1e-03, 8.3e-03, 8.1e-03, 1.65e-02, 3.02e-02, 3.96e-02, 7.98e-02, 1.03e-01, 1.84e-01, 1.79e-01, 7.38e-01, 4.25e-02, 9.7e-02, 1.12e-01, 1.93e-01},
    {1.7e-03, 1.9e-03, 1.5e-03, 1.8e-03, 2.1e-03, 2.5e-03, 2.9e-03, 3.3e-03, 2.3e-03,  3.8e-03,  6.0e-03, 8.9e-03, 1.30e-02, 1.83e-02 , 1.81e-02, 3.63e-02,  6.86e-02, 8.55e-02, 1.56e-01, 2.66e-01, 2.50e-01, 7.14e-01, 5.06e-02, 7.99e-02, 1.49e-01, 2.08e-01, 2.77e-01}};
    float raavalsyst[6][27] = {{2.62e-02, 2.70e-02, 2.84e-02, 2.96e-02, 3.05e-02, 3.07e-02, 3.07e-02, 3.00e-02, 2.65e-02, 1.93e-02, 1.38e-02, 1.09e-02, 1.000e-02, 1.02e-02, 1.12e-02, 1.37e-02, 1.70e-02, 2.22e-02, 3.12e-02, 3.53e-02, 4.06e-02, 5.35e-02, 5.72e-02, 6.10e-02, 7.73e-02, 6.47e-02, 6.36e-02},
    {2.69e-02, 2.77e-02, 2.92e-02 , 3.05e-02,  3.14e-02, 3.19e-02, 3.18e-02,  3.14e-02  , 2.81e-02   , 2.07e-02, 1.55e-02, 1.26e-02, 1.19e-02, 1.20e-02, 1.33e-02, 1.58e-02, 1.88e-02, 2.37e-02,  3.45e-02,  3.26e-02,  3.77e-02, 4.96e-02, 8.6e-02, 6.42e-02, 7.51e-02 ,7.38e-02,  6.79e-02 },
    {2.84e-02,  2.93e-02,  3.08e-02, 3.22e-02, 3.32e-02, 3.36e-02, 3.37e-02, 3.33e-02, 3.05e-02,  2.38e-02 , 1.88e-02, 1.60e-02, 1.53e-02, 1.54e-02, 1.69e-02, 1.93e-02, 2.30e-02, 2.93e-02, 3.80e-02, 3.87e-02, 3.86e-02, 5.54e-02, 8.2e-02, 7.68e-02, 8.15e-02,  6.23e-02, 7.8e-02},
    {3.34e-02,  3.44e-02, 3.60e-02, 3.73e-02, 3.83e-02, 3.86e-02, 3.89e-02, 3.85e-02, 3.67e-02, 3.10e-02, 2.68e-02, 2.45e-02, 2.43e-02, 2.47e-02, 2.58e-02, 3.00e-02, 3.34e-02, 4.04e-02, 4.81e-02 , 4.76e-02, 5.76e-02 , 6.9e-02, 8.9e-02, 1.19e-01, 1.01e-01 , 1.9e-01, 8.6e-02},
    {3.95e-02, 4.03e-02, 4.15e-02, 4.26e-02, 4.31e-02 , 4.36e-02 , 4.38e-02, 4.39e-02,  4.33e-02, 3.94e-02, 3.72e-02, 3.63e-02, 3.67e-02, 3.65e-02, 3.87e-02 , 4.20e-02, 4.65e-02, 5.04e-02, 5.77e-02, 4.7e-02,  7.5e-02, 5.8e-02, 2.15e-01, 7.89e-02, 1.11e-01, 7.9e-02, 9.8e-02},
    {4.30e-02, 4.33e-02, 4.37e-02, 4.44e-02, 4.47e-02, 4.50e-02, 4.55e-02, 4.59e-02, 4.61e-02, 4.54e-02, 4.45e-02, 4.31e-02, 4.35e-02, 4.41e-02, 4.59e-02, 4.69e-02, 5.47e-02 , 5.35e-02, 5.0e-02, 5.6e-02, 3.6e-02, 9.4e-02, 4.78e-02, 7.67e-02, 8.8e-02, 8.5e-02,  7.3e-02 }};
    
    //turning off some points that are outliars
    //raaval[4][22]=0; raavalstat[4][22]=0; raavalsyst[4][22]=0; 
    //raaval[3][25]=0; raavalstat[3][25]=0; raavalsyst[3][25]=0; 
        for(int i = 0; i<27; i++){
            CMS276->SetPoint(i,(raaaxis[i]+raaaxis[i+1])/2.0,raaval[tempCentralityBin][i]);
            CMS276->SetPointError(i,0,raavalstat[tempCentralityBin][i]);
     }
     //if(tempCentralityBin==4)CMS276->RemovePoint(22);
     //if(tempCentralityBin==3)CMS276->RemovePoint(25);
     for(int i = 1; i<28; i++){
      if(tempCentralityBin==3 && i-1==25) continue;
      if(tempCentralityBin==4 && i-1==22) continue;
      boxes[i-1]->SetFillStyle(0);
      boxes[i-1]->SetLineColor(kRed);
      boxes[i-1]->SetLineWidth(1);
      boxes[i-1]->SetX1(raaaxis[i-1]);
      boxes[i-1]->SetX2(raaaxis[i]);
      boxes[i-1]->SetY1(raaval[tempCentralityBin][i-1]-raavalsyst[tempCentralityBin][i-1]);
      boxes[i-1]->SetY2(raaval[tempCentralityBin][i-1]+raavalsyst[tempCentralityBin][i-1]);
   }
   return;
}

void getAtlas276(TGraphErrors * Atlas276, int centBin ){
        Atlas276->SetMarkerColor(kBlue);
        Atlas276->SetMarkerStyle(32);
        Atlas276->SetLineColor(kBlue);
        Atlas276->SetMarkerSize(1.4);
double p8800_d40x1y1_xval[] = { 0.5365, 0.615, 0.7050000000000001, 0.808, 0.9259999999999999, 1.0594999999999999, 1.21, 1.385, 1.5899999999999999, 
    1.825, 2.095, 2.4050000000000002, 2.755, 3.1550000000000002, 3.62, 4.15, 4.755, 5.455, 6.255, 
    7.17, 8.219999999999999, 9.39, 10.75, 12.35, 14.149999999999999, 16.2, 18.6, 21.35, 24.450000000000003, 
    28.05, 33.85, 42.6, 53.65, 67.55, 85.05, 106.9, 134.5 };
  double p8800_d40x1y1_xerrminus[] = { 0.03649999999999998, 0.04200000000000004, 0.04800000000000004, 0.05500000000000005, 0.06299999999999994, 0.0704999999999999, 0.08000000000000007, 0.09499999999999997, 0.10999999999999988, 
    0.125, 0.14500000000000024, 0.16500000000000004, 0.18500000000000005, 0.2150000000000003, 0.25, 0.28000000000000025, 0.3250000000000002, 0.375, 0.4249999999999998, 
    0.4900000000000002, 0.5599999999999987, 0.6100000000000012, 0.75, 0.8499999999999996, 0.9499999999999993, 1.0999999999999996, 1.3000000000000007, 1.4500000000000028, 1.6500000000000021, 
    1.9499999999999993, 3.8500000000000014, 4.899999999999999, 6.149999999999999, 7.75, 9.75, 12.100000000000009, 15.5 };
  double p8800_d40x1y1_xerrplus[] = { 0.03649999999999998, 0.04200000000000004, 0.04799999999999993, 0.05499999999999994, 0.06300000000000006, 0.07050000000000001, 0.08000000000000007, 0.09499999999999997, 0.1100000000000001, 
    0.125, 0.14500000000000002, 0.1649999999999996, 0.18500000000000005, 0.21499999999999986, 0.25, 0.27999999999999936, 0.3250000000000002, 0.375, 0.4249999999999998, 
    0.4900000000000002, 0.5600000000000005, 0.6099999999999994, 0.75, 0.8499999999999996, 0.9500000000000011, 1.1000000000000014, 1.2999999999999972, 1.4499999999999993, 1.6499999999999986, 
    1.9499999999999993, 3.8500000000000014, 4.899999999999999, 6.149999999999999, 7.75, 9.75, 12.099999999999994, 15.5 };
  double p8800_d40x1y1_yval[] = { 0.271, 0.289, 0.311, 0.336, 0.361, 0.388, 0.413, 0.434, 0.452, 
    0.452, 0.451, 0.427, 0.388, 0.336, 0.277, 0.221, 0.18, 0.154, 0.142, 
    0.141, 0.15, 0.164, 0.179, 0.201, 0.225, 0.245, 0.276, 0.322, 0.348, 
    0.351, 0.43, 0.51, 0.545, 0.562, 0.609, 0.614, 0.645 };
  double p8800_d40x1y1_yerrminus[] = { 0.013550227909083301, 0.013872270934854176, 0.014306324856094243, 0.01545642218553828, 0.01696752574610983, 0.01823669384637468, 0.019411929667099048, 0.019965226955774882, 0.02124573099330781, 
    0.02079447651772941, 0.02074962533106562, 0.019647615175639, 0.017467284258292703, 0.015130785486550263, 0.012479887553980606, 0.00996562982505371, 0.008309705409940835, 0.006824980471766933, 0.006475822777068562, 
    0.0063590999999999995, 0.006738323233564861, 0.007493920736170086, 0.00786989536143906, 0.00896423945463306, 0.009902556488099425, 0.011100960544025009, 0.012590562179664576, 0.016180299750004633, 0.017744587907302888, 
    0.017897558492710672, 0.022757524030526697, 0.029549208111216793, 0.039004641390480696, 0.045396939324143876, 0.057652714784301354, 0.0751115513087035, 0.11254279685968356 };
  double p8800_d40x1y1_yerrplus[] = { 0.013550227909083301, 0.013872270934854176, 0.014306324856094243, 0.01545642218553828, 0.01696752574610983, 0.01823669384637468, 0.019411929667099048, 0.019965226955774882, 0.02124573099330781, 
    0.02079447651772941, 0.02074962533106562, 0.019647615175639, 0.017467284258292703, 0.015130785486550263, 0.012479887553980606, 0.00996562982505371, 0.008309705409940835, 0.006824980471766933, 0.006475822777068562, 
    0.0063590999999999995, 0.006738323233564861, 0.007493920736170086, 0.00786989536143906, 0.00896423945463306, 0.009902556488099425, 0.011100960544025009, 0.012590562179664576, 0.016180299750004633, 0.017744587907302888, 
    0.017897558492710672, 0.022757524030526697, 0.029549208111216793, 0.039004641390480696, 0.045396939324143876, 0.057652714784301354, 0.0751115513087035, 0.11254279685968356 };
  double p8800_d40x1y1_ystatminus[] = { 7.859E-5, 8.669999999999999E-5, 9.641E-5, 1.1424000000000002E-4, 1.3356999999999998E-4, 1.5908E-4, 1.8998E-4, 2.2134E-4, 2.712E-4, 
    3.2091999999999993E-4, 3.8785999999999996E-4, 4.696999999999999E-4, 5.044E-4, 5.712000000000001E-4, 6.094000000000001E-4, 6.409E-4, 7.02E-4, 8.162E-4, 0.0010508, 
    0.0013958999999999998, 0.0019500000000000001, 0.002952, 0.002327, 0.0030150000000000003, 0.0036000000000000008, 0.004164999999999999, 0.005520000000000001, 0.008372000000000001, 0.007656, 
    0.007722, 0.0086, 0.01071, 0.010355, 0.016860000000000003, 0.029841000000000003, 0.052804, 0.094815 };
  double p8800_d40x1y1_ystatplus[] = { 7.859E-5, 8.669999999999999E-5, 9.641E-5, 1.1424000000000002E-4, 1.3356999999999998E-4, 1.5908E-4, 1.8998E-4, 2.2134E-4, 2.712E-4, 
    3.2091999999999993E-4, 3.8785999999999996E-4, 4.696999999999999E-4, 5.044E-4, 5.712000000000001E-4, 6.094000000000001E-4, 6.409E-4, 7.02E-4, 8.162E-4, 0.0010508, 
    0.0013958999999999998, 0.0019500000000000001, 0.002952, 0.002327, 0.0030150000000000003, 0.0036000000000000008, 0.004164999999999999, 0.005520000000000001, 0.008372000000000001, 0.007656, 
    0.007722, 0.0086, 0.01071, 0.010355, 0.016860000000000003, 0.029841000000000003, 0.052804, 0.094815 };
  int p8800_d40x1y1_numpoints = 37;
        for(int i = 0; i<p8800_d40x1y1_numpoints; i++){
            Atlas276->SetPoint(i,p8800_d40x1y1_xval[i],p8800_d40x1y1_yval[i]);
            Atlas276->SetPointError(i,0,p8800_d40x1y1_ystatminus[i]);
        }
   return;
}

void getAlice276(TGraphErrors * Alice276, int centBin){
        Alice276->SetMarkerStyle(27);
        Alice276->SetMarkerColor(kGreen+2);
        Alice276->SetMarkerSize(1.7);
        Alice276->SetLineColor(kGreen+2);

        double p8210_d16x1y1_xval[] = { 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 
    0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.05, 1.15, 
    1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.1, 2.3, 
    2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, 4.25, 4.75, 
    5.25, 5.75, 6.25, 6.75, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 
    13.5, 14.5, 15.5, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 
    31.0, 33.0, 35.0, 38.0, 42.5, 47.5 };
     
        double p8210_d16x1y1_yval[] = { 0.1921, 0.2001, 0.2076, 0.2163, 0.2254, 0.2366, 0.248, 0.2589, 0.27, 
    0.281, 0.2927, 0.3038, 0.3137, 0.3255, 0.3355, 0.3462, 0.3552, 0.3673, 0.385, 
    0.3968, 0.4063, 0.4153, 0.4227, 0.4284, 0.4323, 0.4392, 0.4393, 0.43, 0.4158, 
    0.4037, 0.3804, 0.3532, 0.3248, 0.2983, 0.2791, 0.2565, 0.2339, 0.2032, 0.1721, 
    0.1494, 0.139, 0.1336, 0.1323, 0.1334, 0.1431, 0.1535, 0.1614, 0.1743, 0.1847, 
    0.1965, 0.2063, 0.2246, 0.238, 0.2789, 0.2715, 0.3038, 0.3466, 0.3443, 0.3779, 
    0.3463, 0.4304, 0.3841, 0.4235, 0.3933, 0.416 };
        double p8210_d16x1y1_ystatminus[] = { 4.0E-4, 3.0E-4, 2.0E-4, 2.0E-4, 3.0E-4, 3.0E-4, 3.0E-4, 3.0E-4, 4.0E-4, 
    4.0E-4, 4.0E-4, 5.0E-4, 6.0E-4, 6.0E-4, 7.0E-4, 8.0E-4, 9.0E-4, 6.0E-4, 7.0E-4, 
    7.0E-4, 8.0E-4, 0.001, 0.001, 0.0012, 0.0012, 0.0015, 0.0015, 0.0012, 0.0012, 
    0.0014, 0.0014, 0.002, 0.0016, 0.0015, 0.0015, 0.0014, 0.0014, 0.0013, 0.0012, 
    0.001, 9.0E-4, 9.0E-4, 0.001, 0.001, 0.0014, 0.0018, 0.0023, 0.003, 0.0038, 
    0.0048, 0.0059, 0.0073, 0.0069, 0.01, 0.0127, 0.017, 0.0226, 0.0276, 0.035, 
    0.04, 0.0529, 0.0584, 0.0541, 0.0628, 0.0873 };
          double p8210_d17x1y1_yval[] = { 0.2002, 0.2064, 0.2135, 0.2225, 0.2319, 0.2432, 0.2551, 0.2669, 0.2783, 
    0.2899, 0.3023, 0.3136, 0.3235, 0.3364, 0.3471, 0.3575, 0.3679, 0.3812, 0.3985, 
    0.4128, 0.4233, 0.4319, 0.4402, 0.4467, 0.4511, 0.4594, 0.4594, 0.4507, 0.4378, 
    0.4264, 0.4044, 0.378, 0.351, 0.3249, 0.3063, 0.283, 0.2618, 0.2293, 0.1972, 
    0.1742, 0.164, 0.1563, 0.1564, 0.1566, 0.1655, 0.1793, 0.1889, 0.1994, 0.2124, 
    0.2179, 0.2394, 0.2312, 0.278, 0.2848, 0.31, 0.3346, 0.3339, 0.3364, 0.4297, 
    0.4296, 0.3144, 0.4334, 0.4193, 0.4502, 0.6208 };
          double p8210_d17x1y1_ystatminus[] = { 4.0E-4, 3.0E-4, 2.0E-4, 3.0E-4, 3.0E-4, 3.0E-4, 3.0E-4, 3.0E-4, 4.0E-4, 
    4.0E-4, 5.0E-4, 5.0E-4, 6.0E-4, 6.0E-4, 7.0E-4, 8.0E-4, 9.0E-4, 6.0E-4, 7.0E-4, 
    8.0E-4, 8.0E-4, 0.001, 0.001, 0.0013, 0.0013, 0.0016, 0.0016, 0.0013, 0.0013, 
    0.0015, 0.0015, 0.0022, 0.0017, 0.0017, 0.0016, 0.0016, 0.0016, 0.0015, 0.0014, 
    0.0011, 0.0011, 0.0011, 0.0013, 0.0012, 0.0016, 0.0022, 0.0028, 0.0036, 0.0046, 
    0.0056, 0.0071, 0.0082, 0.0084, 0.0112, 0.0152, 0.02, 0.0249, 0.0306, 0.0421, 
    0.0503, 0.0507, 0.07, 0.0605, 0.0757, 0.1204 };

        for(int i = 0; i<65; i++){
          if(centBin==0){
            Alice276->SetPoint(i,p8210_d16x1y1_xval[i],p8210_d16x1y1_yval[i]);
            Alice276->SetPointError(i,0,p8210_d16x1y1_ystatminus[i]);
          } else if(centBin==1){
            Alice276->SetPoint(i,p8210_d16x1y1_xval[i],p8210_d17x1y1_yval[i]);
            Alice276->SetPointError(i,0,p8210_d17x1y1_ystatminus[i]);
          } 
        }
        return;
}

#endif
