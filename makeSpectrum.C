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
//#include "ppPlotting.C"
#include "PbPbPlotting.C"
#include "prettyPlotting.C"

void makeSpectrum()
{
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  Settings s;

//output hist for jet triggered
  s.pp = new TH1D("ppTrackSpectrum",";p_{T} (GeV);E#frac{d^{3}#sigma}{d^{3}p} (mb/GeV^{2})",s.ntrkBins,s.xtrkbins);
  for(int i = 0; i<s.nTriggers; i++) s.ppByTrigger[i] = new TH1D(Form("ppTrackSpectrumByTrigger%d",i),"",s.ntrkBins,s.xtrkbins);
  for(int i = 0; i<s.nTriggers; i++) s.ppUsedByTrigger[i] = new TH1D(Form("ppUsedTrackSpectrumByTrigger%d",i),"",s.ntrkBins,s.xtrkbins);
  s.ppJets = new TH1D("ppJetSpectrum",";Leading Jet P_{T} (GeV);#sigma (mb)",s.njetBins,0,s.maxJetBin);
  for(int i = 0; i<s.nTriggers; i++) s.ppJetsByTrigger[i] = new TH1D(Form("ppJetSpectrumByTrigger%d",i),"",s.njetBins,0,s.maxJetBin);

//same for PbPb
  for(int c = 0 ; c<s.nCentBins; c++){
    s.HI[c] = new TH1D(Form("PbPbTrackSpectrum_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]),";p_{T} (GeV);E#frac{d^{3}#sigma}{d^{3}p} (mb/GeV^{2})",s.ntrkBins,s.xtrkbins);
    for(int i = 0; i<s.HInTriggers; i++) s.HIByTrigger[i][c] = new TH1D(Form("PbPbTrackSpectrumByTrigger%d_%d_%d",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.ntrkBins,s.xtrkbins);
    for(int i = 0; i<s.HInTriggers; i++) s.HIUsedByTrigger[i][c] = new TH1D(Form("PbPbUsedTrackSpectrumByTrigger%d_%d_%d",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.ntrkBins,s.xtrkbins);
    s.HIJets[c] = new TH1D(Form("PbPbJetSpectrum_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]),";Leading Jet P_{T} (GeV);#sigma (mb)",s.njetBins,0,s.maxJetBin);
    for(int i = 0; i<s.HInTriggers; i++) s.HIJetsByTrigger[i][c] = new TH1D(Form("PbPbJetSpectrumByTrigger%d_%d_%d",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.njetBins,0,s.maxJetBin);
  }


  //loading files
  TFile * inFile = TFile::Open("ppInput.root","read");
  for(int i = 0; i<s.nTriggers; i++)
  {
    s.spec[i] = (TH2D*) inFile->Get(Form("spectrum_trigger%d",i));
    s.evtCount[i] = (TH1D*) inFile->Get(Form("evtCount%d",i));
    s.spec[i]->SetDirectory(0);
    s.evtCount[i]->SetDirectory(0);
  }
  s.nVtxMB = (TH1D*) inFile->Get("nVtxMB");
  s.nVtxMB->SetDirectory(0);
  inFile->Close();
  
  inFile = TFile::Open("CurrentInput.root","read");
  for(int c = 0; c<20; c++)
  {
    for(int i = 0; i<s.HInTriggers; i++)
    {
      if(s.doBetterHITrig && c<6 && i==1){//replaces a few triggers for better stats
          s.HIspec[i][c] = (TH2D*) inFile->Get(Form("HI_spectrum_trigger%d_cent%d",i-1,c));
          s.HIevtCount[i][c] = (TH1D*) inFile->Get(Form("HI_evtCount%d_cent%d",i-1,c));
      /*}else if(s.doBetterHITrig && c>=6 && c<10 && i==s.HInTriggers-1){
          s.HIspec[i][c] = (TH2D*) inFile->Get(Form("HI_spectrum_trigger%d_cent%d",s.HInTriggers-2,c));
          s.HIevtCount[i][c] = (TH1D*) inFile->Get(Form("HI_evtCount%d_cent%d",s.HInTriggers-2,c));
      }else if(s.doBetterHITrig && c>=10 && i>=s.HInTriggers-2){
        s.HIspec[i][c] = (TH2D*) inFile->Get(Form("HI_spectrum_trigger%d_cent%d",s.HInTriggers-3,c));
        s.HIevtCount[i][c] = (TH1D*) inFile->Get(Form("HI_evtCount%d_cent%d",s.HInTriggers-3,c));*/ //removes other unprescaled triggers, should be equivilent to below (maybe changes if have failed jobs)
      }else{//behavior without better trigger handling
        s.HIspec[i][c] = (TH2D*) inFile->Get(Form("HI_spectrum_trigger%d_cent%d",i,c));
        s.HIevtCount[i][c] = (TH1D*) inFile->Get(Form("HI_evtCount%d_cent%d",i,c));
      }
      s.HIspec[i][c]->SetDirectory(0);
      s.HIevtCount[i][c]->SetDirectory(0);
    }
    s.HInVtxMB[c] = (TH1D*) inFile->Get(Form("HI_nVtxMB_%d",c));
    s.HInVtxMB[c]->SetDirectory(0);
  }
  inFile->Close();

  //calculation of overlaps
  float scale[s.nTriggers];
  float HIscale[s.HInTriggers][3];
  s.h_scale = new TH2D("h_scale","h_scale",10,0,10,10,0,10);
  s.h_HIscale = new TH2D("h_HIscale","h_HIscale",10,0,10,10,0,10);
  s.h_normErr = new TH2D("h_normErr","h_normErr",10,0,10,10,0,10);
  s.h_HInormErr = new TH2D("h_HInormErr","h_HInormErr",10,0,10,10,0,10);
  //*****************************************************************************************************************************************************
  //pp jet triggers
  for(int i = 0; i<s.nTriggers; i++)
  {
    scale[i] = 1;//using 68mb as inelastic pp xsection
    for(int j = 0; j<i; j++){
      scale[i] = scale[i]*s.evtCount[j]->Integral(s.evtCount[j]->FindBin(s.triggerOverlapBins[j+1]),s.evtCount[j]->FindBin(s.maxJetBin))/(double)s.evtCount[j+1]->Integral(s.evtCount[j+1]->FindBin(s.triggerOverlapBins[j+1]),s.evtCount[j+1]->FindBin(s.maxJetBin));
      if(i>0) s.h_normErr->SetBinContent(i,1,TMath::Power(1.0/s.evtCount[j]->Integral(s.evtCount[j]->FindBin(s.triggerOverlapBins[j+1]),s.evtCount[j]->FindBin(s.maxJetBin))+1.0/s.evtCount[j+1]->Integral(s.evtCount[j+1]->FindBin(s.triggerOverlapBins[j+1]),s.evtCount[j+1]->FindBin(s.maxJetBin)),0.5));//error on last overlap
    }
    std::cout << scale[i] << std::endl;
    s.spec[i]->Scale(scale[i]);

    //total spectrum
    for(int j = s.evtCount[i]->FindBin(s.triggerBins[i]); j<s.evtCount[i]->FindBin(s.triggerBins[i+1]); j++)
    {
      for(int k = 1; k<s.pp->GetSize()+1; k++)
      {
        s.pp->SetBinContent(k,s.pp->GetBinContent(k)+s.spec[i]->GetBinContent(j,k)); 
        s.pp->SetBinError(k,TMath::Power(TMath::Power(s.pp->GetBinError(k),2)+TMath::Power(s.spec[i]->GetBinError(j,k),2),0.5)); 
        s.ppUsedByTrigger[i]->SetBinContent(k,s.ppUsedByTrigger[i]->GetBinContent(k)+s.spec[i]->GetBinContent(j,k)); 
        s.ppUsedByTrigger[i]->SetBinError(k,TMath::Power(TMath::Power(s.ppUsedByTrigger[i]->GetBinError(k),2)+TMath::Power(s.spec[i]->GetBinError(j,k),2),0.5)); 
      }
      s.ppJets->SetBinContent(j,s.evtCount[i]->GetBinContent(j)*scale[i]);
      s.ppJets->SetBinError(j,s.evtCount[i]->GetBinError(j)*scale[i]);
    }
   
    //spectrum for each trigger
    for(int j = s.evtCount[i]->FindBin(0); j<s.evtCount[i]->FindBin(s.maxJetBin); j++)
    {
      for(int k = 1; k<s.ppByTrigger[i]->GetSize()+1; k++)
      {
        s.ppByTrigger[i]->SetBinContent(k,s.ppByTrigger[i]->GetBinContent(k)+s.spec[i]->GetBinContent(j,k)); 
        s.ppByTrigger[i]->SetBinError(k,TMath::Power(TMath::Power(s.ppByTrigger[i]->GetBinError(k),2)+TMath::Power(s.spec[i]->GetBinError(j,k),2),0.5)); 
      }
      s.ppJetsByTrigger[i]->SetBinContent(j,s.evtCount[i]->GetBinContent(j)*scale[i]);
      s.ppJetsByTrigger[i]->SetBinError(j,s.evtCount[i]->GetBinError(j)*scale[i]);
    }
  }
  

  TH1D * tempEvtCount[s.HInTriggers][3];
  int combinationCentUpperBoundary[3] = {5,9,19};
  int combinationCentLowerBoundary[3] = {0,6,10};
  for(int i = 0; i<s.HInTriggers; i++)
  {
    for(int m = 0; m<3; m++){
      HIscale[i][m] = 1;//using 68mb as inelastic pp xsection
      tempEvtCount[i][m] = (TH1D*)s.HIevtCount[i][0]->Clone(Form("HItempEvtCount%d",m));
      tempEvtCount[i][m]->Reset();
      for(int c = combinationCentUpperBoundary[m]; c>=combinationCentLowerBoundary[m]; c--) tempEvtCount[i][m]->Add(s.HIevtCount[i][c]);
      for(int j = 0; j<i; j++){
        HIscale[i][m] = HIscale[i][m]*tempEvtCount[j][m]->Integral(tempEvtCount[j][m]->FindBin(s.HItriggerOverlapBins[j+1]),tempEvtCount[j][m]->FindBin(s.maxJetBin))/(double)tempEvtCount[j+1][m]->Integral(tempEvtCount[j+1][m]->FindBin(s.HItriggerOverlapBins[j+1]),tempEvtCount[j+1][m]->FindBin(s.maxJetBin));
        if(i>0) s.h_HInormErr->SetBinContent(i,m+1,TMath::Power(1.0/tempEvtCount[j][m]->Integral(tempEvtCount[j][m]->FindBin(s.HItriggerOverlapBins[j+1]),tempEvtCount[j][m]->FindBin(s.maxJetBin))+1.0/tempEvtCount[j+1][m]->Integral(tempEvtCount[j+1][m]->FindBin(s.HItriggerOverlapBins[j+1]),tempEvtCount[j+1][m]->FindBin(s.maxJetBin)),0.5));
        if(s.doBetterHITrig && i==1 && m==0) s.h_HInormErr->SetBinContent(i,m+1,0);//remove jet40
        if(s.doBetterHITrig && i==s.HInTriggers-1 && m>0)  s.h_HInormErr->SetBinContent(i-1,m+1,s.h_HInormErr->GetBinContent(i,m+1));//remove jet100
        if(s.doBetterHITrig && i==s.HInTriggers-2 && m>1)  s.h_HInormErr->SetBinContent(i-1,m+1,s.h_HInormErr->GetBinContent(i,m+1));//remove jet80
      }
      std::cout <<"PbPb scale: Trigger and cent region "<< i<<" "<<m<<" "<<HIscale[i][m] << std::endl;
    }

    for(int c = 0; c<20; c++){
      if(c<6){
        s.HIspec[i][c]->Scale(HIscale[i][0]);
        s.HIevtCount[i][c]->Scale(HIscale[i][0]);
      }
      else if(c<10){
        s.HIspec[i][c]->Scale(HIscale[i][1]);
        s.HIevtCount[i][c]->Scale(HIscale[i][1]);
      }else{
        s.HIspec[i][c]->Scale(HIscale[i][2]);
        s.HIevtCount[i][c]->Scale(HIscale[i][2]);
      } 
    }
    for(int c = 20; c<s.nCentBins; c++){
      s.HIspec[i][c] = (TH2D*)s.HIspec[i][0]->Clone(Form("HI_spectrum_trigger%d_cent%d",i,c));
      s.HIspec[i][c]->Reset();
      s.HIevtCount[i][c] = (TH1D*)s.HIevtCount[i][0]->Clone(Form("HI_evtCount%d_cent%d",i,c));
      s.HIevtCount[i][c]->Reset();
      for(int cc = s.lowCentBin[c]; cc<s.highCentBin[c]; cc++){
        s.HIspec[i][c]->Add(s.HIspec[i][cc]);
        s.HIevtCount[i][c]->Add(s.HIevtCount[i][cc]);
      }
    }
     
    //total spectrum
    for(int c = 0; c<s.nCentBins; c++){
      for(int j = s.HIevtCount[i][c]->FindBin(s.HItriggerBins[i]); j<s.HIevtCount[i][c]->FindBin(s.HItriggerBins[i+1]); j++)
      {
        for(int k = 1; k<s.HI[c]->GetSize()+1; k++)
        {
          s.HI[c]->SetBinContent(k,s.HI[c]->GetBinContent(k)+s.HIspec[i][c]->GetBinContent(j,k)); 
          s.HI[c]->SetBinError(k,TMath::Power(TMath::Power(s.HI[c]->GetBinError(k),2)+TMath::Power(s.HIspec[i][c]->GetBinError(j,k),2),0.5)); 
          s.HIUsedByTrigger[i][c]->SetBinContent(k,s.HIUsedByTrigger[i][c]->GetBinContent(k)+s.HIspec[i][c]->GetBinContent(j,k)); 
          s.HIUsedByTrigger[i][c]->SetBinError(k,TMath::Power(TMath::Power(s.HIUsedByTrigger[i][c]->GetBinError(k),2)+TMath::Power(s.HIspec[i][c]->GetBinError(j,k),2),0.5)); 
        }
        s.HIJets[c]->SetBinContent(j,s.HIevtCount[i][c]->GetBinContent(j));
        s.HIJets[c]->SetBinError(j,s.HIevtCount[i][c]->GetBinError(j));
      }
     
      //spectrum for each trigger
      for(int j = s.HIevtCount[i][c]->FindBin(0); j<s.HIevtCount[i][c]->FindBin(s.maxJetBin); j++)
      {
        for(int k = 1; k<s.HIByTrigger[i][c]->GetSize()+1; k++)
        {
          s.HIByTrigger[i][c]->SetBinContent(k,s.HIByTrigger[i][c]->GetBinContent(k)+s.HIspec[i][c]->GetBinContent(j,k)); 
          s.HIByTrigger[i][c]->SetBinError(k,TMath::Power(TMath::Power(s.HIByTrigger[i][c]->GetBinError(k),2)+TMath::Power(s.HIspec[i][c]->GetBinError(j,k),2),0.5)); 
        }
        s.HIJetsByTrigger[i][c]->SetBinContent(j,s.HIevtCount[i][c]->GetBinContent(j));
        s.HIJetsByTrigger[i][c]->SetBinError(j,s.HIevtCount[i][c]->GetBinError(j));
      }
    }//cent loop closed
  }

  //*********************************************************************************************************************************************88
  //Divide by bin width and jacobian stuff
  for(int i = 1; i<s.pp->GetSize()+1; i++)
  {
    s.pp->SetBinContent(i,s.pp->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
    s.pp->SetBinError(i,s.pp->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
    for(int j = 0; j<s.nTriggers; j++)
    {
      s.ppByTrigger[j]->SetBinContent(i,s.ppByTrigger[j]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.ppByTrigger[j]->SetBinError(i,s.ppByTrigger[j]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      s.ppUsedByTrigger[j]->SetBinContent(i,s.ppUsedByTrigger[j]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.ppUsedByTrigger[j]->SetBinError(i,s.ppUsedByTrigger[j]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
    } 
  }

  for(int c = 0; c<s.nCentBins; c++)
  {
    for(int i = 1; i<s.HI[c]->GetSize()+1; i++)
    {
      s.HI[c]->SetBinContent(i,s.HI[c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.HI[c]->SetBinError(i,s.HI[c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      for(int j = 0; j<s.HInTriggers; j++)
      {
        s.HIByTrigger[j][c]->SetBinContent(i,s.HIByTrigger[j][c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
        s.HIByTrigger[j][c]->SetBinError(i,s.HIByTrigger[j][c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
        s.HIUsedByTrigger[j][c]->SetBinContent(i,s.HIUsedByTrigger[j][c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
        s.HIUsedByTrigger[j][c]->SetBinError(i,s.HIUsedByTrigger[j][c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      } 
    }
  }
 
  //calculate total number of verticies from MB events
  int nVtx = 0, nVtx_trk = 0;
  for(int i = 1; i<s.nVtxMB->GetSize()+1;i++) nVtx = nVtx+i*s.nVtxMB->GetBinContent(s.nVtxMB->FindBin(i));
  
  TFile * outF = TFile::Open("Spectra.root","recreate");
  //propagating nomalization errors
  //**************************************************************SYSTEMATICS***************************************************
  for(int m = 1; m<4; m++){
    for(int i = 1; i<10; i++){
      s.h_normErr->SetBinContent(i,m,TMath::Power(TMath::Power(s.h_normErr->GetBinContent(i,m),2)+TMath::Power(s.h_normErr->GetBinContent(i-1,m),2),0.5));
      s.h_HInormErr->SetBinContent(i,m,TMath::Power(TMath::Power(s.h_HInormErr->GetBinContent(i,m),2)+TMath::Power(s.h_HInormErr->GetBinContent(i-1,m),2),0.5));
    }
    //invert pp because lumi normalized
    for(int i = 0; i<10; i++){
      s.h_normErr->SetBinContent(i,m,TMath::Power(-TMath::Power(s.h_normErr->GetBinContent(i,m),2)+TMath::Power(s.h_normErr->GetBinContent(s.nTriggers,m),2),0.5));
    }
  }
  s.h_normErr->Write();
  s.h_HInormErr->Write();
 
  //adding up normalization systematics 
  s.h_normSyst = (TH1D*)s.pp->Clone("h_normSyst");
  s.h_normSyst->Reset();
  for(int c = 0; c<s.nCentBins; c++){
    s.h_HInormSyst[c] = (TH1D*) s.HI[c]->Clone(Form("h_HInormSyst_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]));
    s.h_HInormSyst[c]->Reset();
  }
  for(int j = 1; j<s.h_normSyst->GetSize()-1; j++){
    for(int i = 0; i<s.nTriggers; i++) s.h_normSyst->SetBinContent(j,s.h_normSyst->GetBinContent(j)+s.h_normErr->GetBinContent(i,1)*s.ppUsedByTrigger[i]->GetBinContent(j)/s.pp->GetBinContent(j)); 
  }
  s.h_normSyst->Write();
  s.h_normSyst->SetDirectory(0);
 
  for(int c = 0; c<s.nCentBins; c++){ 
    for(int cc = s.lowCentBin[c]; cc<s.highCentBin[c]; cc++){ 
      for(int j = 1; j<s.h_HInormSyst[c]->GetSize()-1; j++){
        for(int i = 0; i<s.HInTriggers; i++) if(s.HI[c]->GetBinContent(j)!=0) s.h_HInormSyst[c]->SetBinContent(j,s.h_HInormSyst[c]->GetBinContent(j)+s.h_HInormErr->GetBinContent(i,(cc<6)?1:((cc<10)?2:3))*s.HIUsedByTrigger[i][cc]->GetBinContent(j)/s.HI[c]->GetBinContent(j));
      }
    }
    s.h_HInormSyst[c]->Write();
    s.h_HInormSyst[c]->SetDirectory(0);
  }
  
  //*************************OUTPUT**********************************************************
  
  s.pp_perMBTrigger = (TH1D*)s.pp->Clone("pp_NotperMBTrigger");
  s.pp_perMBTrigger->Scale(1.0/((scale[s.nTriggers-1])*(27.391*1000000000)));//25.775 pb GOLDEN JSON LUMI for jet80 trigger

  s.pp_perMBTrigger->Write();
  s.pp_perMBTrigger->SetDirectory(0);
  s.pp->Scale(1/(float)nVtx);
  s.pp->Write();
  
  s.ppJets->SetMarkerSize(0);
  s.ppJets->Scale(1/(float)nVtx);
  s.ppJets->Write();
  
  for(int j = 0; j<s.nTriggers; j++)
  {
    s.ppByTrigger[j]->Scale(1/(float)nVtx);
    s.ppUsedByTrigger[j]->Scale(1/(float)nVtx);
    s.ppJetsByTrigger[j]->Scale(1/(float)nVtx);
    s.ppByTrigger[j]->Write();
    s.ppUsedByTrigger[j]->Write();
    s.ppJetsByTrigger[j]->Write();
  }
  for(int c = 0; c<s.nCentBins; c++){
    s.HI_perMBTrigger[c] = (TH1D*)s.HI[c]->Clone(Form("PbPb_NotperMBTrigger_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]));
    s.HI_perMBTrigger[c]->Write();
   
    double nMBInCentRange = 0;
    for(int cc = s.lowCentBin[c]; cc<s.highCentBin[c]; cc++) nMBInCentRange += s.HInVtxMB[cc]->GetEntries();
    s.HI[c]->Scale(1/(nMBInCentRange));
    s.HIJets[c]->Scale(1/(nMBInCentRange));
    s.HI[c]->Write();
    s.HIJets[c]->Write();
    for(int j = 0; j<s.HInTriggers; j++){
      s.HIByTrigger[j][c]->Scale(1/(nMBInCentRange));
      s.HIUsedByTrigger[j][c]->Scale(1/(nMBInCentRange));
      s.HIJetsByTrigger[j][c]->Scale(1/(nMBInCentRange));
      s.HIByTrigger[j][c]->Write();
      s.HIUsedByTrigger[j][c]->Write();
      s.HIJetsByTrigger[j][c]->Write();
    }
  } 
  
  for(int i = 0; i<s.nTriggers; i++)      s.h_scale->SetBinContent(i+1,1,scale[i]);
  for(int i = 0; i<3; i++){
    for(int j = 0; j<s.HInTriggers; j++)     s.h_HIscale->SetBinContent(j+1,i+1,HIscale[j][i]);
  }
  s.h_scale->Write();
  s.h_HIscale->Write();
 
  makePlotsPbPb(s);
  outF->Close();
  //makePlotsPP(s);
  prettyPlotting(s);
}
