#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "getTrkCorr.h"
#include "TrkSettings.h"
#include "TComplex.h"
#include "Settings.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void countTracks(std::vector<std::string> inputFiles, int jobNum)
{

  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  bool doEvtPlane = true;
  float evtPlaneLow = 0;//2*(TMath::Pi()/2.0)/3;//between 0 and 90 degrees
  float evtPlaneHigh = (TMath::Pi()/2.0);//between 0 and 90 degrees
 

  bool doDebug = false;
  bool dotrkcorr = true;
  bool doNoiseFilter = false;
  float caloMatchValue = 0.5;
  float caloMatchStart = 20;
  float jetEtaSelection = 2;
  
  Settings s; 

  TH1D * hiHFDist = new TH1D("hiHFDist",";hiHF;Counts",200,0,10000);
  s.h_evtPlanePsi =  new TH1D("hiEvtPlanePsi","",100,-TMath::Pi()/2.0,TMath::Pi()/2.0); 
  s.h_Q2Mag =  new TH1D("hiQ2Mag","",200,0,200); 
  for(int i = 0; i<s.HInTriggers; i++)
  {
    for(int j = 0; j<20; j++){
      s.HIspec[i][j] = new TH2D(Form("HI_spectrum_trigger%d_cent%d",i,j),"",s.njetBins,0,s.maxJetBin,s.ntrkBins,s.xtrkbins);
      s.HIevtCount[i][j] = new TH1D(Form("HI_evtCount%d_cent%d",i,j),";max jet p_{T};N",s.njetBins,0,s.maxJetBin);
      s.HIevtCount[i][j]->SetMarkerColor(i);
      s.HIevtCount_JetVars[i][j] = new TH2D(Form("HI_evtCount_JetVars%d_cent%d",i,j),"max jet #eta;max jet p_{T};N",10,-2,2,16,40,120);
    }
  }
  for(int j = 0; j<20; j++) s.HInVtxMB[j] = new TH1D(Form("HI_nVtxMB_%d",j),"nVtx;N Events",12,0,12);
//******************************************************************************************************************************
//******************************************************************************************************************************
//******************************************************************************************************************************
  int nTrk;
  int nVtx;
  int nTrkTimesnVtx;
  bool highPurity[50000];
  float trkPt[50000];
  float trkPtError[50000];
  float trkEta[50000];
  float trkPhi[50000];
  float trkMVA[50000];
  float trkDxy1[50000];
  float trkDxyError1[50000];
  float trkDz1[50000];
  float trkDzError1[50000];
  float trkDzOverDzError[500000];
  float trkDxyOverDxyError[500000];
  float pfEcal[50000];
  float pfHcal[50000];
  float trkChi2[50000];
  float zVtx[20];
  unsigned char trkNHit[50000];
  unsigned char trkNlayer[50000];
  unsigned char trkNdof[50000];
  unsigned char trkAlgo[50000];
  unsigned char trkOriginalAlgo[50000];
  int trkCharge[50000];

  unsigned int run=0;
  unsigned int lumi=0;
  unsigned long long evt=0;

  int pVtx;
  int pBeamScrape;
  int NoiseFilter; 
  int pclusterCompatibilityFilter; 
  int pprimaryVertexFilter;  
  int phfCoincFilter3;
  int hiBin = 1;
  float hiHF = 0;

  int nref;
  float jtpt[200];
  float jteta[200];
  float jtphi[200];
  float rawpt[200];
  float chargedSum[200];
  float ecalSum[200];
  float hcalSum[200];

  int MB[20]={0};
  int j40=0;
  int j60=0;
  int j80=0;
  int t18=0;
  int t24=0;
  int t34=0;
  int t45=0;
  int t53=0;
  int HIMB[21]={0};
  int HIj40_v1=0, HIj40_v2=0;
  int HIj40=0, HIj40_c30=0, HIj40_c50=0;
  int HIj60=0, HIj60_c30=0, HIj60_c50=0;
  int HIj80=0, HIj80_c30=0, HIj80_c50=0;
  int HIj100=0,HIj100_c30=0,HIj100_c50=0;
  int HIt12=0,              HIt12_c30=0;
  int HIt18=0,              HIt18_c30=0;
  int HIt24=0,              HIt24_c30=0;
  int HIt34=0,              HIt34_c30=0;
  int HI_Muon_L2Mu20=0;

  int nPFpart;
  std::vector<float> * pfEnergy = 0;
  std::vector<float> * pfEta = 0;
  std::vector<float> * pfPhi = 0;
  std::vector<int>   * pfId = 0;
 
  TrkCorr* trkCorr;
  trkCorr = new TrkCorr("TrkCorr_May6_Iterative_pp/");

  TFile * inputFile;
  TTree * trkCh;
  TTree * jetCh;
  TTree * evtCh;
  TTree * hltCh;
  TTree * hiCh;
  TTree * pfCh;

  //for documenting which PD a file comes out of to avoid overlaps between PDs
  //0 is MB, 1 is jet40/60, 2 is jet80
  int PDindx[5000];
  int MBPDindx[5000] = {0};
  for(unsigned int i = 0; i<inputFiles.size(); i++)
  {
    if((inputFiles.at(i).find("MinimumBias") != std::string::npos) || (inputFiles.at(i).find("MinBias") != std::string::npos)) PDindx[i]=0;
    else if(inputFiles.at(i).find("HIHardProbesPeriperhal") != std::string::npos) PDindx[i]=1;
    else PDindx[i]=-1;
  }

  if(doDebug) std::cout << "\nstarting file loop \n" << std::endl;
  for(int nFile = 0; nFile<inputFiles.size(); nFile++){
    inputFile = TFile::Open(inputFiles.at(nFile).c_str(),"read");
    trkCh = (TTree*)inputFile->Get("ppTrack/trackTree");
    
    trkCh->SetBranchAddress("nTrk",&nTrk);
    trkCh->SetBranchAddress("nVtx",&nVtx);
    trkCh->SetBranchAddress("trkPt",&trkPt);
    trkCh->SetBranchAddress("trkEta",&trkEta);
    trkCh->SetBranchAddress("trkPhi",&trkPhi);
    trkCh->SetBranchAddress("highPurity",&highPurity);
    trkCh->SetBranchAddress("trkMVA",&trkMVA);
    trkCh->SetBranchAddress("trkNHit",&trkNHit);
    trkCh->SetBranchAddress("trkPtError",&trkPtError);
    trkCh->SetBranchAddress("pfHcal",&pfHcal);
    trkCh->SetBranchAddress("pfEcal",&pfEcal);
    trkCh->SetBranchAddress("trkDxy1",&trkDxy1);
    trkCh->SetBranchAddress("trkDxyError1",&trkDxyError1);
    trkCh->SetBranchAddress("trkDz1",&trkDz1);
    trkCh->SetBranchAddress("trkDzError1",&trkDzError1);
    trkCh->SetBranchAddress("trkChi2",&trkChi2);
    trkCh->SetBranchAddress("trkNlayer",&trkNlayer);
    trkCh->SetBranchAddress("trkNdof",&trkNdof);
    trkCh->SetBranchAddress("trkAlgo",&trkAlgo);
    trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo);
    trkCh->SetBranchAddress("zVtx",&zVtx);
    trkCh->SetBranchAddress("trkCharge",&trkCharge);
    trkCh->SetBranchAddress("nTrkTimesnVtx",&nTrkTimesnVtx);
    trkCh->SetBranchAddress("trkDzOverDzError",&trkDzOverDzError);
    trkCh->SetBranchAddress("trkDxyOverDxyError",&trkDxyOverDxyError); 
  
    jetCh = (TTree*)inputFile->Get("ak4CaloJetAnalyzer/t");
    jetCh->SetBranchAddress("nref",&nref);
    jetCh->SetBranchAddress("jtpt",&jtpt);
    jetCh->SetBranchAddress("jteta",&jteta);  
    jetCh->SetBranchAddress("jtphi",&jtphi);  
    jetCh->SetBranchAddress("rawpt",&rawpt);
    jetCh->SetBranchAddress("chargedSum",&chargedSum);  
    jetCh->SetBranchAddress("ecalSum",&ecalSum);
    jetCh->SetBranchAddress("hcalSum",&hcalSum);  
    trkCh->AddFriend(jetCh);
 
    //FIXME FIXME 
    evtCh = (TTree*)inputFile->Get("skimanalysis/HltTree");
    evtCh->SetBranchAddress("pPAprimaryVertexFilter",&pVtx);
    evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrape);
    evtCh->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&NoiseFilter);
  
    //PbPb-style evt sel
    /*evtCh->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);  
      evtCh->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);  
      evtCh->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3); 
    */   

    hiCh = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
    hiCh->SetBranchAddress("hiBin",&hiBin);
    hiCh->SetBranchAddress("hiHF",&hiHF);
    hiCh->SetBranchAddress("run",&run);
    hiCh->SetBranchAddress("lumi",&lumi);
    hiCh->SetBranchAddress("evt",&evt);
    evtCh->AddFriend(hiCh);
  
    if(doEvtPlane){
      pfCh = (TTree*)inputFile->Get("pfcandAnalyzer/pfTree");
      pfCh->SetBranchAddress("nPFpart",&nPFpart);
      pfCh->SetBranchAddress("pfEnergy",&pfEnergy);
      pfCh->SetBranchAddress("pfEta",&pfEta);
      pfCh->SetBranchAddress("pfPhi",&pfPhi);
      pfCh->SetBranchAddress("pfId",&pfId);
      trkCh->AddFriend(pfCh);
    }
 
    hltCh = (TTree*)inputFile->Get("hltanalysis/HltTree");
    if(PDindx[nFile]!=1){ for(int i = 0; i<3; i++) hltCh->SetBranchAddress(Form("HLT_HIL1Centralityext30100MinimumumBiasHF2AND_part%d_v1",i+1),&(HIMB[i]));}
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1",&HIj40_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1",&HIj60_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1",&HIj80_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1",&HIj100_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1",&HIj40_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1",&HIj60_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1",&HIj80_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1",&HIj100_c30);
    trkCh->AddFriend(hltCh);
  //***********************************************************************************
  //***********************************************************************
    std::cout << "starting event loop" << std::endl;
    std::cout << trkCh->GetEntries() << std::endl;
    for(int i = 0; i<trkCh->GetEntries(); i++)
    {
      for(int trig = 0; trig<21; trig++) HIMB[trig]=0;

      //if(i%1000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<" "<<std::endl;
      evtCh->GetEntry(i);
      if(!NoiseFilter && doNoiseFilter) continue;
      if(!pVtx || !pBeamScrape) continue;
      //if(!isPP && (!pclusterCompatibilityFilter || !pprimaryVertexFilter || !phfCoincFilter3)) continue; //PbPb style evt sel

      if(doDebug) std::cout << "Getting track tree" << std::endl;
      trkCh->GetEntry(i);
      if(hiBin<60) continue;
 
      bool hasGoodVtx = 0; 
      for(int vtx = 0; vtx<nVtx; vtx++){
        if(TMath::Abs(zVtx[vtx])<15) hasGoodVtx = true;
      }
      if(hasGoodVtx==false) continue;

      bool MinBias = 0;
      for(int j = 0; j<3; j++) MinBias = MinBias || (HIMB[j]==1);
      if(!MinBias && !HIj40_c50 && !HIj60_c50 && !HIj80_c50 && !HIj100_c50 && !HIj40_c30 && !HIj60_c30 && !HIj80_c30 && !HIj100_c30) continue;
  
      TComplex  Q2 = TComplex(0,0);
      float evtPlanePsi = 0;
      float Q2Mag = 0;
      int nTowersAbove3GeV[2] = {0};
      if(doEvtPlane){
        for(int q = 0; q<nPFpart; q++){
          if(!((pfId->at(q))==6 || (pfId->at(q))==7)) continue;//select HF cands
          if(TMath::Abs(pfEta->at(q)) < 3) continue; //just to be safe that we are in HF fiducial region
          if((pfEnergy->at(q)) > 3) nTowersAbove3GeV[(int)((pfEta->at(q))<0)]++;// 0 index is positive eta, 1 index is negative eta; for hf3coinc
          TComplex tempC = TComplex(pfEnergy->at(q)/TMath::CosH(pfEta->at(q)),2*pfPhi->at(q),true); 
          Q2 += tempC; 
        }
      }
      if(nTowersAbove3GeV[0]<3 || nTowersAbove3GeV[1]<3) continue;//hf3coincidence calculated from pfcands
      evtPlanePsi = 0.5*Q2.Theta();
      Q2Mag =  Q2.Rho();
      if(MinBias){
        s.h_evtPlanePsi->Fill(evtPlanePsi);
        s.h_Q2Mag->Fill(Q2Mag);
      }

      //**************************************************
      //for trigger combination with jet triggers
      float maxJtPt = 0;
      float maxJtEta = -99;
      float maxJtPhi = -99;
      for(int j=0; j<nref; j++)
      {
        if((chargedSum[j]/rawpt[j]<0.01 || TMath::Abs(jteta[j])>jetEtaSelection)) continue;
        if(jtpt[j]>maxJtPt){
          maxJtPt = jtpt[j];
          maxJtEta = jteta[j];
          maxJtPhi = jtphi[j];
        }
      }//end maxJt
 
      int PD = PDindx[nFile];
      if(MinBias==1 && PD==0)
      {
        s.HIevtCount[0][hiBin/10]->Fill(maxJtPt);
        //s.HIevtCount_JetVars[0][hiBin/10]->Fill(maxJtEta,maxJtPt);
        s.HInVtxMB[hiBin/10]->Fill(nVtx);
        hiHFDist->Fill(hiHF);
      }
      if(HIj40_c30  &&  PD==1 && hiBin>=60)  s.HIevtCount[1][hiBin/10]->Fill(maxJtPt);  
      if(HIj60_c30  &&  PD==1 && hiBin>=60)  s.HIevtCount[2][hiBin/10]->Fill(maxJtPt);  
      if(HIj80_c30  &&  PD==1 && hiBin>=60)  s.HIevtCount[3][hiBin/10]->Fill(maxJtPt);  
      if(HIj100_c30 &&  PD==1 && hiBin>=60)  s.HIevtCount[4][hiBin/10]->Fill(maxJtPt);  
      if(HIj40_c50  && !HIj40_c30 && PD==1 && hiBin>=100)  s.HIevtCount[1][hiBin/10]->Fill(maxJtPt);  
      if(HIj60_c50  && !HIj60_c30 && PD==1 && hiBin>=100)  s.HIevtCount[2][hiBin/10]->Fill(maxJtPt);  
      if(HIj80_c50  && !HIj80_c30 && PD==1 && hiBin>=100)  s.HIevtCount[3][hiBin/10]->Fill(maxJtPt);  
      if(HIj100_c50 && !HIj100_c30 && PD==1 && hiBin>=100)  s.HIevtCount[4][hiBin/10]->Fill(maxJtPt);  
      
     if(doDebug) std::cout << "Starting Track Loop" << std::endl;
     for(int j = 0; j<nTrk; j++)
     { 
       if(trkPt[j]<0.5 || trkPt[j]>=400) continue;
       if(TMath::Abs(trkEta[j])>1) continue;
       if(highPurity[j]!=1) continue; 
       if( trkPtError[j]/trkPt[j]>0.3) continue;       
      
       bool isCompatibleWithVertex = false;
       for(int v = 0; v<nVtx; v++){
         if(TMath::Abs(zVtx[v])>15) continue;
         if(TMath::Abs(trkDxyOverDxyError[j*nVtx+v])<3 && TMath::Abs(trkDzOverDzError[j*nVtx+v])<3){
           isCompatibleWithVertex = true;
           break;
         }
       }          
       if(!isCompatibleWithVertex) continue;
       if(doEvtPlane){
         float dPhi1 = TMath::ACos(TMath::Cos(trkPhi[j]-evtPlanePsi)); 
         float dPhi2 = TMath::ACos(TMath::Cos(trkPhi[j]-(evtPlanePsi+TMath::Pi())));//reflected part of the evt plane (for tracks that are pi away from evt plane angle)
         float mindPhi = TMath::Min(dPhi1,dPhi2); 
         if(mindPhi<evtPlaneLow || mindPhi>evtPlaneHigh) continue;
       }
 
       //if(TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
       //if(trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>0.15) continue; 
       //if(trkPtError[j]/trkPt[j]>0.1) continue;       
       //if(trkNHit[j]<11 && trkPt[j]>0.7) continue; 
   
       float rmin=999;
       for(int jt=0; jt<nref; jt++)
       {
         if((chargedSum[jt]/rawpt[jt]<0.01 || TMath::Abs(jteta[jt])>2)) continue;
         if(jtpt[jt]<50) continue;
         float R = TMath::Power(jteta[jt]-trkEta[j],2) + TMath::Power(TMath::ACos(TMath::Cos(jtphi[jt]-trkPhi[j])),2);
         if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
       }
  
       if(doDebug) std::cout << "before Corr: "<< hiBin <<" "<< trkEta[j] <<" "<< trkPhi[j]<<" " << trkPt[j] <<" " << rmin << std::endl;
       float correction=1;
       correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],1,rmin);
       if(doDebug) std::cout << "after Corr" << std::endl;
 
       float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
       if(!(trkPt[j]<caloMatchStart || (Et>caloMatchValue*trkPt[j]))){
         continue; //Calo Matching
       }
       
       //dividing by pt at bin center instead of track by track pt (just a convention)
       float binCenter;
       binCenter = s.HIspec[0][0]->GetYaxis()->GetBinCenter(s.HIspec[0][0]->GetYaxis()->FindBin(trkPt[j]));
       
       if(!dotrkcorr) correction=1;
       if(MinBias==1 && PD==0) s.HIspec[0][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
       if(HIj40_c30  && PD==1 && hiBin>=60)   s.HIspec[1][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
       if(HIj60_c30  && PD==1 && hiBin>=60)   s.HIspec[2][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
       if(HIj80_c30  && PD==1 && hiBin>=60)   s.HIspec[3][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
       if(HIj100_c30 && PD==1 && hiBin>=60) s.HIspec[4][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);  
       if(HIj40_c50  && !HIj40_c30 && PD==1 && hiBin>=100)   s.HIspec[1][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
       if(HIj60_c50  && !HIj60_c30 && PD==1 && hiBin>=100)   s.HIspec[2][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
       if(HIj80_c50  && !HIj80_c30 && PD==1 && hiBin>=100)   s.HIspec[3][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
       if(HIj100_c50 && !HIj100_c30 && PD==1 && hiBin>=100) s.HIspec[4][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);  
      } //end trk loop
    }//end event loop
    inputFile->Close();
  }//end file loop

  if(doDebug) std::cout << "Starting Output" << std::endl;
  TFile * outF;
  outF = TFile::Open(Form("PbPb_output_%d.root",jobNum),"recreate");
  outF->cd();
  hiHFDist->Write();
  s.h_evtPlanePsi->Write();
  s.h_Q2Mag->Write();
  for(int i = 0; i<s.HInTriggers; i++)
  {
    for(int j = 0; j<20; j++){
      s.HIspec[i][j]->Write();
      s.HIevtCount[i][j]->Write();
      //s.HIevtCount_JetVars[i][j]->Write();
    }
  }
  for(int i = 0; i<s.HInTriggers_trk; i++)
  {
    for(int j = 0; j<20; j++){
      //s.HIspec_trk[i][j]->Write();
      //s.HIevtCount_trk[i][j]->Write();
    }
  }
  for(int j = 0; j<20; j++){
    s.HInVtxMB[j]->Write();
    //s.HInVtxMB_trk[j]->Write();
  }
  outF->Close(); 
}



//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
int main(int argc, const char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: countTracks <fileList>  <job>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  int job = std::atoi(argv[2]);
  int totalJobs = std::atoi(argv[3]);
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      if(line%totalJobs==job) listOfFiles.push_back(buffer);
      line++;
    }
  }
   
  countTracks(listOfFiles,job);
  return 0; 
}
