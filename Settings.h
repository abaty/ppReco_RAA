#ifndef RAAINPUTSETTINGS
#define RAAINPUTSETTINGS

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

class Settings {
  public:
  /*static const int ntrkBins = 87;
  double xtrkbins[ntrkBins+1] = { 1,2,3,4,5,6,7,8,9,10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80,83,86,89, 92, 95, 98, 101, 104, 107, 110, 115, 120, 125, 130, 135, 140, 145, 150};
  double xt_xtrkbins[ntrkBins+1] = { 1*2/5020.0,2*2/5020.0,3*2/5020.0,4*2/5020.0,5*2/5020.0,6*2/5020.0,7*2/5020.0,8*2/5020.0,9*2/5020.0,10*2/5020.0, 11*2/5020.0, 12*2/5020.0, 13*2/5020.0, 14*2/5020.0, 15*2/5020.0, 16*2/5020.0, 17*2/5020.0, 18*2/5020.0, 19*2/5020.0, 20*2/5020.0, 21*2/5020.0, 22*2/5020.0, 23*2/5020.0, 24*2/5020.0, 25*2/5020.0, 26*2/5020.0, 27*2/5020.0, 28*2/5020.0, 29*2/5020.0, 30*2/5020.0, 31*2/5020.0, 32*2/5020.0, 33*2/5020.0, 34*2/5020.0, 35*2/5020.0, 36*2/5020.0, 37*2/5020.0, 38*2/5020.0, 39*2/5020.0, 40*2/5020.0, 41*2/5020.0, 42*2/5020.0, 43*2/5020.0, 44*2/5020.0, 45*2/5020.0, 46*2/5020.0, 47*2/5020.0, 48*2/5020.0, 49*2/5020.0, 50*2/5020.0, 51*2/5020.0, 52*2/5020.0, 53*2/5020.0, 54*2/5020.0, 55*2/5020.0, 56*2/5020.0, 57*2/5020.0, 58*2/5020.0, 59*2/5020.0, 60*2/5020.0, 62*2/5020.0, 64*2/5020.0, 66*2/5020.0, 68*2/5020.0, 70*2/5020.0, 72*2/5020.0, 74*2/5020.0, 76*2/5020.0, 78*2/5020.0, 80*2/5020.0,83*2/5020.0,86*2/5020.0,89*2/5020.0, 92*2/5020.0, 95*2/5020.0, 98*2/5020.0, 101*2/5020.0, 104*2/5020.0, 107*2/5020.0, 110*2/5020.0, 115*2/5020.0, 120*2/5020.0, 125*2/5020.0, 130*2/5020.0, 135*2/5020.0, 140*2/5020.0, 145*2/5020.0, 150*2/5020.0};
 */

  //for D comparison
  /*static const int ntrkBins = 14;
  double xtrkbins[ntrkBins+1] = {2.,3.,4.,5.,6.,8.,10.,12.5,15.0,20.,25.,30.,40.,60.,100};
  double xt_xtrkbins[ntrkBins+1] = {2.*2/5020.0,3.*2/5020.0,4.*2/5020.0,5.*2/5020.0,6.*2/5020.0,8.*2/5020.0,10.*2/5020.0,12.5*2/5020.0,15.0*2/5020.0,20.*2/5020.0,25.*2/5020.0,30.*2/5020.0,40.*2/5020.0,60.*2/5020.0,100*2/5020.0};
 */
   //for RAA paper
  static const int ntrkBins = 37;
  double xtrkbins[ntrkBins+1] = { 0.5,0.6, 0.7 , 0.8 , 0.9 , 1.0 , 1.1 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 3.2 , 4.0 , 4.8 , 5.6 , 6.4 , 7.2 , 9.6 , 12.0, 14.4,19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8,73.6,86.4,103.6,120.8,140,165,250,400};
  double xt_xtrkbins[ntrkBins+1] = { 0.5*2/5020.0,0.6*2/5020.0, 0.7*2/5020.0, 0.8*2/5020.0, 0.9*2/5020.0, 1.0*2/5020.0, 1.1*2/5020.0, 1.2*2/5020.0, 1.4*2/5020.0, 1.6*2/5020.0, 1.8*2/5020.0, 2.0*2/5020.0, 2.2*2/5020.0, 2.4*2/5020.0, 3.2*2/5020.0, 4.0*2/5020.0, 4.8*2/5020.0, 5.6*2/5020.0, 6.4*2/5020.0, 7.2*2/5020.0, 9.6*2/5020.0, 12.0*2/5020.0, 14.4*2/5020.0,19.2*2/5020.0, 24.0*2/5020.0, 28.8*2/5020.0, 35.2*2/5020.0, 41.6*2/5020.0, 48.0*2/5020.0, 60.8*2/5020.0,73.6*2/5020.0,86.4*2/5020.0,103.6*2/5020.0,120.8*2/5020.0,140*2/5020.0,165*2/5020.0,250*2/5020.0,400*2/5020.0};
 
  static const int nEvtPlaneBins = 9; 

  static const int nTriggers = 4;
  
  double triggerBins[nTriggers+1] = {0,60,80,100,1200};
  double triggerOverlapBins[nTriggers] = {0,60,80,100};
  /*double triggerBins[nTriggers+1] = {0,55,75,95,1200};
  double triggerOverlapBins[nTriggers] = {0,55,75,95};*/

  bool   doBetterHITrig = true;//removes jet40 from 0-30% and removes jet80/100 for 50-100 and 30-100%
  static const int HInTriggers = 5;
  static const int offsetFor30_50_pprecoPbPb = 60;
  double HItriggerBins[HInTriggers+1] = {0,60,80,100,120,1200};
  double HItriggerOverlapBins[HInTriggers] = {0,60,80,100,120};
  /*double HItriggerBins[HInTriggers+1] = {0,55,75,95,115,1200};
  double HItriggerOverlapBins[HInTriggers] = {0,55,75,95,115};*/
  
  static const int njetBins = 240;
  static const int maxJetBin = 1200;
  
  static const int nTriggers_trk = 6;
  double triggerBins_trk[nTriggers_trk+1] = {0,20,26,36,47,55,500};
  double triggerOverlapBins_trk[nTriggers_trk] = {0,20,26,36,47,55};
  static const int HInTriggers_trk = 5;
  double HItriggerBins_trk[HInTriggers_trk+1] = {0,14,20,35,50,500};
  double HItriggerOverlapBins_trk[HInTriggers_trk] = {0,14,20,35,50};

  static const int nTrktriggerBins = 500;
  static const int maxTrktriggerBin = 500;
  
  static const int nCentBins = 33;
  int lowCentBin[nCentBins] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,0,6,10,2,6,10,14,16,0,0,14,0,10};
  int highCentBin[nCentBins] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,20,20,20,6,10,14,16,18,6,10,18,2,18};
  double nColl[nCentBins] = {1819,1433,1127,882,685.9,526.5,399.3,297.5,217.1,155.1,107.9,73.51,48.76,31.46,19.69,12.02,7.042,3.974,2.12,1.164,392.4,98.33,30.74,805.7,267.25,65.415,15.87,5.502,1079.1,754,10.686, 1626,38.04};
  double TAA[nCentBins] = {25.98,20.46,1127./70.0,882./70.0,685.9/70.0,526.5/70.0,399.3/70.0,297.5/70.0,217.1/70.0,155.1/70.0,107.9/70.0,73.51/70.0,48.76/70.0,31.46/70.0,19.69/70.0,12.02/70.0,7.042/70.0,3.974/70.0,2.12/70.0,1.164/70.0,5.607,98.33/70.0,30.74/70.0,11.51,3.819,0.9345,15.87/70.0,5.502/70.0,1079.1/70.0,754/70.0,0.1525, 23.22, 0.5435};//Ncoll/70.0for TAA in most places, otherwise we are using offical numbers from cent group if no 70.0is there
  //double TAAuncert[nCentBins] = {1.7,1.7,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,8.9,5,5,2.4,4.9,9.6,5,5,5,5,16.,1.7,12.8};//assume 5% uncert for 'unofficial' values for now (5.0 is official, 5 is not in this list)
  double TAAuncert[nCentBins] = {1.8,1.9,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2.8,5,5,2.7,5.5,10.3,5,5,5,5,16.1,1.9,12.8};//assume 5% uncert for 'unofficial' values for now (5.0 is official, 5 is not in this list)
  double TAAuncertlo[nCentBins] = {3.0,3.0,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3.4,5,5,3.4,5.4,9.5,5,5,5,5,14.1,3.0,12.8};//assume 5% uncert for 'unofficial' values for now (5.0 is official, 5 is not in this list)
  /*static const int nCentBins = 31;
  int lowCentBin[nCentBins] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,0,6,10,2,6,10,14,16,0,0,14};
  int highCentBin[nCentBins] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,20,20,20,6,10,14,16,18,6,10,18};
  double nColl[nCentBins] = {1819,1433,1127,882,685.9,526.5,399.3,297.5,217.1,155.1,107.9,73.51,48.76,31.46,19.69,12.02,7.042,3.974,2.12,1.164,392.4,98.33,30.74,805.7,267.25,65.415,15.87,5.502,1079.1,754,10.686};
  double TAA[nCentBins] = {25.98,20.46,1127./70.0,882./70.0,685.9/70.0,526.5/70.0,399.3/70.0,297.5/70.0,217.1/70.0,155.1/70.0,107.9/70.0,73.51/70.0,48.76/70.0,31.46/70.0,19.69/70.0,12.02/70.0,7.042/70.0,3.974/70.0,2.12/70.0,1.164/70.0,392.4/70.0,98.33/70.0,30.74/70.0,11.51,3.819,0.9345,15.87/70.0,5.502/70.0,1079.1/70.0,754/70.0,0.1525};//Ncoll/70.0for TAA in most places, otherwise we are using offical numbers from cent group if no 70.0is there
  double TAAuncert[nCentBins] = {1.7,1.8,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2.3,5.0,9.5,5,5,5,5,16.};//assume 5% uncert for 'unofficial' values for now (5.0 is official, 5 is not in this list)*/

  TH2D *spec[nTriggers],               *HIspec[HInTriggers][nCentBins][nEvtPlaneBins+1];
  TH1D *evtCount[nTriggers],           *HIevtCount[HInTriggers][nCentBins];
  TH2D *evtCount_JetVars[nTriggers],   *HIevtCount_JetVars[HInTriggers][nCentBins];
  TH1D *nVtxMB,                        *HInVtxMB[nCentBins];
  TH2D *spec_trk[nTriggers_trk],       *HIspec_trk[HInTriggers_trk][nCentBins];
  TH1D *evtCount_trk[nTriggers_trk],   *HIevtCount_trk[HInTriggers_trk][nCentBins];
  TH1D *nVtxMB_trk,                    *HInVtxMB_trk[nCentBins];

  TH1D * pp,                           *HI[nCentBins][nEvtPlaneBins+1];
  TH1D * ppByTrigger[nTriggers],       *HIByTrigger[HInTriggers][nCentBins][nEvtPlaneBins+1];
  TH1D * ppUsedByTrigger[nTriggers],   *HIUsedByTrigger[HInTriggers][nCentBins][nEvtPlaneBins+1];
  TH1D * ppJets,                       *HIJets[nCentBins];
  TH1D * ppJetsByTrigger[nTriggers],   *HIJetsByTrigger[HInTriggers][nCentBins];
  TH1D * pp_perMBTrigger,              *HI_perMBTrigger[nCentBins];               

  TH1D * pp_trk,                              *HI_trk[nCentBins];                                  
  TH1D * ppByTrigger_trk[nTriggers_trk],      *HIByTrigger_trk[HInTriggers_trk][nCentBins];            
  TH1D * ppUsedByTrigger_trk[nTriggers_trk],  *HIUsedByTrigger_trk[HInTriggers_trk][nCentBins];        
  TH1D * ppMaxtrk,                            *HIMaxtrk[nCentBins];                              
  TH1D * ppMaxtrkByTrigger[nTriggers_trk],    *HIMaxtrkByTrigger[HInTriggers_trk][nCentBins];        
  TH1D * pp_perMBTrigger_trk,                 *HI_perMBTrigger_trk[nCentBins];                     

  TH1D * RAA[nCentBins];
  TH1D * RAA_trk[nCentBins];
  TH1D * nTrkOffline[20];

  TH2D *h_scale, *h_scale_trk, *h_HIscale, *h_HIscale_trk;
  TH2D *h_normErr, *h_normErr_trk, *h_HInormErr, *h_HInormErr_trk;
  TH1D *h_normSyst, *h_normSyst_trk, *h_HInormSyst[nCentBins], *h_HInormSyst_trk[nCentBins];
  TH1D *RAA_totSyst[nCentBins], *PbPb_totSyst[nCentBins], *pp_totSyst, *RCP_totSyst[nCentBins]; 
 
  TH1D * h_evtPlanePsi;
  TH1D * h_Q2Mag;

  Settings();
};
 
Settings::Settings()
{
  std::cout << "Getting setting.." << std::endl;
  return;
}

#endif
