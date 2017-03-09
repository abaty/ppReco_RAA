#include "Settings.h"
#include "TFrame.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
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

void makePlotsPbPb(Settings s)
{
  for(int c = 0; c<s.nCentBins; c++){
    for(int j = 0; j<s.HInTriggers; j++){
      s.HIByTrigger[j][c]->SetLineColor(j+1);
      s.HIByTrigger[j][c]->SetLineWidth(1);
      s.HIByTrigger[j][c]->SetMarkerColor(j+1);
      s.HIByTrigger[j][c]->SetFillColor(j+1);
      s.HIUsedByTrigger[j][c]->SetLineColor(kBlack);
      s.HIUsedByTrigger[j][c]->SetLineWidth(2);
      s.HIUsedByTrigger[j][c]->SetMarkerColor(j+1);
      s.HIUsedByTrigger[j][c]->SetMarkerSize(0);
      s.HIUsedByTrigger[j][c]->SetFillColor(j+1);
      s.HIJetsByTrigger[j][c]->SetLineColor(j+1);
      s.HIJetsByTrigger[j][c]->SetLineWidth(1);
      s.HIJetsByTrigger[j][c]->SetMarkerColor(j+1);
      s.HIJetsByTrigger[j][c]->SetMarkerSize(0.8);
      s.HIJetsByTrigger[j][c]->SetFillColor(j+1);
      if(j==4){
        s.HIJetsByTrigger[j][c]->SetMarkerColor(kOrange-3);
        s.HIJetsByTrigger[j][c]->SetLineColor(kOrange-3);
      }
    }
  }
  
//******************************************************************************************************************
//************************************************JET TRIGGER PLOTS**********************************************************
//******************************************************************************************************************
  std::cout << "jet plotting" << std::endl;
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  TLine * l[8];
  TLatex * lat = new TLatex(1,1,"test");
  TLegend * leg;
  c1->SetLogy();
  for(int c = 0; c<s.nCentBins; c++){
    if(s.lowCentBin[c]*5<50) continue; 
    gStyle->SetLegendBorderSize(0); 
    leg = new TLegend(0.5,0.6,0.9,0.9);
    s.HIJets[c]->Scale(100);
    s.HIJets[c]->SetMarkerSize(0);
    float Ymin = 0.00000000001;
    float Ymax = 10;
    c1->Clear();
    s.HIJets[c]->GetYaxis()->SetRangeUser(Ymin,Ymax);
    s.HIJets[c]->Draw("h");
    for(int i = 0; i<s.HInTriggers; i++){
      s.HIJetsByTrigger[i][c]->Draw("same");
    }
    std::cout << "jet plotting" << std::endl;
    leg->Clear();
    leg->AddEntry(s.HIJets[c],"Jet Spectrum (x100)","l");
    leg->AddEntry(s.HIJetsByTrigger[0][c],"MB (I)","p");
    leg->AddEntry(s.HIJetsByTrigger[1][c],"jet 40 (II)","p");
    leg->AddEntry(s.HIJetsByTrigger[2][c],"jet 60 (III)","p");
    leg->AddEntry(s.HIJetsByTrigger[3][c],"jet 80 (IV)","p");
    leg->AddEntry(s.HIJetsByTrigger[4][c],"jet 100 (V)","p");
    leg->AddEntry((TObject*)0,"akPu4Calo Jets, |#eta|<2","");
    leg->AddEntry((TObject*)0,Form("%d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    std::cout << "jet plotting" << std::endl;
    c1->SaveAs(Form("plots/png/HIJets_FullSpectrum_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c1->SaveAs(Form("plots/pdf/HIJets_FullSpectrum_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c1->Clear();
    std::cout << "jet plotting" << std::endl;
    
    leg->SetX1NDC(0.5);
    leg->SetX2NDC(0.8);
    s.HIJets[c]->GetXaxis()->SetRangeUser(20,200);
    Ymin = 0.00000001;
    Ymax = 100;
    s.HIJets[c]->GetYaxis()->SetRangeUser(Ymin,Ymax);
    s.HIJets[c]->Draw("h");
    for(int i = 0; i<s.HInTriggers; i++) s.HIJetsByTrigger[i][c]->Draw("same");
    leg->Draw("same"); 
    for(int i = 0; i<4; i++) 
    {
      l[i] = new TLine(s.HItriggerBins[i+1],Ymin,s.HItriggerBins[i+1],Ymax); 
      l[i]->SetLineWidth(2);
      l[i]->SetLineStyle(2);
      l[i]->SetLineColor(1);
      l[i]->Draw("same");
    }
    lat->DrawLatex(45,Ymin*3,"I");
    lat->DrawLatex(65,Ymin*3,"II");
    lat->DrawLatex(85,Ymin*3,"III");
    lat->DrawLatex(105,Ymin*3,"IV");
    lat->DrawLatex(125,Ymin*3,"V");
    c1->SaveAs(Form("plots/png/HIJets_FullSpectrum_XZoom_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c1->SaveAs(Form("plots/pdf/HIJets_FullSpectrum_XZoom_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->Clear();
    for(int i = 0; i<4; i++) delete l[i]; 
    delete leg;
  }
  for(int c = 0; c<s.nCentBins;c++){
    if(s.lowCentBin[c]*5<50) continue; 
    c1->SetLogy(0);
    for(int i = 0; i<s.HInTriggers-1; i++) s.HIJetsByTrigger[s.HInTriggers-1-i][c]->Divide(s.HIJetsByTrigger[s.HInTriggers-2-i][c]);
    s.HIJetsByTrigger[1][c]->GetYaxis()->SetRangeUser(0,2);
    s.HIJetsByTrigger[2][c]->GetYaxis()->SetRangeUser(0,2);
    s.HIJetsByTrigger[1][c]->GetXaxis()->SetRangeUser(20,140);
    s.HIJetsByTrigger[2][c]->GetXaxis()->SetRangeUser(20,140);
    s.HIJetsByTrigger[1][c]->GetXaxis()->SetTitle("Leading jet p_{T}");
    s.HIJetsByTrigger[2][c]->GetXaxis()->SetTitle("Leading jet p_{T}");
    s.HIJetsByTrigger[1][c]->Draw();
    s.HIJetsByTrigger[2][c]->Draw("same");
    if(!s.doBetterHITrig && s.lowCentBin[c]<10) s.HIJetsByTrigger[3][c]->Draw("same");
    if(!s.doBetterHITrig && s.lowCentBin[c]<6)  s.HIJetsByTrigger[4][c]->Draw("same");
    gStyle->SetLegendBorderSize(0); 
    leg = new TLegend(0.1779449,0.652819,0.6077694,0.8916914);
    leg->AddEntry((TObject*)0,"Anti-K_{T} Jets, |#eta|<2","");
    if(!s.doBetterHITrig || s.lowCentBin[c]>5){ 
      leg->AddEntry(s.HIJetsByTrigger[1][c],"Jet40/Minimum Bias","p");
      leg->AddEntry(s.HIJetsByTrigger[2][c],"Jet60/Jet40","p");
    }else{
      leg->AddEntry(s.HIJetsByTrigger[2][c],"Jet60/Minimum Bias","p");
    }
    if(!s.doBetterHITrig || s.lowCentBin[c]<10)  leg->AddEntry(s.HIJetsByTrigger[3][c],"Jet80/Jet60","p");
    if(!s.doBetterHITrig || s.lowCentBin[c]<6)   leg->AddEntry(s.HIJetsByTrigger[4][c],"Jet100/Jet80","p");
    leg->AddEntry((TObject*)0,Form("%d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c1->SaveAs(Form("plots/png/HI_JetRelativeTurnOnes_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/HI_JetRelativeTurnOnes_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->Clear();
    //start pretty plot
    setTDRStyle();
    int W = 800;
    int H = 700;//700
    int H_ref = 700;//700
    int W_ref = 800;
    float T = 0.08*H_ref;
    float B = 0.12*H_ref; 
    float L = 0.15*W_ref;
    float R = 0.04*W_ref;
    
    TCanvas* canv = new TCanvas("RAA","RAA",50,50,W,H);
    canv->SetFillColor(0);
    canv->SetBorderMode(0);
    canv->SetFrameFillStyle(0);
    canv->SetFrameBorderMode(0);
    canv->SetLeftMargin( L/W );
    canv->SetRightMargin( R/W );
    canv->SetTopMargin( T/H );
    canv->SetBottomMargin( B/H );
    canv->SetTickx(0);
    canv->SetTicky(0);
    for(int i=1; i<5; i++) s.HIJetsByTrigger[i][c]->SetMarkerSize(1);
    s.HIJetsByTrigger[2][c]->GetXaxis()->CenterTitle();
    s.HIJetsByTrigger[2][c]->GetXaxis()->SetTitle("Offline Leading Jet p_{T} (GeV)");
    s.HIJetsByTrigger[2][c]->GetXaxis()->SetLabelSize(0.04);
    s.HIJetsByTrigger[2][c]->GetXaxis()->SetTitleSize(0.05);
    s.HIJetsByTrigger[2][c]->GetYaxis()->CenterTitle();
    s.HIJetsByTrigger[2][c]->GetYaxis()->SetLabelSize(0.04);
    s.HIJetsByTrigger[2][c]->GetYaxis()->SetTitleSize(0.05);
    s.HIJetsByTrigger[2][c]->GetYaxis()->SetTitleOffset(1.3);
    s.HIJetsByTrigger[2][c]->GetYaxis()->SetTitle("Triggered Spectrum Ratio");
    gStyle->SetErrorX(0.);
    s.HIJetsByTrigger[2][c]->Draw("p");
    TLine * line1 = new TLine(20,1,140,1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw("same");
    if(!s.doBetterHITrig || s.lowCentBin[c]>5) s.HIJetsByTrigger[1][c]->Draw("same p");
    s.HIJetsByTrigger[2][c]->Draw("same p");
    if(!s.doBetterHITrig || s.lowCentBin[c]<10) s.HIJetsByTrigger[3][c]->Draw("same p");
    if(!s.doBetterHITrig || s.lowCentBin[c]<6)  s.HIJetsByTrigger[4][c]->Draw("same p");
    s.HIJetsByTrigger[2][c]->Draw("sameaxis");
    s.HIJetsByTrigger[2][c]->Draw("same");
    if(!s.doBetterHITrig || s.lowCentBin[c]>5) s.HIJetsByTrigger[1][c]->Draw("same");
    leg->Draw("same");
    
    int iPeriod = 0;
    lumi_sqrtS = "404 #mub^{-1} (5.02 TeV PbPb)";
    writeExtraText = false;  
    extraText  = "Preliminary";
    //extraText  = "Unpublished";
    CMS_lumi( canv, iPeriod, 33 );
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();   
    canv->SaveAs(Form("plots/png/HI_PrettyJetRelativeTurnOnes_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    canv->SaveAs(Form("plots/pdf/HI_PrettyJetRelativeTurnOnes_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    canv->SaveAs(Form("plots/png/HI_PrettyJetRelativeTurnOnes_%d_%d.C",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete canv;
    //end pretty plot
    delete leg;
    delete line1;
  }
  
  c1->SetLogy();
  c1->SetLogx();
  for(int c = 0; c<s.nCentBins;c++){
    if(s.lowCentBin[c]*5<50) continue; 
    s.HI[c]->SetMarkerSize(0.8);
    s.HI[c]->GetYaxis()->SetRangeUser(TMath::Max(s.HI[20]->GetMinimum()/200.0,1e-15),s.HI[20]->GetMaximum()*1000);
    s.HI[c]->GetXaxis()->SetRangeUser(0.7,400);
    s.HI[c]->Draw();
    s.HIUsedByTrigger[0][c]->SetFillColor(kGray);
    //s.HIUsedByTrigger[4][c]->SetFillColor(kCyan+2);
    s.HIUsedByTrigger[4][c]->SetFillColor(kCyan-10);
    s.HIUsedByTrigger[3][c]->Add(s.HIUsedByTrigger[4][c]);
    s.HIUsedByTrigger[3][c]->SetFillColor(kCyan+2);
    s.HIUsedByTrigger[2][c]->Add(s.HIUsedByTrigger[3][c]);
    s.HIUsedByTrigger[1][c]->Add(s.HIUsedByTrigger[2][c]);
    s.HIUsedByTrigger[0][c]->Add(s.HIUsedByTrigger[1][c]);
    /*if(s.doBetterHITrig && s.lowCentBin[c]<6){  s.HIUsedByTrigger[1][c]->SetFillColor(kGray); s.HIUsedByTrigger[1][c]->SetLineWidth(-1);}
    if(s.doBetterHITrig && s.lowCentBin[c]>5){  s.HIUsedByTrigger[4][c]->SetFillColor(kCyan+2);s.HIUsedByTrigger[3][c]->SetLineWidth(-1);}
    if(s.doBetterHITrig && s.lowCentBin[c]>9){  s.HIUsedByTrigger[3][c]->SetFillColor(2); s.HIUsedByTrigger[4][c]->SetFillColor(2); s.HIUsedByTrigger[2][c]->SetLineWidth(0); s.HIUsedByTrigger[3][c]->SetLineWidth(0);}*/
    for(int i = 0; i<s.HInTriggers; i++){
      if(s.doBetterHITrig && s.lowCentBin[c]<6 && i==1) continue;
      if(s.doBetterHITrig && s.lowCentBin[c]>5 && i==4) continue;
      if(s.doBetterHITrig && s.lowCentBin[c]>9 && (i==3 || i==4)) continue;
      s.HIUsedByTrigger[i][c]->Draw("HIST same");
    }
    s.HI[c]->Draw("sameaxis");
    s.HI[c]->Draw("same");
    leg = new TLegend(0.51,0.61,0.95,0.9,NULL,"brNDC");
    leg->AddEntry(s.HI[c],"PbPb Uncorrected Spectrum","p");
    leg->AddEntry(s.HIUsedByTrigger[0][c],"Minimum Bias","f");
    if(!s.doBetterHITrig || s.lowCentBin[c]>5) leg->AddEntry(s.HIUsedByTrigger[1][c],"Jet40 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger[2][c],"Jet60 trigger","f");
    if(!s.doBetterHITrig && s.lowCentBin[c]<10) leg->AddEntry(s.HIUsedByTrigger[3][c],"Jet80 trigger","f");
    if(!s.doBetterHITrig && s.lowCentBin[c]<6) leg->AddEntry(s.HIUsedByTrigger[4][c],"Jet100 trigger","f");
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c1->SaveAs(Form("plots/png/HITrack_FullSpectrum_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/HITrack_FullSpectrum_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));

    //start pretty plot
    setTDRStyle();
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
    canv->SetTickx(0);
    canv->SetTicky(0);
    s.HI[c]->GetXaxis()->SetRangeUser(0.7,350);
    s.HI[c]->GetYaxis()->SetRangeUser(1e-14,999);
    s.HI[c]->GetYaxis()->SetNdivisions(508);
    s.HI[c]->SetMarkerSize(1.2);
    s.HI[c]->GetXaxis()->CenterTitle();
    s.HI[c]->GetXaxis()->SetTitle("p_{T} (GeV)");
    s.HI[c]->GetXaxis()->SetLabelSize(0.04);
    s.HI[c]->GetYaxis()->CenterTitle();
    s.HI[c]->GetYaxis()->SetLabelSize(0.04);
    s.HI[c]->GetYaxis()->SetTitleSize(0.04);
    s.HI[c]->GetYaxis()->SetTitleOffset(1.6);
    s.HI[c]->GetYaxis()->SetTitle("#frac{1}{N_{evt}} E#frac{d^{3}N_{trk}}{dp^{3}} (GeV)^{-2}");
    s.HI[c]->Draw();
    for(int i = 0; i<s.HInTriggers; i++){
      if(s.doBetterHITrig && s.lowCentBin[c]<6 && i==1) continue;
      if(s.doBetterHITrig && s.lowCentBin[c]>5 && i==4) continue;
      if(s.doBetterHITrig && s.lowCentBin[c]>9 && (i==3 || i==4)) continue;
      s.HIUsedByTrigger[i][c]->Draw("HIST same");
    }
    s.HI[c]->Draw("sameaxis");
    s.HI[c]->Draw("same");
    leg->SetTextSize(0.03);
    leg->Draw("same");
    
    int iPeriod = 0;
    lumi_sqrtS = "404 #mub^{-1} (5.02 TeV PbPb)";
    writeExtraText = false;  
    extraText  = "Preliminary";
    //extraText  = "Unpublished";
    CMS_lumi( canv, iPeriod, 11, true );
    canv->SetLogy();
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();   
    canv->SaveAs(Form("plots/png/HITrack_PrettyFullSpectrum_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    canv->SaveAs(Form("plots/pdf/HITrack_PrettyFullSpectrum_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    canv->SaveAs(Form("plots/png/HITrack_PrettyFullSpectrum_%d_%d.C",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    delete canv;
    //end pretty plot
   
    c1->SetLogy(0);
    for(int i = 0; i<s.HInTriggers; i++) s.HIUsedByTrigger[i][c]->Divide(s.HI[c]);
    s.HIUsedByTrigger[0][c]->GetYaxis()->SetRangeUser(0,2);
    s.HIUsedByTrigger[0][c]->GetYaxis()->SetTitle("Relative Contribution to Bin");
    s.HIUsedByTrigger[0][c]->Draw("HIST");
    for(int i = 1; i<s.HInTriggers; i++){
      if(s.doBetterHITrig && s.lowCentBin[c]<6 && i==1) continue;
      if(s.doBetterHITrig && s.lowCentBin[c]>5 && i==4) continue;
      if(s.doBetterHITrig && s.lowCentBin[c]>9 && (i==3 || i==4)) continue;
      s.HIUsedByTrigger[i][c]->Draw("HIST same");
    }
    s.HIUsedByTrigger[0][c]->Draw("sameaxis");
    leg->Draw("same");
    c1->SaveAs(Form("plots/png/HITrack_FullSpectrum_%d_%d_relativeContribution.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/HITrack_FullSpectrum_%d_%d_relativeContribution.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SetLogy();
    delete leg;
  }

  c1->SetLogy(0);
  for(int c = 0; c<s.nCentBins;c++){
    if(s.lowCentBin[c]*5<50) continue; 
    c1->Clear();
    s.RAA[c] = (TH1D*)s.HI[c]->Clone(Form("RAA_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    //s.RAA[c]->Scale(1/s.nColl[c]);//for not using lumi
    //s.RAA[c]->Divide(s.pp);//for not using lumi
    s.RAA[c]->Scale(70.0/s.nColl[c]);//TAA
    s.RAA[c]->Divide(s.pp_perMBTrigger);//for not using lumi
    s.RAA[c]->Write();

    s.RAA[c]->GetYaxis()->SetRangeUser(0,1.2); 
    s.RAA[c]->GetYaxis()->SetTitle("R_{AA}"); 
    s.RAA[c]->Draw();
    leg = new TLegend(0.2,0.7,0.4,0.9);
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->AddEntry((TObject*)0,"#sqrt{s_{NN}} = 5.02 TeV","");
    leg->Draw("same");
    TLine *line = new TLine(0.5,1,400,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("same");
    c1->SaveAs(Form("plots/png/RAA_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/RAA_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
    delete line;
  }
 
  /*TFile * KKRAA = TFile::Open("Krajczar_RAA_070316.root","read");
  TH1D * KKRAA_h[21];
  for(int c = 0; c<21; c++){
    if(c<20) KKRAA_h[c] = (TH1D*)KKRAA->Get(Form("hForPlotting_%d_%d",c*5,5+c*5));
    else     KKRAA_h[c] = (TH1D*)KKRAA->Get("hTrackPt_trkCorr_PbPb_copy1");
    c1->Clear();
    s.RAA[c]->Draw();
    KKRAA_h[c]->SetMarkerStyle(25); 
    KKRAA_h[c]->Draw("same"); 
    leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->AddEntry(s.RAA[c],"AB's RAA","p");
    leg->AddEntry(KKRAA_h[c],"KK's RAA","p");
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->AddEntry((TObject*)0,"#sqrt{s_{NN}} = 5.02 TeV","");
    leg->Draw("same");
    TLine *line = new TLine(0.5,1,400,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("same");
    c1->SaveAs(Form("plots/png/RAA_Comparison_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/RAA_Comparison_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
    delete line;
  }
  KKRAA->Close();
  */
 
  delete lat;
  delete c1;
  
//************************************************************************************************************
//***********************************************TRACK TRIGGER PLOTS******************************************
//************************************************************************************************************
  TCanvas * c2 = new TCanvas("c2","c2",800,600);
/*  for(int c = 0; c<s.nCentBins; c++){
    c2->SetLogy();
    c2->Clear();
    s.HIMaxtrk[c]->SetMarkerSize(0);
    s.HIMaxtrkByTrigger[4][c]->SetMarkerColor(kOrange);
    s.HIMaxtrkByTrigger[4][c]->SetLineColor(kOrange);
    s.HIMaxtrk[c]->Scale(100);
    float Ymin = 0.000000001;
    float Ymax = 1000;
    s.HIMaxtrk[c]->GetYaxis()->SetRangeUser(Ymin,Ymax);
    s.HIMaxtrk[c]->Draw("h");
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIMaxtrkByTrigger[i][c]->Draw("same");
    leg = new TLegend(0.5,0.6,0.8,0.9);
    leg->AddEntry(s.HIMaxtrk[c],"Leading Trk p_{T}Spectrum (x100)","l");
    leg->AddEntry(s.HIMaxtrkByTrigger[0][c],"MB (I)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[1][c],"Track 12 (II)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[2][c],"Track 18 (III)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[3][c],"Track 24 (IV)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[4][c],"Track 34 (V)","p");
    leg->AddEntry((TObject*)0,Form("%d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c2->SaveAs(Form("plots/png/HIMaxtrk_FullSpectrum_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c2->SaveAs(Form("plots/pdf/HIMaxtrk_FullSpectrum_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c2->Clear();
 
    s.HIMaxtrk[c]->GetXaxis()->SetRangeUser(0,70);
    Ymin = 0.00000001;
    Ymax = 1000;
    s.HIMaxtrk[c]->GetYaxis()->SetRangeUser(Ymin,Ymax);
    s.HIMaxtrk[c]->Draw("h");
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIMaxtrkByTrigger[i][c]->Draw("same");
    leg->Draw("same"); 
    for(int i = 0; i<5; i++) 
    {
      l[i] = new TLine(s.HItriggerBins_trk[i+1],Ymin,s.HItriggerBins_trk[i+1],Ymax); 
      l[i]->SetLineWidth(2);
      l[i]->SetLineStyle(2);
      l[i]->SetLineColor(1);
      l[i]->Draw("same");
    }
    lat->DrawLatex(8,Ymin*3,"I");
    lat->DrawLatex(16,Ymin*3,"II");
    lat->DrawLatex(22,Ymin*3,"III");
    lat->DrawLatex(28,Ymin*3,"IV");
    lat->DrawLatex(38,Ymin*3,"V");
    c2->SaveAs(Form("plots/png/HIMaxtrk_FullSpectrum_XZoom_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c2->SaveAs(Form("plots/pdf/HIMaxtrk_FullSpectrum_XZoom_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    for(int i = 0; i<5; i++) delete l[i]; 
    delete leg;
  
    c2->SetLogy(0);
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIMaxtrkByTrigger[i][c]->Rebin(2);
    for(int i = 0; i<s.HInTriggers_trk-1; i++) s.HIMaxtrkByTrigger[s.HInTriggers_trk-1-i][c]->Divide(s.HIMaxtrkByTrigger[s.HInTriggers_trk-2-i][c]);
    s.HIMaxtrkByTrigger[1][c]->GetYaxis()->SetRangeUser(0,2);
    s.HIMaxtrkByTrigger[1][c]->GetXaxis()->SetRangeUser(10,90);
    s.HIMaxtrkByTrigger[1][c]->GetXaxis()->SetTitle("Leading Track p_{T}");
    s.HIMaxtrkByTrigger[1][c]->Draw();
    s.HIMaxtrkByTrigger[2][c]->Draw("same");
    s.HIMaxtrkByTrigger[3][c]->Draw("same");
    s.HIMaxtrkByTrigger[4][c]->Draw("same");
    leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->AddEntry(s.HIMaxtrkByTrigger[1][c],"Trk12/MB","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[2][c],"Trk18/Trk12","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[3][c],"Trk24/Trk18","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[4][c],"Trk34/Trk24","p");
    leg->AddEntry((TObject*)0,Form("%d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c2->SaveAs(Form("plots/png/HITrackRelativeTurnOnes_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SaveAs(Form("plots/pdf/HITrackRelativeTurnOnes_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
  }

  c2->SetLogy();
  c2->SetLogx();
  for(int c = 0; c<s.nCentBins; c++){
    c2->Clear();
    s.HI_trk[c]->SetMarkerSize(0.8);
    s.HI_trk[c]->GetYaxis()->SetRangeUser(TMath::Max(s.HI_trk[20]->GetMinimum()/200.0,1e-15),s.HI_trk[20]->GetMaximum()*10);
    s.HI_trk[c]->GetXaxis()->SetRangeUser(0.7,400);
    s.HI_trk[c]->Draw();
    s.HIUsedByTrigger_trk[0][c]->SetFillColor(kGray);
    s.HIUsedByTrigger_trk[3][c]->SetFillColor(kBlue);
    s.HIUsedByTrigger_trk[4][c]->SetFillColor(kCyan+2);
    s.HIUsedByTrigger_trk[3][c]->Add(s.HIUsedByTrigger_trk[4][c]);
    s.HIUsedByTrigger_trk[2][c]->Add(s.HIUsedByTrigger_trk[3][c]);
    s.HIUsedByTrigger_trk[1][c]->Add(s.HIUsedByTrigger_trk[2][c]);
    s.HIUsedByTrigger_trk[0][c]->Add(s.HIUsedByTrigger_trk[1][c]);
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIUsedByTrigger_trk[i][c]->Draw("HIST same");
    leg = new TLegend(0.5,0.6,0.8,0.9);
    s.HI_trk[c]->Draw("sameaxis");
    s.HI_trk[c]->Draw("same");
    leg->Clear();
    leg->AddEntry(s.HI_trk[c],"PbPb track Spectrum","p");
    leg->AddEntry(s.HIUsedByTrigger_trk[0][c],"MB trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[1][c],"Track12 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[2][c],"Track18 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[3][c],"Track24 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[4][c],"Track34 trigger","f");
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
     
    c2->SaveAs(Form("plots/png/HITrack_FullSpectrum_trk_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SaveAs(Form("plots/pdf/HITrack_FullSpectrum_trk_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SetLogy(0);
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIUsedByTrigger_trk[i][c]->Divide(s.HI_trk[c]);
    s.HIUsedByTrigger_trk[0][c]->GetYaxis()->SetRangeUser(0,2);
    s.HIUsedByTrigger_trk[0][c]->GetYaxis()->SetTitle("Relative Contribution to Bin");
    s.HIUsedByTrigger_trk[0][c]->Draw("HIST");
    for(int i = 1; i<s.HInTriggers_trk; i++) s.HIUsedByTrigger_trk[i][c]->Draw("HIST same");
    s.HIUsedByTrigger_trk[0][c]->Draw("sameaxis");
    leg->Draw("same");
    c2->SaveAs(Form("plots/png/HITrack_FullSpectrum_trk_%d_%d_relativeContribution.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SaveAs(Form("plots/pdf/HITrack_FullSpectrum_trk_%d_%d_relativeContribution.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SetLogy();
    delete leg;
  }
  
  c2->SetLogy(0);
  for(int c = 0; c<s.nCentBins;c++){
    c2->Clear();
    s.RAA_trk[c] = (TH1D*)s.HI_trk[c]->Clone(Form("RAA_trk_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    //s.RAA_trk[c]->Scale(1/s.nColl[c]);//for using lumi
    //s.RAA_trk[c]->Divide(s.pp_trk); //for using lumi
    s.RAA_trk[c]->Scale(1.0/s.TAA[c]);//TAA
    s.RAA_trk[c]->Divide(s.pp_perMBTrigger_trk);//for not using lumi
    //s.RAA_trk[c]->Write();

    s.RAA_trk[c]->GetYaxis()->SetRangeUser(0,1.2); 
    s.RAA_trk[c]->GetYaxis()->SetTitle("R_{AA}"); 
    s.RAA_trk[c]->SetMarkerStyle(25); 
    s.RAA_trk[c]->Draw();
    s.RAA[c]->Draw("same");
    leg = new TLegend(0.2,0.7,0.4,0.9);
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->AddEntry((TObject*)0,"#sqrt{s_{NN}} = 5.02 TeV","");
    leg->AddEntry(s.RAA[c],"Jet Triggers","p");
    leg->AddEntry(s.RAA_trk[c],"Track Triggers","p");
    leg->Draw("same");
    TLine *line = new TLine(0.5,1,400,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("same");
    c1->SaveAs(Form("plots/png/RAA_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/RAA_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
    delete line;
  }

  c2->SetLogy(0);
  TH1D * HIjtVsTrk[s.nCentBins];
  for(int c = 0; c<s.nCentBins;c++){
    HIjtVsTrk[c] = (TH1D*)s.HI[c]->Clone(Form("HIJetVsTrk%d",c));
    HIjtVsTrk[c]->Divide(s.HI_trk[c]);
    HIjtVsTrk[c]->GetYaxis()->SetTitle("PbPb Jet triggers/PbPb Track triggers");
    HIjtVsTrk[c]->GetYaxis()->SetRangeUser(0.5,1.5);
    HIjtVsTrk[c]->Draw();
    leg = new TLegend(0.2,0.7,0.4,0.9);
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c2->SaveAs(Form("plots/png/HI_jetVsTrkTriggers_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SaveAs(Form("plots/pdf/HI_JetVsTrkTriggers_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    delete leg;
  }
*/  
  c2->SetLogx();
  for(int c = 0; c<s.nCentBins; c++)
  {
    if(s.lowCentBin[c]*5<50) continue; 
    s.h_HInormSyst[c]->GetXaxis()->SetRangeUser(0.7,400);
    s.h_HInormSyst[c]->GetXaxis()->SetTitle("p_{T}");
    s.h_HInormSyst[c]->GetYaxis()->SetRangeUser(0,0.04);
    s.h_HInormSyst[c]->GetYaxis()->SetTitle("Normalization Uncertainty");
    leg = new TLegend(0.2,0.7,0.4,0.9);
    leg->AddEntry((TObject*)0,Form("%d_%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->AddEntry((TObject*)0,"Jet Triggers","");
    s.h_HInormSyst[c]->Draw();
    //s.h_HInormSyst[c]->Print("All");
    leg->Draw("same");
    c2->SaveAs(Form("plots/png/HI_normalizationError_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SaveAs(Form("plots/pdf/HI_normalizationError_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
  }
  s.h_normSyst->GetXaxis()->SetRangeUser(0.7,400);
  s.h_normSyst->GetXaxis()->SetTitle("p_{T}");
  s.h_normSyst->GetYaxis()->SetRangeUser(0,0.04);
  s.h_normSyst->GetYaxis()->SetTitle("Normalization Uncertainty");
  leg = new TLegend(0.2,0.7,0.4,0.9);
  leg->AddEntry(s.h_normSyst,"pp Jet Triggers","");
  s.h_normSyst->Draw();
  leg->Draw("same");
  c2->SaveAs("plots/png/pp_normalizationError.png");
  c2->SaveAs("plots/pdf/pp_normalizationError.pdf");
  delete leg;

  delete c2;

  return;
}
