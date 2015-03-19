#include "setTDRStyle.C"
#include <TMath.h>

TH1F* readHist(TString nameHist,TString nameFile, int rebin)
{
 TFile* file = new TFile(nameFile);

 TH1F* hist = (TH1F*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 if(rebin>0) hist->Rebin(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TH1F* readHistFromTree(TString nameTree, TString nameHist, TString nameFile, int rebin, int nbins, double xmin, double xmax)
{
  
 TFile* file = new TFile(nameFile);
 TTree *tree = (TTree*)file->Get(nameTree);

 TH1F *h1=new TH1F("h1","",nbins,xmin,xmax);
 TString h2=nameHist+"_h1";
 tree->Draw(nameHist+">> h1");
 h1->Rebin(rebin);

 return h1;
}

TH2D* readHist2D(TString nameHist,TString nameFile, int rebin)
{
 TFile* file = new TFile(nameFile);

 TH2D* hist = (TH2D*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 // if(rebin>0) hist->RebinX(rebin); hist->RebinY(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TH3D* readHist3D(TString nameHist,TString nameFile, int rebin)
{
 TFile* file = new TFile(nameFile);

 TH3D* hist = (TH3D*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 // if(rebin>0) hist->RebinX(rebin); hist->RebinY(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TCanvas* getaCanvas(TString name)
{

  TCanvas* aCanvas = new TCanvas(name,"",92,55,500,500);//,"",181,237,1575,492);

  aCanvas->SetFillColor(0);
  aCanvas->SetBottomMargin(0.125);
  aCanvas->SetLeftMargin(0.125);
  aCanvas->SetFrameFillColor(0);
  aCanvas->SetFrameBorderMode(0);
  aCanvas->SetFrameLineWidth(2);
  return aCanvas;
}

TLegend *legend() {

 TLegend *leg2 = new TLegend(0.52,0.67,0.92,0.90);
 leg2->SetFillStyle(0);
 leg2->SetBorderSize(0);
 leg2->SetTextSize(0.05);
 leg2->SetTextFont(42);

 return leg2;

}

TPad *getaPad_up(TString name){

  TPad *pad1 = new TPad("name", "The pad with the function",0.05,0.4,0.95,0.1);
  //   pad1->Draw();

   //   pad1->Range(-112.6742,-73.17708,1143.438,551.3021);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetGridx();
   pad1->SetGridy();
   pad1->SetLeftMargin(0.1271439);
   pad1->SetRightMargin(0.07307979);
   pad1->SetTopMargin(0.08215179);
   pad1->SetBottomMargin(0.117181);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderMode(0);
  
   return pad1;

}

TPad *getaPad_dn(TString name){

  
 TPad *pad2 = new TPad("pad2", "The pad with the histogram",0.05,0.95,0.95,0.5);
 //   pad2->Draw();

   //   pad2->Range(-153.3652,-2.142584,1185.427,5.367464);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetGridx();
   pad2->SetGridy();
   pad2->SetLeftMargin(0.1215511);
   pad2->SetRightMargin(0.07867263);
   pad2->SetTopMargin(0.04892967);
   pad2->SetBottomMargin(0.1521407);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderMode(0);

   return pad2;

}

TH1D *hTemp(TString ifile, TString dirname, TString hname) {

  TH1D *h=readHist(dirname+hname,ifile,2);

  TH1D *h_temp=h->Clone();
  h_temp->Reset();

  Double_t con, err;

  Int_t nB=h->GetNbinsX();
  for (Int_t i=0; i<nB; i++) {

    con=h->IntegralAndError(i,nB+1,err);

    h_temp->SetBinContent(i,con);
    h_temp->SetBinError(i,err);
  }

  return h_temp;

}

void drawPlot(TString mode, TString cname,TH1F *he, TH1F *hs, TString xname, double ymin, double ymax) {

  if (mode=="histos") TCanvas *c1=getaCanvas(cname);
  // gPad->SetGridx(); gPad->SetGridy();
  
  TLegend *leg = legend();
  //  leg->SetHeader("");
  leg->AddEntry(he,"emulator","L");
  leg->AddEntry(hs,"producer","LF");

  he->SetStats(kTRUE);
  he->Draw("EHIST"); 
  //  he->SetFillColor(4); he->SetFillStyle(3001);
  he->SetLineColor(kBlue-2); he->SetLineWidth(2);
  //  he->SetLineColor(kRed+1); he->SetLineWidth(2);
  hs->Draw("EHISTSAMES");
  hs->SetLineColor(1); hs->SetLineWidth(2);
  //  hs->SetFillColor(1); hs->SetFillStyle(3001);
  he->Draw("EHISTSAMES");

  he->GetYaxis()->SetRangeUser(ymin,ymax);
  he->SetTitle(";"+xname+";a.u.");
  leg->Draw("SAME");

  //  if (mode=="histos") c1->SaveAs(cname+"_uncalib.png");

}


void compare( TString mode="histos", Bool_t debug=false) {
  
  // OR else mode="printfile" or mode="histos"
  Bool_t noPUS=false;
  Bool_t calib=false;

  TString efile, sfile, edata, sdata;

  if (noPUS) { 
    if (calib) { efile="l1t_calib_NoPus.root"; }
    else {efile="l1t_uncalib_NoPus.root";} 
  }  else {
    if (calib) { efile="l1t_calib_withPus.root"; }
    else {efile="l1t_uncalib_withPus.root";} 
  }
  sfile="ttbar_test.root";//"emulator_1k.root";

  edata="emulator";
  sdata="producer";

  TString edirname="l1tStage2CaloAnalyzer/";
  TString sdirname="makeTestTree/tree";
 
  // Stage2 Producer histograms
  TH1F *s_towerPt=readHistFromTree(sdirname,"towerPt",sfile,0,101, -0.5, 100.5);
  TH1F *s_towerEta=readHistFromTree(sdirname,"towerEta",sfile,4,166,-82.5, 82.5); // 83, -41.5, 41.5);
  TH1F *s_towerPhi=readHistFromTree(sdirname,"towerPhi",sfile,4,145, 0.5, 145.5); //73, -0.5, 72.5);

  TH1F *s_towerEtEm=readHistFromTree(sdirname,"towerEtEm",sfile,0,101, -0.5, 100.5);
  TH1F *s_towerEtHad=readHistFromTree(sdirname,"towerEtHad",sfile,0,101, -0.5, 100.5);

  if (noPUS) {

    if (calib) {
      
      TH1F *s_Njet=readHistFromTree(sdirname,"NjetNoPus",sfile,0,10, -0.5, 9.5);
      TH1F *s_jetPt=readHistFromTree(sdirname,"jetNoPusPt",sfile,4,101, -0.5, 100.5);
      TH1F *s_jetEta=readHistFromTree(sdirname,"jetNoPusEta",sfile,4, 166,-82.5, 82.5);//83, -41.5, 41.5);
      TH1F *s_jetPhi=readHistFromTree(sdirname,"jetNoPusPhi",sfile,4,145, 0.5, 145.5); 

    } else {

      TH1F *s_Njet=readHistFromTree(sdirname,"NjetNoPusUncalib",sfile,0,10, -0.5, 9.5);
      TH1F *s_jetPt=readHistFromTree(sdirname,"jetNoPusUncalibPt",sfile,4,101, -0.5, 100.5);
      TH1F *s_jetEta=readHistFromTree(sdirname,"jetNoPusUncalibEta",sfile,4, 166,-82.5, 82.5);//83, -41.5, 41.5);
      TH1F *s_jetPhi=readHistFromTree(sdirname,"jetNoPusUncalibPhi",sfile,4,145, 0.5, 145.5);
    }

    // MHT, HT is always with Calibrated Jets from the Producer
    TH1F *s_mht=readHistFromTree(sdirname,"mhtNoPus",sfile,5,201,-0.5,1000.5);
    TH1F *s_mhtPhi=readHistFromTree(sdirname,"mhtPhiNoPus",sfile,4,145, 0.5, 145.5);
    TH1F *s_ht=readHistFromTree(sdirname,"htNoPus",sfile,10,201,-0.5,1000.5);

  } else {

    if (calib) {

      TH1F *s_Njet=readHistFromTree(sdirname,"Njet",sfile,0,10, -0.5, 9.5);
      TH1F *s_jetPt=readHistFromTree(sdirname,"jetPt",sfile,4,101, -0.5, 100.5);
      TH1F *s_jetEta=readHistFromTree(sdirname,"jetEta",sfile,4, 166,-82.5, 82.5);//83, -41.5, 41.5);
      TH1F *s_jetPhi=readHistFromTree(sdirname,"jetPhi",sfile,4,145, 0.5, 145.5);

    } else {

      TH1F *s_Njet=readHistFromTree(sdirname,"NjetUncalib",sfile,0,10, -0.5, 9.5);
      TH1F *s_jetPt=readHistFromTree(sdirname,"jetUncalibPt",sfile,4,101, -0.5, 100.5);
      TH1F *s_jetEta=readHistFromTree(sdirname,"jetUncalibEta",sfile,4, 166,-82.5, 82.5);//83, -41.5, 41.5);
      TH1F *s_jetPhi=readHistFromTree(sdirname,"jetUncalibPhi",sfile,4,145, 0.5, 145.5);
    }

    TH1F *s_mht=readHistFromTree(sdirname,"mht",sfile,5,201,-0.5,1000.5);
    TH1F *s_mhtPhi=readHistFromTree(sdirname,"mhtPhi",sfile,4,145, 0.5, 145.5);
    TH1F *s_ht=readHistFromTree(sdirname,"ht",sfile,10,201,-0.5,1000.5);

  }

  
  TH1F *s_met=readHistFromTree(sdirname,"met",sfile,5,201,-0.5,1000.5);
  TH1F *s_metPhi=readHistFromTree(sdirname,"metPhi",sfile,4,145, 0.5, 145.5);
  TH1F *s_et=readHistFromTree(sdirname,"et",sfile,10,201,-0.5,1000.5);

  // Emulator histograms

  TH1F *e_towerPt=readHist(edirname+"tower/"+"et",efile,0);
  TH1F *e_towerEta=readHist(edirname+"tower/"+"eta",efile,4);
  TH1F *e_towerPhi=readHist(edirname+"tower/"+"phi",efile,4);
  TH1F *e_towerEtEm=readHist(edirname+"tower/"+"em",efile,0);
  TH1F *e_towerEtHad=readHist(edirname+"tower/"+"had",efile,0);

  TH1F *e_Njet=readHist(edirname+"jet/"+"Njet",efile,0);
  TH1F *e_jetPt=readHist(edirname+"jet/"+"et",efile,4);
  TH1F *e_jetEta=readHist(edirname+"jet/"+"eta",efile,4);
  TH1F *e_jetPhi=readHist(edirname+"jet/"+"phi",efile,4);
  
  TH1F *e_et=readHist(edirname+"sum1/"+"et",efile,10);
  TH1F *e_ht=readHist(edirname+"sum2/"+"et",efile,10);
  
  TH1F *e_met=readHist(edirname+"sum3/"+"et",efile,5);
  TH1F *e_metPhi=readHist(edirname+"sum3/"+"phi",efile,4);
  TH1F *e_mht=readHist(edirname+"sum4/"+"et",efile,5);
  TH1F *e_mhtPhi=readHist(edirname+"sum4/"+"phi",efile,4);

  // Draw comparison plots
  TString fname="validation";

  TString cname="canvas";
  TCanvas *c1=getaCanvas(cname);
  // gPad->SetGridx(); gPad->SetGridy();
  gPad->SetLogy(1);

  // Jet pT/eta/phi
  drawPlot(mode,"nJet",e_Njet,s_Njet,"jet multiplicity",0.1,10000.);
  if (mode=="printfile") c1->Print(fname+"_plots.ps(","test");
   
  //  gPad->SetLogy(1);

   // Towers pT/eta/phi
   drawPlot(mode, "towerPt",e_towerPt,s_towerPt,"tower hwPt (GeV)",1.1,1000000.);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   // drawPlot(mode, "towerEtEm",e_towerEtEm,s_towerEtEm,"tower hwEtEm (GeV)",1.1,100000.);
   // if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   // drawPlot(mode, "towerEtHad",e_towerEtHad,s_towerEtHad,"tower hwEtHad (GeV)",1.1,100000.);
   // if (mode=="printfile") c1->Print(fname+"_plots.ps","test");

   gPad->SetLogy(0);
   drawPlot(mode,"towerEta",e_towerEta,s_towerEta,"tower #eta",1.1,1000000.);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   drawPlot(mode,"towerPhi",e_towerPhi,s_towerPhi,"tower #phi",1.1,1000000.);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");

   //   return;
   drawPlot(mode,"jetPt",e_jetPt,s_jetPt,"jet p_{T} (GeV)",1.1,500.);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   drawPlot(mode,"jetEta",e_jetEta,s_jetEta,"jet #eta",1.1,500.);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   drawPlot(mode,"jetPhi",e_jetPhi,s_jetPhi,"jet #phi",1.1,500.);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");

   //return;

   if (!calib) return;

  // Energy sums
  double sum_max=300.;
  
  drawPlot(mode,"et",e_et,s_et,"sum E_{T} (GeV)",1.1,sum_max);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   drawPlot(mode,"ht",e_ht,s_ht,"H_{T} (GeV)",1.1,sum_max);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   drawPlot(mode,"met",e_met,s_met,"MET (GeV)",1.1,sum_max);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   drawPlot(mode,"metphi",e_metPhi,s_metPhi,"MET #phi",1.1,sum_max);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   drawPlot(mode,"mht",e_mht,s_mht,"MHT (GeV)",1.1,sum_max);
   if (mode=="printfile") c1->Print(fname+"_plots.ps","test");
   drawPlot(mode,"mhtphi",e_mhtPhi,s_mhtPhi,"MHT #phi",1.1,sum_max);
   if (mode=="printfile") c1->Print(fname+"_plots.ps)","test");

}
