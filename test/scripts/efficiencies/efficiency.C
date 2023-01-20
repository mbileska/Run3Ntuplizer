#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <new>
#include<climits>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"

#include "TStyle.h"
#include "TAttFill.h"
#include "TPaveStats.h"
#include "TMinuit.h"
#include "TPostScript.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"

using namespace std;

void efficiency(){
TCanvas *c1 =new TCanvas("c1", " ", 0, 0,600,700);
c1->Range(0,0,1,1);
c1->SetFillColor(0);
c1->SetBorderMode(0);
c1->SetBorderSize(2);
c1->SetFrameBorderMode(0);
c1->Draw();
c1->SetGrid();
gStyle->SetOptStat(0);

TFile *g1 =TFile::Open("out_version1.root");
TH1F *h1 = (TH1F*)g1->Get("recojetpt_eff_den1");
h1->Rebin();
TH1F *h2 = (TH1F*)g1->Get("recojetpt_eff_num1");
h2->Rebin();

TFile *g2 =TFile::Open("out_version2.root");
TH1F *h3 = (TH1F*)g2->Get("recojetpt_eff_den1");
h3->Rebin();
TH1F *h4 = (TH1F*)g2->Get("recojetpt_eff_num1");
h4->Rebin();

h1->SetLineColorAlpha(kWhite, 1.);
h1->SetTitle("");
h1->GetXaxis()->SetTitle("Offline p_{T} [GeV]");
h1->GetXaxis()->SetTitleSize(0.045);
h1->GetXaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitle("Efficiency");
h1->GetYaxis()->SetTitleSize(0.045);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetRangeUser(0.,1.2);
h1->Draw();

c1->Update();
c1->RedrawAxis();
TLine l;
l.DrawLine(c1->GetUxmin(), c1->GetUymax(), c1->GetUxmax(), c1->GetUymax());
l.DrawLine(c1->GetUxmax(), c1->GetUymin(), c1->GetUxmax(), c1->GetUymax());

TEfficiency* pEff1 = 0;

if(TEfficiency::CheckConsistency(*h2,*h1))
{
  pEff1 = new TEfficiency(*h2,*h1);
  pEff1->SetLineWidth(3.);
  pEff1->SetLineColor(kBlack);
  pEff1->Draw("same");
}

TEfficiency* pEff2 = 0;

if(TEfficiency::CheckConsistency(*h4,*h3))
{
  pEff2 = new TEfficiency(*h4,*h3);
  pEff2->SetLineWidth(3.);
  pEff2->SetLineColor(kRed);
  pEff2->Draw("same");
}

TLegend *legend1 = new TLegend(0.45, 0.3, 0.8, 0.45);
legend1->SetTextFont(42);
legend1->SetLineColor(0);
legend1->SetTextSize(0.045);
legend1->SetFillColor(0);
legend1->AddEntry(pEff1, "emulator v1", "l");
legend1->AddEntry(pEff2, "emulator v2", "l");
legend1->Draw("same");

TLatex *t1 = new TLatex(0.65,0.45," L1 boosted p_{T} > 120 GeV");
t1->SetNDC();
t1->SetTextFont(42);
t1->SetTextSize(0.045);
t1->SetTextAlign(20);
t1->Draw("same");

TLatex *t2 = new TLatex(0.5,0.9," #bf{CMS} #it{Simulation Preliminary}         (13 TeV)");
t2->SetNDC();
t2->SetTextFont(42);
t2->SetTextSize(0.045);
t2->SetTextAlign(20);
t2->Draw("same");

TLatex *t3 = new TLatex(0.7,0.8," #bf{Boosted SM H #rightarrow bb} ");
t3->SetNDC();
t3->SetTextFont(42);
t3->SetTextSize(0.035);
t3->SetTextAlign(20);
t3->Draw("same");

c1->SaveAs("efficiencies.pdf");
}
