# Usage: ./make_eff_subleadingjet.sh sample date
sample=$1
date=$2
echo "Making efficiency plots..."
if [ -d "${date}" ]; then
        echo "Directory exists..."
else
        echo "Making directory..."
        mkdir ${date}
fi

for etarange in barrel endcap forward all; do

if [ ${etarange} == 'barrel' ]; then 
	rangeprint="|#eta^{off}| < 1.474"
elif [ ${etarange} == 'endcap' ]; then
	rangeprint="1.474 < |#eta^{off}| < 3"
elif [ ${etarange} == 'forward' ]; then
	rangeprint="3 < |#eta^{off}| < 5"
elif [ ${etarange} == 'all' ]; then
        rangeprint=""
fi

for num in recojetpt_eff_num2 recojetpt_eff_num2_s2; do

if [ ${num} == 'recojetpt_eff_num2' ]; then
        labelprint="Stage-3"
elif [ ${num} == 'recojetpt_eff_num2_s2' ]; then
        labelprint="Stage-2"
fi
	
cat>dummy.C<<EOF

{
TCanvas *c1 =new TCanvas("c1", " ", 0, 0,700,1000);
c1->Range(0,0,1,1);
c1->SetFillColor(0);
c1->SetBorderMode(0);
c1->SetBorderSize(2);
c1->SetFrameBorderMode(0);
c1->Draw();
c1->SetGrid();
gStyle->SetOptStat(0);

TFile *g1 =TFile::Open("plots_${sample}_${date}_l1_35_${etarange}.root");
TH1F *h1 = (TH1F*)g1->Get("recojetpt_eff_den2");
TH1F *h2 = (TH1F*)g1->Get("${num}");

TFile *g2 =TFile::Open("plots_${sample}_${date}_l1_90_${etarange}.root");
TH1F *h3 = (TH1F*)g2->Get("recojetpt_eff_den2");
TH1F *h4 = (TH1F*)g2->Get("${num}");

TFile *g3 =TFile::Open("plots_${sample}_${date}_l1_120_${etarange}.root");
TH1F *h5 = (TH1F*)g3->Get("recojetpt_eff_den2");
TH1F *h6 = (TH1F*)g3->Get("${num}");

TFile *g4 =TFile::Open("plots_${sample}_${date}_l1_180_${etarange}.root");
TH1F *h7 = (TH1F*)g4->Get("recojetpt_eff_den2");
TH1F *h8 = (TH1F*)g4->Get("${num}");

TEfficiency* pEff1 = 0;
TEfficiency* pEff2 = 0;
TEfficiency* pEff3 = 0;
TEfficiency* pEff4 = 0;

h1->SetLineColorAlpha(kWhite, 1.);
h1->SetTitle("Sub-leading jet");
h1->GetXaxis()->SetTitle("Offline p_{T} [GeV]");
h1->GetXaxis()->SetTitleSize(0.04);
h1->GetXaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitle("Efficiency");
h1->GetYaxis()->SetTitleSize(0.04);
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetRangeUser(0.,1.5);
h1->Draw();

c1->Update();
c1->RedrawAxis();
TLine l;
l.DrawLine(c1->GetUxmin(), c1->GetUymax(), c1->GetUxmax(), c1->GetUymax());
l.DrawLine(c1->GetUxmax(), c1->GetUymin(), c1->GetUxmax(), c1->GetUymax());

if(TEfficiency::CheckConsistency(*h2,*h1))
{ 
  pEff1 = new TEfficiency(*h2,*h1);
  //pEff1->SetTitle("Leading jet;Offline p_{T} [GeV];Efficiency");
  pEff1->SetLineWidth(2.);
  pEff1->SetLineColor(kBlue-4);
  pEff1->Draw("same");
}

if(TEfficiency::CheckConsistency(*h4,*h3))
{
  pEff2 = new TEfficiency(*h4,*h3);
  pEff2->SetLineWidth(2.);
  pEff2->SetLineColor(kCyan);
  pEff2->Draw("same");
}

if(TEfficiency::CheckConsistency(*h6,*h5))
{
  pEff3 = new TEfficiency(*h6,*h5);
  pEff3->SetLineWidth(2.);
  pEff3->SetLineColor(kSpring-9);
  pEff3->Draw("same");
}

if(TEfficiency::CheckConsistency(*h8,*h7))
{
  pEff4 = new TEfficiency(*h8,*h7);
  pEff4->SetLineWidth(2.);
  pEff4->SetLineColor(kPink-3);
  pEff4->Draw("same");
}

TLegend *legend1 = new TLegend(0.15, 0.65, 0.35, 0.85);
legend1->SetTextFont(42);
legend1->SetLineColor(0);
legend1->SetTextSize(0.04);
legend1->SetFillColor(0);
legend1->AddEntry(pEff1, "${labelprint} jet p_{T} > 35 GeV", "l");
legend1->AddEntry(pEff2, "${labelprint} jet p_{T} > 90 GeV", "l");
legend1->AddEntry(pEff3, "${labelprint} jet p_{T} > 120 GeV", "l");
legend1->AddEntry(pEff4, "${labelprint} jet p_{T} > 180 GeV", "l");
legend1->Draw();

TLatex *t2b = new TLatex(0.75,0.75,"${rangeprint}");
t2b->SetNDC();
t2b->SetTextFont(42);
t2b->SetTextSize(0.04);
t2b->SetTextAlign(20);
t2b->Draw("same");
c1->SaveAs("${date}/efficiencies_2_${etarange}_${sample}_${labelprint}.png");

}

EOF

root -l -q dummy.C
rm -rf dummy.C

done

done

exit 0
