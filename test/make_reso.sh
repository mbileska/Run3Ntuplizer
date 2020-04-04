# Usage: ./make_reso.sh date
date=$1
echo "Making resolution plots..."
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

for type in pt eta phi; do

for index in 1 2; do

if [ ${index} = '1' ]; then
	labelprint="Leading jet"
elif [ ${index} = '2' ]; then
	labelprint="Sub-leading jet"
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

TFile *g1 =TFile::Open("plots_ZeroBias_${date}_l1_20_${etarange}.root");
TH1F *h1 = (TH1F*)g1->Get("reso_${type}${index}");
h1->Scale(1./h1->Integral());
h1->SetLineWidth(2.);
h1->SetLineColor(kPink-9);
h1->GetXaxis()->SetTitleSize(0.045);
h1->GetXaxis()->SetTitleOffset(1.0);
h1->GetYaxis()->SetTitle("a. u.");
h1->GetYaxis()->SetTitleOffset(1.1);
h1->GetYaxis()->SetTitleSize(0.045);
h1->SetTitle("${labelprint}");
double a = h1->GetMaximum();

TH1F *h2 = (TH1F*)g1->Get("reso_${type}${index}_s2");
h2->Scale(1./h2->Integral());
h2->SetLineWidth(2.);
h2->SetLineColor(kAzure);
double b = h2->GetMaximum();

TFile *g2 =TFile::Open("plots_VBF_${date}_l1_20_${etarange}.root");
TH1F *h3 = (TH1F*)g2->Get("reso_${type}${index}");
h3->Scale(1./h3->Integral());
h3->SetLineWidth(2.);
h3->SetLineStyle(2);
h3->SetLineColor(kPink-9);
double c = h3->GetMaximum();

TH1F *h4 = (TH1F*)g2->Get("reso_${type}${index}_s2");
h4->Scale(1./h4->Integral());
h4->SetLineWidth(2.);
h4->SetLineStyle(2);
h4->SetLineColor(kAzure);
double d = h4->GetMaximum();

double maximum = std::max(a, std::max(b, std::max(c, d)));
h1->SetMaximum(1.2*maximum);
h1->Draw();
h2->Draw("same");
h3->Draw("same");
h4->Draw("same");

TLegend *legend1 = new TLegend(0.15, 0.65, 0.35, 0.85);
//TLegend *legend1 = new TLegend(0.15, 0.75, 0.4, 0.85);
legend1->SetTextFont(42);
legend1->SetLineColor(0);
legend1->SetTextSize(0.04);
legend1->SetFillColor(0);
legend1->AddEntry(h1, "Stage-3 (Data)", "l");
legend1->AddEntry(h2, "Stage-2 (Data)", "l");
legend1->AddEntry(h3, "Stage-3 (MC)", "l");
legend1->AddEntry(h4, "Stage-2 (MC)", "l");

legend1->Draw();

TLatex *t2b = new TLatex(0.75,0.75,"${rangeprint}");
t2b->SetNDC();
t2b->SetTextFont(42);
t2b->SetTextSize(0.04);
t2b->SetTextAlign(20);
t2b->Draw("same");

gStyle->SetOptStat(0);

c1->SaveAs("${date}/reso_${index}_${type}_${etarange}.png");

}

EOF

root -l -q dummy.C
rm -rf dummy.C

done

done

done

exit 0
