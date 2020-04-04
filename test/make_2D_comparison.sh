# Usage: ./make_2D_comparison.sh sample date
sample=$1
date=$2
echo "Making 2D comparison plots..."
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

for name in pt_comparison1 pt_comparison2 pt_comparison1_s2 pt_comparison2_s2 eta_comparison1 eta_comparison2 eta_comparison1_s2 eta_comparison2_s2 phi_comparison1 phi_comparison2 phi_comparison1_s2 phi_comparison2_s2; do

cat>dummy.C<<EOF
{
TCanvas *c1 =new TCanvas("c1", " ", 0, 0,1000,700);
c1->Range(0,0,1,1);
c1->SetFillColor(0);
c1->SetBorderMode(0);
c1->SetBorderSize(2);
c1->SetFrameBorderMode(0);
c1->Draw();
gStyle->SetPalette(kLightTemperature);
gStyle->SetOptStat(0);

TFile *g1 =TFile::Open("plots_${sample}_${date}_l1_20_${etarange}.root");

${name}->Draw("colz");
${name}->GetXaxis()->SetTitleSize(0.045);
${name}->GetXaxis()->SetTitleOffset(1.2);
${name}->GetYaxis()->SetTitleSize(0.045);
${name}->GetYaxis()->SetTitleOffset(1.2);

TLatex *t2b = new TLatex(0.75,0.25,"${rangeprint}");
t2b->SetNDC();
t2b->SetTextFont(42);
t2b->SetTextSize(0.04);
t2b->SetTextAlign(20);
t2b->Draw("same");

c1->SaveAs("${date}/${name}_${etarange}_${sample}.png");
}

EOF

root -l -q dummy.C
rm -rf dummy.C

done

done

exit 0
