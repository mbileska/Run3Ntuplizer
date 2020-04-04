#define controlplot_cxx
#include "controlplot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

int main(int argc, char *argv[])
{

  if(argc > 1)
    {
      controlplot t(argv[1], argv[2], argv[3], argv[4]);
      t.Loop(argv[3], argv[4]);
    }
  return 0;
}

using namespace std;

void controlplot::Loop(const char* recoeta, const char* l1pt)
{
//   In a ROOT session, you can do:
//      root> .L controlplot.C
//      root> controlplot t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   std::stringstream sstrm1(recoeta);
   std::string recoetarange;
   sstrm1 >> recoetarange;
   bool recoetacut;
   std::stringstream sstrm2(l1pt);
   double l1ptcut;
   sstrm2 >> l1ptcut;
   double SF = 1.2;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;
	//cout<<l1Pt_1<<"\t"<<l1Pt_2<<endl;
	l1Pt_1 = SF*l1Pt_1; l1Pt_2 = SF*l1Pt_2; // naive scale factor
	if(l1Pt_1 >= 20) { l1jetpt_1->Fill(l1Pt_1); l1jeteta_1->Fill(l1Eta_1); l1jetphi_1->Fill(l1Phi_1); }
	if(l1Pt_2 >= 20) { l1jetpt_2->Fill(l1Pt_2); l1jeteta_2->Fill(l1Eta_2); l1jetphi_2->Fill(l1Phi_2); }
        if(stage2Pt_1 >= 20) { s2jetpt_1->Fill(stage2Pt_1); s2jeteta_1->Fill(stage2Eta_1); s2jetphi_1->Fill(stage2Phi_1); }
        if(stage2Pt_2 >= 20) { s2jetpt_2->Fill(stage2Pt_2); s2jeteta_2->Fill(stage2Eta_2); s2jetphi_2->Fill(stage2Phi_2); }
	recoetacut = false;
	if(recoetarange == "barrel") recoetacut = (abs(recoEta_1) <= 1.474);
	if(recoetarange == "endcap") recoetacut = (abs(recoEta_1) > 1.474 && abs(recoEta_1) <= 3);
	if(recoetarange == "forward") recoetacut = (abs(recoEta_1) > 3 && abs(recoEta_1) <= 5);
	if(recoetarange == "all") recoetacut = true;
	if(recoPt_1 >= 20 && recoetacut) { 
		recojetpt_eff_den1->Fill(recoPt_1); 
		recojeteta_eff_den1->Fill(recoEta_1); 
		recojetphi_eff_den1->Fill(recoPhi_1); 
		if (l1Pt_1 >= l1ptcut) { 
			recojetpt_eff_num1->Fill(recoPt_1);  
			recojeteta_eff_num1->Fill(recoEta_1); 
			recojetphi_eff_num1->Fill(recoPhi_1); 
			l1jetpt_rate1->Fill(l1Pt_1); 
			reso_pt1->Fill((l1Pt_1-recoPt_1)/recoPt_1);
			reso_eta1->Fill((l1Eta_1-recoEta_1)/recoEta_1);
			reso_phi1->Fill((l1Phi_1-recoPhi_1)/recoPhi_1);
			pt_comparison1->Fill(recoPt_1, l1Pt_1);
			eta_comparison1->Fill(recoEta_1, l1Eta_1);
			phi_comparison1->Fill(recoPhi_1, l1Phi_1);
		} 
                if (stage2Pt_1 >= l1ptcut) {
                        recojetpt_eff_num1_s2->Fill(recoPt_1);
                        recojeteta_eff_num1_s2->Fill(recoEta_1);
                        recojetphi_eff_num1_s2->Fill(recoPhi_1);
                        s2jetpt_rate1->Fill(stage2Pt_1);
                        reso_pt1_s2->Fill((stage2Pt_1-recoPt_1)/recoPt_1);
                        reso_eta1_s2->Fill((stage2Eta_1-recoEta_1)/recoEta_1);
                        reso_phi1_s2->Fill((stage2Phi_1-recoPhi_1)/recoPhi_1);
			pt_comparison1_s2->Fill(recoPt_1, stage2Pt_1);
			eta_comparison1_s2->Fill(recoEta_1, stage2Eta_1);
			phi_comparison1_s2->Fill(recoPhi_1, stage2Phi_1);
                }
	}
	if(recoPt_2 >= 20 && recoetacut) { 
		recojetpt_eff_den2->Fill(recoPt_2); 
		recojeteta_eff_den2->Fill(recoEta_2);
		recojetphi_eff_den2->Fill(recoPhi_2);
		if (l1Pt_2 >= l1ptcut) { 
			recojetpt_eff_num2->Fill(recoPt_2);  
			recojeteta_eff_num2->Fill(recoEta_2); 
			recojetphi_eff_num2->Fill(recoPhi_2); 
			l1jetpt_rate2->Fill(l1Pt_2); 
			reso_pt2->Fill((l1Pt_2-recoPt_2)/recoPt_2);
			reso_eta2->Fill((l1Eta_2-recoEta_2)/recoEta_2);
			reso_phi2->Fill((l1Phi_2-recoPhi_2)/recoPhi_2);
			pt_comparison2->Fill(recoPt_2, l1Pt_2);
			eta_comparison2->Fill(recoEta_2, l1Eta_2);
			phi_comparison2->Fill(recoPhi_2, l1Phi_2);
		}
                if (stage2Pt_2 >= l1ptcut) {
                        recojetpt_eff_num2_s2->Fill(recoPt_2);
                        recojeteta_eff_num2_s2->Fill(recoEta_2);
                        recojetphi_eff_num2_s2->Fill(recoPhi_2);
                        s2jetpt_rate2->Fill(stage2Pt_2);
                        reso_pt2_s2->Fill((stage2Pt_2-recoPt_2)/recoPt_2);
                        reso_eta2_s2->Fill((stage2Eta_2-recoEta_2)/recoEta_2);
                        reso_phi2_s2->Fill((stage2Phi_2-recoPhi_2)/recoPhi_2);
			pt_comparison2_s2->Fill(recoPt_2, stage2Pt_2);
			eta_comparison2_s2->Fill(recoEta_2, stage2Eta_2);
			phi_comparison2_s2->Fill(recoPhi_2, stage2Phi_2);
                }

	}
      // if (Cut(ientry) < 0) continue;
   }
}

void controlplot::BookHistos(const char* file2){
        fileName = new TFile(file2, "RECREATE");
        fileName->cd();
        char name[100];
	//Float_t binb[] = { -5.0925, -4.8575, -4.650, -4.470, -4.290, -4.115, -3.940, -3.765, -3.590, -3.415, -3.240, -3.065, -2.825, -2.575, -2.411, -2.247, -2.1075, -1.9865, -1.880, -1.785, -1.7004, -1.6132, -1.526, -1.4388, -1.3516, -1.2644, -1.1772, -1.09, -1.0028, -0.9156, -0.8284, -0.7412, -0.654, -0.5668, -0.4796, -0.3924, -0.3052, -0.218, -0.1308, -0.0436, 0, 0.0436, 0.1308, 0.218, 0.3052, 0.3924, 0.4796, 0.5668, 0.654, 0.7412, 0.8284, 0.9156, 1.0028, 1.09, 1.1772, 1.2644, 1.3516, 1.4388, 1.526, 1.6132, 1.7004, 1.785, 1.880, 1.9865, 2.1075, 2.247, 2.411, 2.575, 2.825, 3.065, 3.240, 3.415, 3.590, 3.765, 3.940, 4.115, 4.290, 4.470, 4.650, 4.8575, 5.0925 };
	Float_t binb[] = { -5., -4.5, -4., -3.5, -3., -2.5, -2.25, -2., -1.75, -1.5, -1.25, -1., -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 3., 3.5, 4., 4.5, 5.};
	Int_t  binnum = 30;
	//l1jeteta_1 = new TH1F (name,"l1jeteta_1", binnum, binb);

        sprintf(name, "l1jetpt_1");
        l1jetpt_1 = new TH1F (name,"l1jetpt_1", 40, 20, 250);
	l1jetpt_1->SetTitle("Leading L1 jet");
        l1jetpt_1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "l1jeteta_1");
        //l1jeteta_1 = new TH1F (name,"l1jeteta_1", 40, -5.,5.);
        l1jeteta_1 = new TH1F (name,"l1jeteta_1", binnum, binb);
        l1jeteta_1->SetTitle("Leading L1 jet");
        l1jeteta_1->GetXaxis()->SetTitle("#eta");

        sprintf(name, "l1jetphi_1");
        l1jetphi_1 = new TH1F (name,"l1jetphi_1", 20, -M_PI, M_PI);
        l1jetphi_1->SetTitle("Leading L1 jet");
        l1jetphi_1->GetXaxis()->SetTitle("#phi");

        sprintf(name, "l1jetpt_2");
        l1jetpt_2 = new TH1F (name,"l1jetpt_2", 40, 20, 250);
        l1jetpt_2->SetTitle("Sub-leading L1 jet");
        l1jetpt_2->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "l1jeteta_2");
        //l1jeteta_2 = new TH1F (name,"l1jeteta_2", 40, -5.,5.);
        l1jeteta_2 = new TH1F (name,"l1jeteta_2", binnum, binb);
        l1jeteta_2->SetTitle("Sub-leading L1 jet");
        l1jeteta_2->GetXaxis()->SetTitle("#eta");

        sprintf(name, "l1jetphi_2");
        l1jetphi_2 = new TH1F (name,"l1jetphi_2", 20, -M_PI, M_PI);
        l1jetphi_2->SetTitle("Sub-leading L1 jet");
        l1jetphi_2->GetXaxis()->SetTitle("#phi");

        sprintf(name, "l1jetpt_rate1");
        l1jetpt_rate1 = new TH1F (name,"l1jetpt_rate1", 40, 20, 250);
        l1jetpt_rate1->SetTitle("Leading L1 jet");
        l1jetpt_rate1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "l1jetpt_rate2");
        l1jetpt_rate2 = new TH1F (name,"l1jetpt_rate2", 40, 20, 250);
        l1jetpt_rate2->SetTitle("Sub-leading L1 jet");
        l1jetpt_rate2->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "s2jetpt_1");
        s2jetpt_1 = new TH1F (name,"s2jetpt_1", 40, 20, 250);
        s2jetpt_1->SetTitle("Leading Stage2 jet");
        s2jetpt_1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "s2jeteta_1");
        //s2jeteta_1 = new TH1F (name,"s2jeteta_1", 40, -5.,5.);
        s2jeteta_1 = new TH1F (name,"s2jeteta_1", binnum, binb);
        s2jeteta_1->SetTitle("Leading Stage2 jet");
        s2jeteta_1->GetXaxis()->SetTitle("#eta");

        sprintf(name, "s2jetphi_1");
        s2jetphi_1 = new TH1F (name,"s2jetphi_1", 20, -M_PI, M_PI);
        s2jetphi_1->SetTitle("Leading Stage2 jet");
        s2jetphi_1->GetXaxis()->SetTitle("#phi");

        sprintf(name, "s2jetpt_2");
        s2jetpt_2 = new TH1F (name,"s2jetpt_2", 40, 20, 250);
        s2jetpt_2->SetTitle("Sub-leading Stage2 jet");
        s2jetpt_2->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "s2jeteta_2");
        //s2jeteta_2 = new TH1F (name,"s2jeteta_2", 40, -5.,5.);
        s2jeteta_2 = new TH1F (name,"s2jeteta_2", binnum, binb);
        s2jeteta_2->SetTitle("Sub-leading Stage2 jet");
        s2jeteta_2->GetXaxis()->SetTitle("#eta");

        sprintf(name, "s2jetphi_2");
        s2jetphi_2 = new TH1F (name,"s2jetphi_2", 20, -M_PI, M_PI);
        s2jetphi_2->SetTitle("Sub-leading Stage2 jet");
        s2jetphi_2->GetXaxis()->SetTitle("#phi");

        sprintf(name, "s2jetpt_rate1");
        s2jetpt_rate1 = new TH1F (name,"s2jetpt_rate1", 40, 20, 250);
        s2jetpt_rate1->SetTitle("Leading Stage2 jet");
        s2jetpt_rate1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "s2jetpt_rate2");
        s2jetpt_rate2 = new TH1F (name,"s2jetpt_rate2", 40, 20, 250);
        s2jetpt_rate2->SetTitle("Sub-leading Stage2 jet");
        s2jetpt_rate2->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojetpt_eff_den1");
        recojetpt_eff_den1 = new TH1F (name,"recojetpt_eff_den1", 40, 20, 250);
        recojetpt_eff_den1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojetpt_eff_den2");
        recojetpt_eff_den2 = new TH1F (name,"recojetpt_eff_den2", 40, 20, 250);
        recojetpt_eff_den2->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojetpt_eff_num1");
        recojetpt_eff_num1 = new TH1F (name,"recojetpt_eff_num1", 40, 20, 250);
        recojetpt_eff_num1->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojetpt_eff_num2");
        recojetpt_eff_num2 = new TH1F (name,"recojetpt_eff_num2", 40, 20, 250);
        recojetpt_eff_num2->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojetpt_eff_num1_s2");
        recojetpt_eff_num1_s2 = new TH1F (name,"recojetpt_eff_num1_s2", 40, 20, 250);
        recojetpt_eff_num1_s2->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojetpt_eff_num2_s2");
        recojetpt_eff_num2_s2 = new TH1F (name,"recojetpt_eff_num2_s2", 40, 20, 250);
        recojetpt_eff_num2_s2->GetXaxis()->SetTitle("p_{T} [GeV]");

        sprintf(name, "recojeteta_eff_den1");
        recojeteta_eff_den1 = new TH1F (name,"recojeteta_eff_den1", 40, -5, 5);
        recojeteta_eff_den1->GetXaxis()->SetTitle("#eta");

        sprintf(name, "recojeteta_eff_den2");
        recojeteta_eff_den2 = new TH1F (name,"recojeteta_eff_den2", 40, -5, 5);
        recojeteta_eff_den2->GetXaxis()->SetTitle("#eta");

        sprintf(name, "recojeteta_eff_num1");
        recojeteta_eff_num1 = new TH1F (name,"recojeteta_eff_num1", 40, -5, 5);
        recojeteta_eff_num1->GetXaxis()->SetTitle("#eta");

        sprintf(name, "recojeteta_eff_num2");
        recojeteta_eff_num2 = new TH1F (name,"recojetpt_eff_num2", 40, -5, 5);
        recojeteta_eff_num2->GetXaxis()->SetTitle("#eta");

        sprintf(name, "recojeteta_eff_num1_s2");
        recojeteta_eff_num1_s2 = new TH1F (name,"recojeteta_eff_num1_s2", 40, -5, 5);
        recojeteta_eff_num1_s2->GetXaxis()->SetTitle("#eta");

        sprintf(name, "recojeteta_eff_num2_s2");
        recojeteta_eff_num2_s2 = new TH1F (name,"recojetpt_eff_num2_s2", 40, -5, 5);
        recojeteta_eff_num2_s2->GetXaxis()->SetTitle("#eta");

        sprintf(name, "recojetphi_eff_den1");
        recojetphi_eff_den1 = new TH1F (name,"recojetphi_eff_den1", 40, -M_PI, M_PI);
        recojetphi_eff_den1->GetXaxis()->SetTitle("#phi");

        sprintf(name, "recojetphi_eff_den2");
        recojetphi_eff_den2 = new TH1F (name,"recojetphi_eff_den2", 40, -M_PI, M_PI);
        recojetphi_eff_den2->GetXaxis()->SetTitle("#phi");

        sprintf(name, "recojetphi_eff_num1");
        recojetphi_eff_num1 = new TH1F (name,"recojetphi_eff_num1", 40, -M_PI, M_PI);
        recojetphi_eff_num1->GetXaxis()->SetTitle("#phi");

        sprintf(name, "recojetphi_eff_num2");
        recojetphi_eff_num2 = new TH1F (name,"recojetphi_eff_num2", 40, -M_PI, M_PI);
        recojetphi_eff_num2->GetXaxis()->SetTitle("#phi");

        sprintf(name, "recojetphi_eff_num1_s2");
        recojetphi_eff_num1_s2 = new TH1F (name,"recojetphi_eff_num1_s2", 40, -M_PI, M_PI);
        recojetphi_eff_num1_s2->GetXaxis()->SetTitle("#phi");

        sprintf(name, "recojetphi_eff_num2_s2");
        recojetphi_eff_num2_s2 = new TH1F (name,"recojetphi_eff_num2_s2", 40, -M_PI, M_PI);
        recojetphi_eff_num2_s2->GetXaxis()->SetTitle("#phi");

        sprintf(name, "reso_pt1");
        reso_pt1 = new TH1F (name,"reso_pt1", 40, -2, 2);
        reso_pt1->SetTitle("Leading L1 jet");
        reso_pt1->GetXaxis()->SetTitle("relative p_{T} difference");

        sprintf(name, "reso_eta1");
        reso_eta1 = new TH1F (name,"reso_eta1", 40, -0.3, 0.3);
        reso_eta1->SetTitle("Leading L1 jet");
        reso_eta1->GetXaxis()->SetTitle("relative #eta difference");

        sprintf(name, "reso_phi1");
        reso_phi1 = new TH1F (name,"reso_phi1", 40, -0.3, 0.3);
        reso_phi1->SetTitle("Leading L1 jet");
        reso_phi1->GetXaxis()->SetTitle("relative #phi difference");

        sprintf(name, "reso_pt2");
        reso_pt2 = new TH1F (name,"reso_pt2", 40, -2, 2);
        reso_pt2->SetTitle("Sub-leading L1 jet");
        reso_pt2->GetXaxis()->SetTitle("relative p_{T} difference");

        sprintf(name, "reso_eta2");
        reso_eta2 = new TH1F (name,"reso_eta2", 40, -0.3, 0.3);
        reso_eta2->SetTitle("Sub-leading L1 jet");
        reso_eta2->GetXaxis()->SetTitle("relative #eta difference");

        sprintf(name, "reso_phi2");
        reso_phi2 = new TH1F (name,"reso_phi2", 40, -0.3, 0.3);
        reso_phi2->SetTitle("Sub-leading L1 jet");
        reso_phi2->GetXaxis()->SetTitle("relative #phi difference");

	sprintf(name, "pt_comparison1");
        pt_comparison1 = new TH2F (name,"pt_comparison1", 50, 0, 500, 50, 0, 500);
        pt_comparison1->SetTitle("Leading jets");
        pt_comparison1->GetYaxis()->SetTitle("Stage3 p_{T} [GeV]");
        pt_comparison1->GetXaxis()->SetTitle("Offline p_{T} [GeV]");

	sprintf(name, "pt_comparison2");
        pt_comparison2 = new TH2F (name,"pt_comparison2", 50, 0, 500, 50, 0, 500);
	pt_comparison2->SetTitle("Sub-leading jets");
        pt_comparison2->GetYaxis()->SetTitle("Stage3 p_{T} [GeV]");
        pt_comparison2->GetXaxis()->SetTitle("Offline p_{T} [GeV]");

        sprintf(name, "eta_comparison1");
        eta_comparison1 = new TH2F (name,"eta_comparison1", 50, -5, 5, 50, -5, 5);
        eta_comparison1->SetTitle("Leading jets");
        eta_comparison1->GetYaxis()->SetTitle("Stage3 #eta");
        eta_comparison1->GetXaxis()->SetTitle("Offline #eta");

        sprintf(name, "eta_comparison2");
        eta_comparison2 = new TH2F (name,"eta_comparison2", 50, -5, 5, 50, -5, 5);
        eta_comparison2->SetTitle("Sub-leading jets");
        eta_comparison2->GetYaxis()->SetTitle("Stage3 #eta");
        eta_comparison2->GetXaxis()->SetTitle("Offline #eta");

        sprintf(name, "phi_comparison1");
        phi_comparison1 = new TH2F (name,"phi_comparison1", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
        phi_comparison1->SetTitle("Leading jets");
        phi_comparison1->GetYaxis()->SetTitle("Stage3 #phi");
        phi_comparison1->GetXaxis()->SetTitle("Offline #phi");

        sprintf(name, "phi_comparison2");
        phi_comparison2 = new TH2F (name,"phi_comparison2", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
        phi_comparison2->SetTitle("Sub-leading jets");
        phi_comparison2->GetYaxis()->SetTitle("Stage3 #phi");
        phi_comparison2->GetXaxis()->SetTitle("Offline #phi");

        sprintf(name, "reso_pt1_s2");
        reso_pt1_s2 = new TH1F (name,"reso_pt1_s2", 40, -2, 2);
        reso_pt1_s2->SetTitle("Leading Stage2 jet");
        reso_pt1_s2->GetXaxis()->SetTitle("relative p_{T} difference");

        sprintf(name, "reso_eta1_s2");
        reso_eta1_s2 = new TH1F (name,"reso_eta1_s2", 40, -0.3, 0.3);
        reso_eta1_s2->SetTitle("Leading Stage2 jet");
        reso_eta1_s2->GetXaxis()->SetTitle("relative #eta difference");

        sprintf(name, "reso_phi1_s2");
        reso_phi1_s2 = new TH1F (name,"reso_phi1_s2", 40, -0.3, 0.3);
        reso_phi1_s2->SetTitle("Leading Stage2 jet");
        reso_phi1_s2->GetXaxis()->SetTitle("relative #phi difference");

        sprintf(name, "reso_pt2_s2");
        reso_pt2_s2 = new TH1F (name,"reso_pt2_s2", 40, -2, 2);
        reso_pt2_s2->SetTitle("Sub-leading Stage2 jet");
        reso_pt2_s2->GetXaxis()->SetTitle("relative p_{T} difference");

        sprintf(name, "reso_eta2_s2");
        reso_eta2_s2 = new TH1F (name,"reso_eta2_s2", 40, -0.3, 0.3);
        reso_eta2_s2->SetTitle("Sub-leading Stage2 jet");
        reso_eta2_s2->GetXaxis()->SetTitle("relative #eta difference");

        sprintf(name, "reso_phi2_s2");
        reso_phi2_s2 = new TH1F (name,"reso_phi2_s2", 40, -0.3, 0.3);
        reso_phi2_s2->SetTitle("Sub-leading Stage2 jet");
        reso_phi2_s2->GetXaxis()->SetTitle("relative #phi difference");

        sprintf(name, "pt_comparison1_s2");
        pt_comparison1_s2 = new TH2F (name,"pt_comparison1_s2", 50, 0, 500, 50, 0, 500);
        pt_comparison1_s2->SetTitle("Leading jets");
        pt_comparison1_s2->GetYaxis()->SetTitle("Stage2 p_{T} [GeV]");
        pt_comparison1_s2->GetXaxis()->SetTitle("Offline p_{T} [GeV]");

        sprintf(name, "pt_comparison2_s2");
        pt_comparison2_s2 = new TH2F (name,"pt_comparison2_s2", 50, 0, 500, 50, 0, 500);
        pt_comparison2_s2->SetTitle("Sub-leading jets");
        pt_comparison2_s2->GetYaxis()->SetTitle("Stage2 p_{T} [GeV]");
        pt_comparison2_s2->GetXaxis()->SetTitle("Offline p_{T} [GeV]");

        sprintf(name, "eta_comparison1_s2");
        eta_comparison1_s2 = new TH2F (name,"eta_comparison1_s2", 50, -5, 5, 50, -5, 5);
        eta_comparison1_s2->SetTitle("Leading jets");
        eta_comparison1_s2->GetYaxis()->SetTitle("Stage2 #eta");
        eta_comparison1_s2->GetXaxis()->SetTitle("Offline #eta");

        sprintf(name, "eta_comparison2_s2");
        eta_comparison2_s2 = new TH2F (name,"eta_comparison2_s2", 50, -5, 5, 50, -5, 5);
        eta_comparison2_s2->SetTitle("Sub-leading jets");
        eta_comparison2_s2->GetYaxis()->SetTitle("Stage2 #eta");
        eta_comparison2_s2->GetXaxis()->SetTitle("Offline #eta");

        sprintf(name, "phi_comparison1_s2");
        phi_comparison1_s2 = new TH2F (name,"phi_comparison1_s2", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
        phi_comparison1_s2->SetTitle("Leading jets");
        phi_comparison1_s2->GetYaxis()->SetTitle("Stage2 #phi");
        phi_comparison1_s2->GetXaxis()->SetTitle("Offline #phi");

        sprintf(name, "phi_comparison2_s2");
        phi_comparison2_s2 = new TH2F (name,"phi_comparison2_s2", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
        phi_comparison2_s2->SetTitle("Sub-leading jets");
        phi_comparison2_s2->GetYaxis()->SetTitle("Stage2 #phi");
        phi_comparison2_s2->GetXaxis()->SetTitle("Offline #phi");
}
