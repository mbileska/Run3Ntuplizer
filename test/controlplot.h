//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 18 15:14:21 2020 by ROOT version 6.14/09
// from TTree efficiencyTree/Reco Matched Jet Tree 
// found on file: l1TNtuple-ZeroBias_test.root
//////////////////////////////////////////////////////////

#ifndef controlplot_h
#define controlplot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMap.h>
#include <vector>
//#ifdef __MAKECINT__
//#pragma link C++ class vector<std::vector<int> >+;
//#pragma link C++ class vector<std::vector<std::string> >+;
//#pragma link C++ class vector<std::vector<float> >+;
//#pragma link C++ class vector<std::vector<bool> >+;
//#pragma extra_include "vector";
//#endif

using namespace std;

// Header file for the classes stored in the TTree if any.

class controlplot {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
 
   TFile *fileName;
   TH1F *l1jetpt_1, *l1jetpt_2, *l1jeteta_1, *l1jeteta_2, *l1jetphi_1, *l1jetphi_2, *recojetpt_eff_den1, *recojetpt_eff_den2, *recojetpt_eff_num1, *recojetpt_eff_num2, *reso_pt1, *reso_eta1, *reso_phi1, *reso_pt2, *reso_eta2, *reso_phi2;
   TH1F *l1jetpt_rate1, *l1jetpt_rate2, *recojeteta_eff_den1, *recojeteta_eff_den2, *recojeteta_eff_num1, *recojeteta_eff_num2, *recojetphi_eff_den1, *recojetphi_eff_den2, *recojetphi_eff_num1, *recojetphi_eff_num2;
   TH1F *s2jetpt_1, *s2jetpt_2, *s2jeteta_1, *s2jeteta_2, *s2jetphi_1, *s2jetphi_2, *recojetpt_eff_num1_s2, *recojetpt_eff_num2_s2, *reso_pt1_s2, *reso_eta1_s2, *reso_phi1_s2, *reso_pt2_s2, *reso_eta2_s2, *reso_phi2_s2, *recojeteta_eff_num1_s2, *recojeteta_eff_num2_s2, *recojetphi_eff_num1_s2, *recojetphi_eff_num2_s2, *s2jetpt_rate1, *s2jetpt_rate2;
   TH2F *pt_comparison1, *pt_comparison2, *eta_comparison1, *eta_comparison2, *phi_comparison1, *phi_comparison2, *pt_comparison1_s2, *pt_comparison2_s2, *eta_comparison1_s2, *eta_comparison2_s2, *phi_comparison1_s2, *phi_comparison2_s2;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Double_t        recoPt_1;
   Double_t        recoEta_1;
   Double_t        recoPhi_1;
   Int_t           recoNthJet_1;
   Double_t        recoPt_2;
   Double_t        recoEta_2;
   Double_t        recoPhi_2;
   Int_t           recoNthJet_2;
   Double_t        recoDeltaEta;
   Double_t        recoDeltaPhi;
   Double_t        recoDeltaR;
   Double_t        recoMass;
   Double_t        l1Pt_1;
   Double_t        l1Eta_1;
   Double_t        l1Phi_1;
   Int_t           l1NthJet_1;
   Double_t        l1Pt_2;
   Double_t        l1Eta_2;
   Double_t        l1Phi_2;
   Int_t           l1NthJet_2;
   Double_t        l1DeltaEta;
   Double_t        l1DeltaPhi;
   Double_t        l1DeltaR;
   Double_t        l1Mass;
   Int_t           l1Matched_1;
   Int_t           l1Matched_2;
   Int_t           nRecoJets;
   Int_t           nL1Jets;
   Double_t        stage2Pt_1;
   Double_t        stage2Eta_1;
   Double_t        stage2Phi_1;
   Int_t           stage2NthJet_1;
   Double_t        stage2Pt_2;
   Double_t        stage2Eta_2;
   Double_t        stage2Phi_2;
   Int_t           stage2NthJet_2;
   Double_t        stage2DeltaEta;
   Double_t        stage2DeltaPhi;
   Double_t        stage2DeltaR;
   Double_t        stage2Mass;
   Int_t           stage2Matched_1;
   Int_t           stage2Matched_2;
   Int_t           nStage2Jets;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_recoPt_1;   //!
   TBranch        *b_recoEta_1;   //!
   TBranch        *b_recoPhi_1;   //!
   TBranch        *b_recoNthJet_1;   //!
   TBranch        *b_recoPt_2;   //!
   TBranch        *b_recoEta_2;   //!
   TBranch        *b_recoPhi_2;   //!
   TBranch        *b_recoNthJet_2;   //!
   TBranch        *b_recoDeltaEta;   //!
   TBranch        *b_recoDeltaPhi;   //!
   TBranch        *b_recoDeltaR;   //!
   TBranch        *b_recoMass;   //!
   TBranch        *b_l1Pt_1;   //!
   TBranch        *b_l1Eta_1;   //!
   TBranch        *b_l1Phi_1;   //!
   TBranch        *b_l1NthJet_1;   //!
   TBranch        *b_l1Pt_2;   //!
   TBranch        *b_l1Eta_2;   //!
   TBranch        *b_l1Phi_2;   //!
   TBranch        *b_l1NthJet_2;   //!
   TBranch        *b_l1DeltaEta;   //!
   TBranch        *b_l1DeltaPhi;   //!
   TBranch        *b_l1DeltaR;   //!
   TBranch        *b_l1Mass;   //!
   TBranch        *b_l1Matched_1;   //!
   TBranch        *b_l1Matched_2;   //!
   TBranch        *b_nRecoJets;   //!
   TBranch        *b_nL1Jets;   //!
   TBranch        *b_stage2Pt_1;   //!
   TBranch        *b_stage2Eta_1;   //!
   TBranch        *b_stage2Phi_1;   //!
   TBranch        *b_stage2NthJet_1;   //!
   TBranch        *b_stage2Pt_2;   //!
   TBranch        *b_stage2Eta_2;   //!
   TBranch        *b_stage2Phi_2;   //!
   TBranch        *b_stage2NthJet_2;   //!
   TBranch        *b_stage2DeltaEta;   //!
   TBranch        *b_stage2DeltaPhi;   //!
   TBranch        *b_stage2DeltaR;   //!
   TBranch        *b_stage2Mass;   //!
   TBranch        *b_stage2Matched_1;   //!
   TBranch        *b_stage2Matched_2;   //!
   TBranch        *b_nStage2Jets;   //!

   controlplot(const char* file1, const char* file2);
   virtual ~controlplot();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     BookHistos(const char* file2);
};

#endif

#ifdef controlplot_cxx
controlplot::controlplot(const char* file1, const char* file2)
{

   BookHistos(file2);
   TChain *chain = new TChain("l1NtupleProducer/Stage3Regions/efficiencyTree");
   chain->Add(file1);
   Init(chain);
}

controlplot::~controlplot()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   fileName->Close();
}

Int_t controlplot::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t controlplot::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void controlplot::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("recoPt_1", &recoPt_1, &b_recoPt_1);
   fChain->SetBranchAddress("recoEta_1", &recoEta_1, &b_recoEta_1);
   fChain->SetBranchAddress("recoPhi_1", &recoPhi_1, &b_recoPhi_1);
   fChain->SetBranchAddress("recoNthJet_1", &recoNthJet_1, &b_recoNthJet_1);
   fChain->SetBranchAddress("recoPt_2", &recoPt_2, &b_recoPt_2);
   fChain->SetBranchAddress("recoEta_2", &recoEta_2, &b_recoEta_2);
   fChain->SetBranchAddress("recoPhi_2", &recoPhi_2, &b_recoPhi_2);
   fChain->SetBranchAddress("recoNthJet_2", &recoNthJet_2, &b_recoNthJet_2);
   fChain->SetBranchAddress("recoDeltaEta", &recoDeltaEta, &b_recoDeltaEta);
   fChain->SetBranchAddress("recoDeltaPhi", &recoDeltaPhi, &b_recoDeltaPhi);
   fChain->SetBranchAddress("recoDeltaR", &recoDeltaR, &b_recoDeltaR);
   fChain->SetBranchAddress("recoMass", &recoMass, &b_recoMass);
   fChain->SetBranchAddress("l1Pt_1", &l1Pt_1, &b_l1Pt_1);
   fChain->SetBranchAddress("l1Eta_1", &l1Eta_1, &b_l1Eta_1);
   fChain->SetBranchAddress("l1Phi_1", &l1Phi_1, &b_l1Phi_1);
   fChain->SetBranchAddress("l1NthJet_1", &l1NthJet_1, &b_l1NthJet_1);
   fChain->SetBranchAddress("l1Pt_2", &l1Pt_2, &b_l1Pt_2);
   fChain->SetBranchAddress("l1Eta_2", &l1Eta_2, &b_l1Eta_2);
   fChain->SetBranchAddress("l1Phi_2", &l1Phi_2, &b_l1Phi_2);
   fChain->SetBranchAddress("l1NthJet_2", &l1NthJet_2, &b_l1NthJet_2);
   fChain->SetBranchAddress("l1DeltaEta", &l1DeltaEta, &b_l1DeltaEta);
   fChain->SetBranchAddress("l1DeltaPhi", &l1DeltaPhi, &b_l1DeltaPhi);
   fChain->SetBranchAddress("l1DeltaR", &l1DeltaR, &b_l1DeltaR);
   fChain->SetBranchAddress("l1Mass", &l1Mass, &b_l1Mass);
   fChain->SetBranchAddress("l1Matched_1", &l1Matched_1, &b_l1Matched_1);
   fChain->SetBranchAddress("l1Matched_2", &l1Matched_2, &b_l1Matched_2);
   fChain->SetBranchAddress("nRecoJets", &nRecoJets, &b_nRecoJets);
   fChain->SetBranchAddress("nL1Jets", &nL1Jets, &b_nL1Jets);
   fChain->SetBranchAddress("stage2Pt_1", &stage2Pt_1, &b_stage2Pt_1);
   fChain->SetBranchAddress("stage2Eta_1", &stage2Eta_1, &b_stage2Eta_1);
   fChain->SetBranchAddress("stage2Phi_1", &stage2Phi_1, &b_stage2Phi_1);
   fChain->SetBranchAddress("stage2NthJet_1", &stage2NthJet_1, &b_stage2NthJet_1);
   fChain->SetBranchAddress("stage2Pt_2", &stage2Pt_2, &b_stage2Pt_2);
   fChain->SetBranchAddress("stage2Eta_2", &stage2Eta_2, &b_stage2Eta_2);
   fChain->SetBranchAddress("stage2Phi_2", &stage2Phi_2, &b_stage2Phi_2);
   fChain->SetBranchAddress("stage2NthJet_2", &stage2NthJet_2, &b_stage2NthJet_2);
   fChain->SetBranchAddress("stage2DeltaEta", &stage2DeltaEta, &b_stage2DeltaEta);
   fChain->SetBranchAddress("stage2DeltaPhi", &stage2DeltaPhi, &b_stage2DeltaPhi);
   fChain->SetBranchAddress("stage2DeltaR", &stage2DeltaR, &b_stage2DeltaR);
   fChain->SetBranchAddress("stage2Mass", &stage2Mass, &b_stage2Mass);
   fChain->SetBranchAddress("stage2Matched_1", &stage2Matched_1, &b_stage2Matched_1);
   fChain->SetBranchAddress("stage2Matched_2", &stage2Matched_2, &b_stage2Matched_2);
   fChain->SetBranchAddress("nStage2Jets", &nStage2Jets, &b_nStage2Jets);

   Notify();
}

Bool_t controlplot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void controlplot::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t controlplot::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef controlplot_cxx
