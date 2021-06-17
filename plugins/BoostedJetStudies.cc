// system include files
#include <memory>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "L1Trigger/Run3Ntuplizer/plugins/helpers.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

using namespace l1extra;
using namespace std;

bool compareByPt (l1extra::L1JetParticle i, l1extra::L1JetParticle j) { return(i.pt()>j.pt()); };

//
// class declaration
//

class BoostedJetStudies : public edm::EDAnalyzer {
public:
  explicit BoostedJetStudies(const edm::ParameterSet&);
  ~BoostedJetStudies();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void zeroOutAllVariables();

private:
  void analyze(const edm::Event& evt, const edm::EventSetup& es);      
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<vector<reco::CaloJet> > jetSrc_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrcAK8_;
  edm::EDGetTokenT<vector<pat::Jet> > jetSrcCA15_;
  edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;

  edm::EDGetTokenT<BXVector<l1t::Jet>> stage2JetToken_;
  edm::EDGetTokenT<BXVector<l1t::Tau>> stage2TauToken_;
  edm::EDGetTokenT<vector<l1extra::L1JetParticle>> l1BoostedToken_;

  TH1F* nEvents;
//  TH1F* recoJet_pt;
//  TH1F* recoJet_eta;
//  TH1F* recoJet_phi;
//
//  TH1F* recoJetAK8_pt;
//  TH1F* recoJetAK8_eta;
//  TH1F* recoJetAK8_phi;

//  TTree* l1Tree;
  int run, lumi, event;

  double genPt_1, genEta_1, genPhi_1, genM_1, genDR_1;
  double genPt_2, genEta_2, genPhi_2, genM_2, genDR_2;
  int genId_1, genMother_1;
  int genId_2, genMother_2;
  double recoPt, recoEta, recoPhi;
  double recoPt_1, recoEta_1, recoPhi_1, recoMass_1;
  double recoPt_2, recoEta_2, recoPhi_2;
  double l1Pt, l1Eta, l1Phi;
  double l1Pt_1, l1Eta_1, l1Phi_1;
  double l1Pt_2, l1Eta_2, l1Phi_2;
  double seedPt, seedEta, seedPhi;
  double seedPt_1, seedEta_1, seedPhi_1;
  double seedPt_2, seedEta_2, seedPhi_2;
  
  int l1NthJet;
  int l1NthJet_1;
  int l1NthJet_2;
  int seedNthJet;
  int seedNthJet_1;
  int seedNthJet_2;

  double recoPt_;
  std::vector<int> nSubJets, nBHadrons, HFlav;
  std::vector<std::vector<int>> subJetHFlav;
  std::vector<float> tau1, tau2, tau3;

  std::vector<TLorentzVector> *l1Jets  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *seed180  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *tauseed  = new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *ak8Jets  = new std::vector<TLorentzVector>;
  //std::vector<TLorentzVector> *ca15Jets =new std::vector<TLorentzVector>;
  std::vector<TLorentzVector> *subJets  = new std::vector<TLorentzVector>;

  void createBranches(TTree *tree);
  TTree* efficiencyTree;
  edm::Service<TFileService> tfs_;  

};

BoostedJetStudies::BoostedJetStudies(const edm::ParameterSet& iConfig) :
  jetSrc_(    consumes<vector<reco::CaloJet> >(iConfig.getParameter<edm::InputTag>("recoJets"))),
  jetSrcAK8_( consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsAK8"))),
  jetSrcCA15_( consumes<vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsCA15"))),
  genSrc_( consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>( "genParticles"))),
  stage2JetToken_(consumes<BXVector<l1t::Jet>>( edm::InputTag("caloStage2Digis","Jet","RECO"))),
  stage2TauToken_(consumes<BXVector<l1t::Tau>>( edm::InputTag("caloStage2Digis","Tau","RECO"))),
  l1BoostedToken_(consumes<vector<l1extra::L1JetParticle>>( edm::InputTag("uct2016EmulatorDigis","Boosted","")))
{
  // Initialize the Tree

  recoPt_      = iConfig.getParameter<double>("recoPtCut");
  nEvents      = tfs_->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  //recoJet_pt   = tfs_->make<TH1F>( "recoJet_pt" , "p_{t}", 300,  0., 300. );
  //recoJet_eta  = tfs_->make<TH1F>( "recoJet_eta"  , "eta", 100,  -3, 3. );
  //recoJet_phi  = tfs_->make<TH1F>( "recoJet_phi"  , "phi", 100,  -4, 4. );
  //
  //recoJetAK8_pt   = tfs_->make<TH1F>( "recoJetAK8_pt" , "p_{t}", 300,  0., 300. );
  //recoJetAK8_eta  = tfs_->make<TH1F>( "recoJetAK8_eta"  , "eta", 100,  -3, 3. );
  //recoJetAK8_phi  = tfs_->make<TH1F>( "recoJetAK8_phi"  , "phi", 100,  -4, 4. );

  efficiencyTree = tfs_->make<TTree>("efficiencyTree", "efficiencyTree");
  createBranches(efficiencyTree);
}

BoostedJetStudies::~BoostedJetStudies() {
}

//
// member functions
//

// ------------ method called to produce the data  ------------

void BoostedJetStudies::analyze( const edm::Event& evt, const edm::EventSetup& es )
{
  using namespace edm;

  nEvents->Fill(1);
  run = evt.id().run();
  lumi = evt.id().luminosityBlock();
  event = evt.id().event();
  Handle<L1CaloRegionCollection> regions;
   
  std::vector<reco::CaloJet> goodJets;
  std::vector<pat::Jet> goodJetsAK8;
  std::vector<pat::Jet> goodJetsCA15;
  std::vector<l1t::Jet> seeds;

  l1Jets->clear();
  seed180->clear();
  tauseed->clear();
  ak8Jets->clear();
  //ca15Jets->clear();
  subJets->clear();
  nSubJets.clear();
  nBHadrons.clear();
  subJetHFlav.clear();
  tau1.clear();
  tau2.clear();
  tau3.clear();

  // Accessing existing L1 seed stored in MINIAOD
  edm::Handle<BXVector<l1t::Jet>> stage2Jets;
  if(!evt.getByToken(stage2JetToken_, stage2Jets)) cout<<"ERROR GETTING THE STAGE 2 JETS"<<std::endl;
  evt.getByToken(stage2JetToken_, stage2Jets);
  const BXVector<l1t::Jet> &s2j = *stage2Jets;
  for(auto obj : s2j) {
    seeds.push_back(obj);
    TLorentzVector temp;
    temp.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    seed180->push_back(temp);
  }

  edm::Handle<BXVector<l1t::Tau>> stage2Taus;
  if(!evt.getByToken(stage2TauToken_, stage2Taus)) cout<<"ERROR GETTING THE STAGE 2 TAUS"<<std::endl;
  evt.getByToken(stage2TauToken_, stage2Taus);
  const BXVector<l1t::Tau> &s2t = *stage2Taus;
  for(auto obj : s2t) {
    TLorentzVector temp;
    temp.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    tauseed->push_back(temp);
  }

  // Accessing L1boosted collection
  edm::Handle<vector<l1extra::L1JetParticle>> l1Boosted;
  if(!evt.getByToken(l1BoostedToken_, l1Boosted)) cout<<"ERROR GETTING THE L1BOOSTED JETS"<<std::endl;
  evt.getByToken(l1BoostedToken_, l1Boosted);
  const vector<l1extra::L1JetParticle> &l1B = *l1Boosted;
  for(auto obj : l1B) {
    TLorentzVector temp;
    temp.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.pt());
    l1Jets->push_back(temp);
  }

  // Start Runing Analysis


  Handle<vector<reco::CaloJet> > jets;

  if(evt.getByToken(jetSrc_, jets)){//Begin Getting Calo Jets
    for (const reco::CaloJet &jet : *jets) {
      //recoJet_pt->Fill( jet.pt() );
      //recoJet_eta->Fill( jet.eta() );
      //recoJet_phi->Fill( jet.phi() );
      if(jet.pt() > recoPt_ ) {
	goodJets.push_back(jet);
      }
    }
  }
  else
    cout<<"Error getting calo jets"<<std::endl;

  Handle<vector<pat::Jet> > jetsAK8;

  if(evt.getByToken(jetSrcAK8_, jetsAK8)){//Begin Getting AK8 Jets
    for (const pat::Jet &jetAK8 : *jetsAK8) {
      if(jetAK8.pt() > recoPt_ ) {
        nSubJets.push_back(jetAK8.subjets("SoftDropPuppi").size());
        nBHadrons.push_back(jetAK8.jetFlavourInfo().getbHadrons().size());
        TLorentzVector temp ;
        temp.SetPtEtaPhiE(jetAK8.pt(),jetAK8.eta(),jetAK8.phi(),jetAK8.et());
        ak8Jets->push_back(temp);
        //if(jetAK8.subjets("SoftDropPuppi").size() ==  2 && jetAK8.jetFlavourInfo().getbHadrons().size() > 1){
        if(jetAK8.subjets("SoftDropPuppi").size() ==  2){  // for Zprime->qq sample 
          goodJetsAK8.push_back(jetAK8);
        }
      }
    }
  }
  else
    cout<<"Error getting AK8 jets"<<std::endl;

  Handle<vector<pat::Jet> > jetsCA15;

  if(evt.getByToken(jetSrcCA15_, jetsCA15)){//Begin Getting CA15 Jets
    for (const pat::Jet &jetCA15 : *jetsCA15) {
      if(jetCA15.pt() > recoPt_ ) {
        //TLorentzVector temp ;
        //temp.SetPtEtaPhiE(jetCA15.pt(),jetCA15.eta(),jetCA15.phi(),jetCA15.et());
        //ca15Jets->push_back(temp);
        if(jetCA15.subjets("SoftDrop").size() ==  2){
          goodJetsCA15.push_back(jetCA15);
        }
      }
    }
  }
  else
    cout<<"Error getting CA15 jets"<<std::endl;

  zeroOutAllVariables();

  vector<l1extra::L1JetParticle> l1JetsSorted;
  for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1Boosted->begin(); l1Jet != l1Boosted->end(); l1Jet++ ){
    l1JetsSorted.push_back(*l1Jet);
  }
  if(l1JetsSorted.size() > 1){  std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPt);}

  if(goodJets.size() > 0){
    reco::CaloJet recoJet;

    recoPt  = goodJets.at(0).pt();
    recoEta = goodJets.at(0).eta();
    recoPhi = goodJets.at(0).phi();
    recoJet = goodJets.at(0);

    int i = 0;
    int foundL1Jet = 0;
    if(l1JetsSorted.size() > 0){
      for(auto jet : l1JetsSorted){
        if(reco::deltaR(jet, recoJet) < 0.4 && foundL1Jet == 0 ){
          l1Pt  = jet.pt();
          l1Eta = jet.eta();
          l1Phi = jet.phi();
          l1NthJet = i;
          foundL1Jet = 1;
        }
        i++;
      }
    }
    int j = 0;
    int foundSeed = 0;
    if(seeds.size() > 0){
      for(auto seed : seeds){
        if(reco::deltaR(seed, recoJet) < 0.4 && foundSeed == 0 ){
          seedPt  = seed.pt();
          seedEta = seed.eta();
          seedPhi = seed.phi();
          seedNthJet = j;
          foundSeed = 1;
        }
        j++;
      }
    }

  }

  if(goodJetsAK8.size() > 0){

    for(auto jet:goodJetsAK8){
      tau1.push_back(jet.userFloat("NjettinessAK8Puppi:tau1"));
      tau2.push_back(jet.userFloat("NjettinessAK8Puppi:tau2"));
      tau3.push_back(jet.userFloat("NjettinessAK8Puppi:tau3"));
      HFlav.clear();
      for(unsigned int isub=0; isub<((jet.subjets("SoftDropPuppi")).size()); isub++){
        HFlav.push_back(jet.subjets("SoftDropPuppi")[isub]->hadronFlavour());
        TLorentzVector temp;
        temp.SetPtEtaPhiE(jet.subjets("SoftDropPuppi")[isub]->pt(),jet.subjets("SoftDropPuppi")[isub]->eta(),jet.subjets("SoftDropPuppi")[isub]->phi(),jet.subjets("SoftDropPuppi")[isub]->et());
        subJets->push_back(temp);
      }
      subJetHFlav.push_back(HFlav);
      //take more variables from here: https://github.com/gouskos/HiggsToBBNtupleProducerTool/blob/opendata_80X/NtupleAK8/src/FatJetInfoFiller.cc#L215-L217
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
    }

    //Match to boosted jets and see if we can match subjettiness functions...
    //vector<l1extra::L1JetParticle> l1JetsSorted;
    //for( vector<l1extra::L1JetParticle>::const_iterator l1Jet = l1Boosted->begin(); l1Jet != l1Boosted->end(); l1Jet++ ){
    //  l1JetsSorted.push_back(*l1Jet);
    //}
    //if(l1JetsSorted.size() > 1){  std::sort(l1JetsSorted.begin(),l1JetsSorted.end(),compareByPt);}

    pat::Jet recoJet_1;

    recoPt_1  = goodJetsAK8.at(0).pt();
    recoEta_1 = goodJetsAK8.at(0).eta();
    recoPhi_1 = goodJetsAK8.at(0).phi();
    recoMass_1 = goodJetsAK8.at(0).userFloat("ak8PFJetsPuppiSoftDropMass");
    recoJet_1 = goodJetsAK8.at(0);

    int i = 0;
    int foundL1Jet_1 = 0;
    if(l1JetsSorted.size() > 0){
      for(auto jet : l1JetsSorted){
        if(reco::deltaR(jet, recoJet_1) < 0.4 && foundL1Jet_1 == 0 ){
          l1Pt_1  = jet.pt();
          l1Eta_1 = jet.eta();
          l1Phi_1 = jet.phi();
          l1NthJet_1 = i;
          foundL1Jet_1 = 1;
        }
        i++;
      }
    }

    int j = 0;
    int foundSeed_1 = 0;
    if(seeds.size() > 0){
      for(auto seed : seeds){
        if(reco::deltaR(seed, recoJet_1) < 0.4 && foundSeed_1 == 0 ){
          seedPt_1  = seed.pt();
          seedEta_1 = seed.eta();
          seedPhi_1 = seed.phi();
          seedNthJet_1 = j;
          foundSeed_1 = 1;
        }
        j++;
      }
    }

  }

  if(goodJetsCA15.size()>0){
    pat::Jet recoJet_2;

    recoPt_2  = goodJetsCA15.at(0).pt();
    recoEta_2 = goodJetsCA15.at(0).eta();
    recoPhi_2 = goodJetsCA15.at(0).phi();
    recoJet_2 = goodJetsCA15.at(0);

    int i = 0;
    int foundL1Jet_2 = 0;
    if(l1JetsSorted.size() > 0){
      for(auto jet : l1JetsSorted){
        if(reco::deltaR(jet, recoJet_2) < 0.4 && foundL1Jet_2 == 0 ){
          l1Pt_2  = jet.pt();
          l1Eta_2 = jet.eta();
          l1Phi_2 = jet.phi();
          l1NthJet_2 = i;
          foundL1Jet_2 = 1;
        }
        i++;
      }
    }
    int j = 0;
    int foundSeed_2 = 0;
    if(seeds.size() > 0){
      for(auto seed : seeds){
        if(reco::deltaR(seed, recoJet_2) < 0.4 && foundSeed_2 == 0 ){
          seedPt_2  = seed.pt();
          seedEta_2 = seed.eta();
          seedPhi_2 = seed.phi();
          seedNthJet_2 = j;
          foundSeed_2 = 1;
        }
        j++;
      }
    }

  }

  edm::Handle<reco::GenParticleCollection> genParticles;
  if(evt.getByToken(genSrc_, genParticles)){//Begin Getting Gen Particles
    for (reco::GenParticleCollection::const_iterator genparticle = genParticles->begin(); genparticle != genParticles->end(); genparticle++){
      double DR_1 = reco::deltaR(recoEta_1, recoPhi_1, genparticle->eta(), genparticle->phi());
      if (DR_1 < genDR_1 && genparticle->status() > 21 && genparticle->status() < 41){
        genDR_1 = DR_1;
        genId_1 = genparticle->pdgId();
        genMother_1 = genparticle->motherRef(0)->pdgId();
        genPt_1 = genparticle->pt();
        genEta_1 = genparticle->eta();
        genPhi_1 = genparticle->phi();
        genM_1 = genparticle->mass();
      }
      double DR_2 = reco::deltaR(recoEta_2, recoPhi_2, genparticle->eta(), genparticle->phi());
      if (DR_2 < genDR_2 && genparticle->status() > 21 && genparticle->status() < 41){
        genDR_2 = DR_1;
        genId_2 = genparticle->pdgId();
        genMother_2 = genparticle->motherRef(0)->pdgId();
        genPt_2 = genparticle->pt();
        genEta_2 = genparticle->eta();
        genPhi_2 = genparticle->phi();
        genM_2 = genparticle->mass();
      }
    }
  }
  efficiencyTree->Fill();
}

void BoostedJetStudies::zeroOutAllVariables(){
  genPt_1=-99; genEta_1=-99; genPhi_1=-99; genM_1=-99; 
  genDR_1=99; genId_1=-99; genMother_1=-99;
  genPt_2=-99; genEta_2=-99; genPhi_2=-99; genM_2=-99; 
  genDR_2=99; genId_2=-99; genMother_2=-99;
  seedPt=-99; seedEta=-99; seedPhi=-99; seedNthJet=-99;
  seedPt_1=-99; seedEta_1=-99; seedPhi_1=-99; seedNthJet_1=-99;
  seedPt_2=-99; seedEta_2=-99; seedPhi_2=-99; seedNthJet_2=-99;
  recoPt=-99; recoEta=-99; recoPhi=-99;
  recoPt_1=-99; recoEta_1=-99; recoPhi_1=-99; recoMass_1=-99;
  recoPt_2=-99; recoEta_2=-99; recoPhi_2=-99;
  l1Pt=-99; l1Eta=-99; l1Phi=-99; l1NthJet=-99;
  l1Pt_1=-99; l1Eta_1=-99; l1Phi_1=-99; l1NthJet_1=-99;
  l1Pt_2=-99; l1Eta_2=-99; l1Phi_2=-99; l1NthJet_2=-99;
  recoPt_=-99;
}

void BoostedJetStudies::createBranches(TTree *tree){
    tree->Branch("run",     &run,     "run/I");
    tree->Branch("lumi",    &lumi,    "lumi/I");
    tree->Branch("event",   &event,   "event/I");
    tree->Branch("genPt_1",       &genPt_1,      "genPt_1/D");
    tree->Branch("genEta_1",      &genEta_1,     "genEta_1/D");
    tree->Branch("genPhi_1",      &genPhi_1,     "genPhi_1/D");
    tree->Branch("genM_1",        &genM_1,       "genM_1/D");
    tree->Branch("genDR_1",       &genDR_1,      "genDR_1/D");
    tree->Branch("genId_1",       &genId_1,      "genId_1/I");
    tree->Branch("genMother_1",   &genMother_1,  "genMother_1/I");
    tree->Branch("genPt_2",       &genPt_2,      "genPt_2/D");
    tree->Branch("genEta_2",      &genEta_2,     "genEta_2/D");
    tree->Branch("genPhi_2",      &genPhi_2,     "genPhi_2/D");
    tree->Branch("genM_2",        &genM_2,       "genM_2/D");
    tree->Branch("genDR_2",       &genDR_2,      "genDR_2/D");
    tree->Branch("genId_2",       &genId_2,      "genId_2/I");
    tree->Branch("genMother_2",   &genMother_2,  "genMother_2/I");
    tree->Branch("seedPt",        &seedPt,       "seedPt/D");
    tree->Branch("seedEta",       &seedEta,      "seedEta/D");
    tree->Branch("seedPhi",       &seedPhi,      "seedPhi/D");
    tree->Branch("seedNthJet",    &seedNthJet,   "seedNthJet/I");
    tree->Branch("seedPt_1",      &seedPt_1,     "seedPt_1/D");
    tree->Branch("seedEta_1",     &seedEta_1,    "seedEta_1/D");
    tree->Branch("seedPhi_1",     &seedPhi_1,    "seedPhi_1/D");
    tree->Branch("seedNthJet_1",  &seedNthJet_1, "seedNthJet_1/I");
    tree->Branch("seedPt_2",      &seedPt_2,     "seedPt_2/D");
    tree->Branch("seedEta_2",     &seedEta_2,    "seedEta_2/D");
    tree->Branch("seedPhi_2",     &seedPhi_2,    "seedPhi_2/D");
    tree->Branch("seedNthJet_2",  &seedNthJet_2, "seedNthJet_2/I");
    tree->Branch("recoPt",        &recoPt,       "recoPt/D");
    tree->Branch("recoEta",       &recoEta,      "recoEta/D");
    tree->Branch("recoPhi",       &recoPhi,      "recoPhi/D");
    tree->Branch("recoPt_1",      &recoPt_1,     "recoPt_1/D");
    tree->Branch("recoEta_1",     &recoEta_1,    "recoEta_1/D");
    tree->Branch("recoPhi_1",     &recoPhi_1,    "recoPhi_1/D");
    tree->Branch("recoMass_1",    &recoMass_1,   "recoMass_1/D");
    tree->Branch("recoPt_2",      &recoPt_2,     "recoPt_2/D");
    tree->Branch("recoEta_2",     &recoEta_2,    "recoEta_2/D");
    tree->Branch("recoPhi_2",     &recoPhi_2,    "recoPhi_2/D");
    tree->Branch("l1Pt",          &l1Pt,         "l1Pt/D");
    tree->Branch("l1Eta",         &l1Eta,        "l1Eta/D");
    tree->Branch("l1Phi",         &l1Phi,        "l1Phi/D");
    tree->Branch("l1NthJet",      &l1NthJet,     "l1NthJet/I");
    tree->Branch("l1Pt_1",        &l1Pt_1,       "l1Pt_1/D"); 
    tree->Branch("l1Eta_1",       &l1Eta_1,      "l1Eta_1/D");
    tree->Branch("l1Phi_1",       &l1Phi_1,      "l1Phi_1/D");
    tree->Branch("l1NthJet_1",    &l1NthJet_1,   "l1NthJet_1/I");
    tree->Branch("l1Pt_2",        &l1Pt_2,       "l1Pt_2/D");
    tree->Branch("l1Eta_2",       &l1Eta_2,      "l1Eta_2/D");
    tree->Branch("l1Phi_2",       &l1Phi_2,      "l1Phi_2/D");
    tree->Branch("l1NthJet_2",    &l1NthJet_2,   "l1NthJet_2/I");
    tree->Branch("tau1",          &tau1);
    tree->Branch("tau2",          &tau2);
    tree->Branch("tau3",          &tau3);
    tree->Branch("nSubJets",      &nSubJets);
    tree->Branch("subJetHFlav",   &subJetHFlav);
    tree->Branch("nBHadrons",     &nBHadrons);
    tree->Branch("l1Jets", "vector<TLorentzVector>", &l1Jets, 32000, 0);
    tree->Branch("seed180", "vector<TLorentzVector>", &seed180, 32000, 0);
    tree->Branch("tauseed", "vector<TLorentzVector>", &tauseed, 32000, 0);
    tree->Branch("ak8Jets", "vector<TLorentzVector>", &ak8Jets, 32000, 0);
    //tree->Branch("ca15Jets", "vector<TLorentzVector>", &ca15Jets, 32000, 0);
    tree->Branch("subJets", "vector<TLorentzVector>", &subJets, 32000, 0);
  }


// ------------ method called once each job just before starting event loop  ------------
void 
BoostedJetStudies::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BoostedJetStudies::endJob() {
}

// ------------ method called when starting to processes a run  ------------

void
BoostedJetStudies::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}
 
// ------------ method called when ending the processing of a run  ------------
/*
  void
  BoostedJetStudies::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  BoostedJetStudies::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  BoostedJetStudies::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BoostedJetStudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BoostedJetStudies);
