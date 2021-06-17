import FWCore.ParameterSet.Config as cms

l1NtupleProducer = cms.EDAnalyzer("BoostedJetStudies",                                  
                                  recoPtCut               = cms.double(10),
                                  recoJets                = cms.InputTag("slimmedCaloJets"),
                                  recoJetsAK8             = cms.InputTag("slimmedJetsAK8"),
                                  recoJetsCA15            = cms.InputTag("packedPatJetsCA15PFCHSSoftDrop"),
                                  genParticles            = cms.InputTag("genParticles", "", "HLT"),
)
