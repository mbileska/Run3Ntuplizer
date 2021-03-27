import FWCore.ParameterSet.Config as cms

l1NtupleProducer = cms.EDAnalyzer("BoostedJetStudies",                                  
                                  recoPtCut               = cms.double(10),
                                  recoJets                = cms.InputTag("slimmedJets"),
                                  recoJetsAK8             = cms.InputTag("slimmedJetsAK8"),
                                  genParticles            = cms.InputTag("genParticles", "", "HLT"),
)
