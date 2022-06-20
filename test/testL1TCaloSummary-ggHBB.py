import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloSummaryTest")

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing()
options.register('runNumber', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.register('lumis', '1-max', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Lumis')
options.register('dataStream', '/ExpressPhysics/Run2015D-Express-v4/FEVT', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Dataset to look for run in')
options.register('inputFiles', [], VarParsing.multiplicity.list, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('useORCON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Use ORCON for conditions.  This is necessary for very recent runs where conditions have not propogated to Frontier')
options.parseArguments()

def formatLumis(lumistring, run) :
    lumis = (lrange.split('-') for lrange in lumistring.split(','))
    runlumis = (['%d:%s' % (run,lumi) for lumi in lrange] for lrange in lumis)
    return ['-'.join(l) for l in runlumis]

print 'Getting files for run %d...' % options.runNumber

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load('L1Trigger.Configuration.SimL1Emulator_cff')

process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')
process.load('L1Trigger.L1TCaloLayer1.simCaloStage2Layer1Digis_cfi')
process.simCaloStage2Layer1Digis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.simCaloStage2Layer1Digis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")

process.load('L1Trigger.L1TCaloLayer1.uct2016EmulatorDigis_cfi')

process.load("L1Trigger.Run3Ntuplizer.l1BoostedJetStudies_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#				'file:/eos/user/p/pdas/L1Boosted/ggHbb/MiniAOD/RunIIAutumn18MiniAOD_21Dec_0_5300.root'
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part1_757457_762757.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part2_8574_13874.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part3_5358_10658.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part4_4252_9552.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part5_2764_8064.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part6_9684_14984.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part7_1356_6656.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part8_6463_11763.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part9_5369_10669.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part10_9732_15032.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part11_5434_10734.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part12_8654_13954.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/MiniAOD/RunIIAutumn18MiniAOD_part13_9643_14943.root',
#				"root://cms-xrd-global.cern.ch://store/mc/RunIIAutumn18MiniAOD/ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/6340686B-CDEC-5E42-A004-BD829CDD00AF.root"
),
                            secondaryFileNames = cms.untracked.vstring(
#				'file:/eos/user/p/pdas/L1Boosted/ggHbb/DR/RunIIAutumn18DRPremix_step1_21Dec_0_5300.root'
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part1_757457_762757.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part2_8574_13874.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part3_5358_10658.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part4_4252_9552.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part5_2764_8064.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part6_9684_14984.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part7_1356_6656.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part8_6463_11763.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part9_5369_10669.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part10_9732_15032.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part11_5434_10734.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part12_8654_13954.root',
                                'file:/eos/user/p/pdas/L1Boosted/ZH/DR/RunIIAutumn18DRPremix_step1_part13_9643_14943.root',
#				"root://cms-xrd-global.cern.ch://store/mc/RunIIFall18wmLHEGS/ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8/GEN-SIM/102X_upgrade2018_realistic_v11-v1/110000/71D75AEE-FDC7-D843-BAC5-BD0D9B6DDB02.root"
                            )
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:734","1:961","1:966","1:982")
#process.source.eventsToProcess = cms.untracked.VEventRange("1:960110","1:965580")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("l1TFullEvent.root"),
    outputCommands = cms.untracked.vstring('keep *')
    #outputCommands = cms.untracked.vstring('drop *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)


#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple-ZH-ak8.root")
)

process.p = cms.Path(process.l1tCaloLayer1Digis*process.simCaloStage2Layer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)

process.e = cms.EndPath(process.out)

#process.schedule = cms.Schedule(process.p,process.e)
process.schedule = cms.Schedule(process.p)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

# Multi-threading
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)

# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
