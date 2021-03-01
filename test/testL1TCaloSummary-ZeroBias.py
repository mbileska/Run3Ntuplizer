import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TCaloSummaryTest")

#import EventFilter.L1TXRawToDigi.util as util

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
#if len(options.inputFiles) is 0 and options.inputFileList is '' :
#    inputFiles = util.getFilesForRun(options.runNumber, options.dataStream)
#elif len(options.inputFileList) > 0 :
#    with open(options.inputFileList) as f :
#        inputFiles = list((line.strip() for line in f))
#else :
#    inputFiles = cms.untracked.vstring(options.inputFiles)
#if len(inputFiles) is 0 :
#    raise Exception('No files found for dataset %s run %d' % (options.dataStream, options.runNumber))
#print 'Ok, time to analyze'


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v10', '')

process.hcalDigis.saveQIE10DataNSamples = cms.untracked.vint32( 10)
process.hcalDigis.saveQIE10DataTags = cms.untracked.vstring( "MYDATA" )

# To get L1 CaloParams
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_cfi')
# To get CaloTPGTranscoder
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

#process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

process.load("L1Trigger.Run3Ntuplizer.l1BoostedJetStudies_cfi")

#process.l1NtupleProducer.isData = cms.bool(False)
process.l1NtupleProducer.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.l1NtupleProducer.hcalToken = cms.InputTag("l1tCaloLayer1Digis")
process.l1NtupleProducer.useECALLUT = cms.bool(False)
process.l1NtupleProducer.useHCALLUT = cms.bool(False)
process.l1NtupleProducer.useHFLUT   = cms.bool(False)
process.l1NtupleProducer.useLSB     = cms.bool(True)
process.l1NtupleProducer.verbose    = cms.bool(False)
process.l1NtupleProducer.activityFraction12 = cms.double(0.00390625)

process.uct2016EmulatorDigis.useECALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHFLUT = cms.bool(False)
process.uct2016EmulatorDigis.useLSB = cms.bool(True)
process.uct2016EmulatorDigis.verbose = cms.bool(False)
process.uct2016EmulatorDigis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.uct2016EmulatorDigis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#sourceFileList = open("inputFileList-MINI.txt","r")
#secondaryFileList = open("inputFileList.txt","r")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/257FD10B-6D35-7441-BA80-EC6F5205715A.root',
					'/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/45960F97-D686-9C48-9217-18B01B2C2171.root',
					'/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/51C2690F-6EE4-FB4A-9085-59E2EEDB915A.root',
					'/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/C8C9C76B-C6CD-3B41-A187-0F476F68A8F3.root',
					'/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/6541FFB6-9673-2641-9C1D-AAFE0CE36D1E.root',
					'/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/7A4F1CEA-9005-EB47-B4A0-2C9163205D6D.root',
					'/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/A92D9A3F-848E-8949-94AB-51FA607F3095.root',
					'/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/77E4AD09-EFC0-2D41-8A72-D9A58D7B8611.root',
					'/store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/120000/B3A2CADC-75F3-A242-9060-A290CDE14221.root'),
    secondaryFileNames = cms.untracked.vstring('/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/881D06AE-2D89-E811-B2AC-FA163E8F9107.root',
					'/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/14B165C7-2B89-E811-B118-FA163E31F8EA.root',
					'/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/70D101BD-2289-E811-BD8B-FA163E308A6C.root',
					'/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/70032A26-1989-E811-961F-FA163E354E33.root',
					'/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/62996682-2589-E811-A753-02163E017EC6.root',
					'/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/6868D96E-1489-E811-8480-FA163E50A19F.root',
					'/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/E46E586C-1589-E811-ADE3-FA163EF3AACF.root',
					'/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/EC057557-2389-E811-8205-02163E010C3C.root',
					'/store/data/Run2018C/ZeroBias/RAW/v1/000/319/756/00000/DC1A7C29-1389-E811-832B-FA163E98E77A.root')

#                             fileNames = cms.untracked.vstring(sourceFileList),
#                             secondaryFileNames = cms.untracked.vstring(secondaryFileList)
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("")
#process.source.eventsToProcess = cms.untracked.VEventRange("319756:991181466","319756:928807791","319756:608032744","319756:309750363","319756:701995845","319756:181355847","319756:192667981","319756:192601340","319756:644859036","319756:110119890")

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
)


#Output
process.TFileService = cms.Service(
	"TFileService",
	fileName = cms.string("l1TNtuple-ZeroBias.root")
)


process.L1TRawToDigi_Stage2 = cms.Task(process.caloLayer1Digis, process.caloStage2Digis)
process.RawToDigi_short = cms.Sequence(process.L1TRawToDigi_Stage2)
process.p = cms.Path(process.RawToDigi_short*process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)
process.e = cms.EndPath(process.out)

#process.schedule = cms.Schedule(process.p,process.e)
process.schedule = cms.Schedule(process.p)

#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)

# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
