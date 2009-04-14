
import FWCore.ParameterSet.Config as cms

process = cms.Process("L1ReRun")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Analyzers.L1SkimAnalyzer.qcdDijet4_cff")
process.source.skipEvents = cms.untracked.uint32(0)

## Use digis produced by the re-run L1 emulator, not the original emulator
process.source.inputCommands = cms.untracked.vstring(
    'keep *',
    'drop *_simCscTriggerPrimitiveDigis_*_HLT',
    'drop *_simRctDigis_*_HLT',
    'drop *_simGctDigis_*_HLT',
    'drop *_simDtTriggerPrimitiveDigis_*_HLT',
    'drop *_simDttfDigis_*_HLT',
    'drop *_simGmtDigis_*_HLT',
    'drop *_simCsctfDigis_*_HLT',
    'drop *_simRpcTriggerDigis_*_HLT',
    'drop *_simCsctfTrackDigis_*_HLT',
    'drop *_TriggerResults_*_HLT',
    'drop *_l1extraParticles_*_HLT',
    'drop *_hlt*_*_HLT'
)



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

##-------------------------------------------------------------------------------------
## Include geometry and global tag
##-------------------------------------------------------------------------------------

## Normal global tag stuff

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
process.GlobalTag.globaltag = 'IDEAL_V9::All'

## Make sure that GlobalTag doesn't re-write your Jet calibration function
## Jet calibration function is where the HT jet threshold is set
## You need this to do L1 jet threshold studies

process.myPrefer = cms.ESPrefer("L1GctConfigProducers")
#process.SiPrefer = cms.ESPrefer("SiStripPedestalsFakeESSource")

##-------------------------------------------------------------------------------------
## Unpack raw data
##-------------------------------------------------------------------------------------

process.load("Configuration.StandardSequences.RawToDigi_cff")

##-------------------------------------------------------------------------------------
## Run trigger primitive generation on unpacked digis
##-------------------------------------------------------------------------------------

process.load("L1Trigger.Configuration.CaloTriggerPrimitives_cff")
process.load("L1Trigger.Configuration.SimL1Emulator_cff")
process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
process.simHcalTriggerPrimitiveDigis.inputLabel = 'hcalDigis'
process.simDtTriggerPrimitiveDigis.digiTag = 'muonDTDigis'
process.simCscTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag("muonCSCDigis","MuonCSCComparatorDigi")
process.simCscTriggerPrimitiveDigis.CSCWireDigiProducer = cms.InputTag("muonCSCDigis","MuonCSCWireDigi")
process.simRpcTriggerDigis.label = 'muonRPCDigis'

##-------------------------------------------------------------------------------------
## Run L1Extra particles on the emulator digis
##-------------------------------------------------------------------------------------

process.load('L1Trigger.Configuration.L1StartupConfig_cff')

process.L1GctConfigProducers.L1CaloJetZeroSuppressionThresholdInGeV = cms.double(1.0)

# Do L1Extra with re-run digis
process.l1extraParticles = cms.EDProducer("L1ExtraParticlesProd",
    muonSource = cms.InputTag("simGmtDigis","","L1ReRun"),
    etTotalSource = cms.InputTag("simGctDigis","","L1ReRun"),
    nonIsolatedEmSource = cms.InputTag("simGctDigis","nonIsoEm","L1ReRun"),
    etMissSource = cms.InputTag("simGctDigis","","L1ReRun"),
    produceMuonParticles = cms.bool(True),
    forwardJetSource = cms.InputTag("simGctDigis","forJets","L1ReRun"),
    centralJetSource = cms.InputTag("simGctDigis","cenJets","L1ReRun"),
    produceCaloParticles = cms.bool(True),
    tauJetSource = cms.InputTag("simGctDigis","tauJets","L1ReRun"),
    isolatedEmSource = cms.InputTag("simGctDigis","isoEm","L1ReRun"),
    etHadSource = cms.InputTag("simGctDigis","","L1ReRun"),
    centralBxOnly = cms.bool(True)
)

process.load('L1Trigger.L1ExtraFromDigis.l1extraParticleMap_cfi')

##-------------------------------------------------------------------------------------
## Load my analyzer
##-------------------------------------------------------------------------------------

process.load('Analyzers.L1SkimAnalyzer.l1skimanalyzer_Data_cfi')
process.l1SkimAnalyzerData.outSuffix = "_L1EmulatorOnMC_PtBin4_1GeV_job1_10events_TEST"
process.l1SkimAnalyzerData.outPath = cms.untracked.string("./")
process.analyze = cms.Path(process.l1SkimAnalyzerData)

##-------------------------------------------------------------------------------------
## If you need it for debugging, dump everything
##-------------------------------------------------------------------------------------

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring('keep *')
)
process.finalProcess = cms.EndPath(process.out)

##-------------------------------------------------------------------------------------
##Set up my L1 path
##-------------------------------------------------------------------------------------

process.myL1Path = cms.Path(
    process.ecalDigis
    *process.hcalDigis
    *process.muonDTDigis
    *process.muonCSCDigis
    *process.muonRPCDigis
    *process.CaloTriggerPrimitives
    *process.SimL1Emulator
    *process.l1extraParticles
) 

process.myTriggerWordPath = cms.Path(process.l1extraParticleMap)

##-------------------------------------------------------------------------------------
## Load the official HLT analyzer
##-------------------------------------------------------------------------------------

process.load("HLTrigger.HLTanalyzers.HLTAnalyser_cfi")
process.hltAnalyzer = cms.Path( process.hltanalysis )

##-------------------------------------------------------------------------------------
## HLT stuff
##-------------------------------------------------------------------------------------

## Needed for HLT
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

## Load the typical HLT cfg files
process.load("HLTrigger.HLTanalyzers.HLTopen_cff")
from Analyzers.L1SkimAnalyzer import patchToRerunL1Emulator

## Make sure the HLT uses the emulator digis
## patchToRerunL1Emulator.switchToSimGtDigis( process )

##-------------------------------------------------------------------------------------
## Schedule everything
##-------------------------------------------------------------------------------------

process.schedule = cms.Schedule()

process.schedule.append(process.myL1Path)
## process.schedule.append(process.myTriggerWordPath)
## process.schedule.append(process.DoHLTJets)
process.schedule.append(process.analyze)
process.schedule.append(process.finalProcess)
