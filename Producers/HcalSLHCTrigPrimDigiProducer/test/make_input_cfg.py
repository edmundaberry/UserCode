import FWCore.ParameterSet.Config as cms

process = cms.Process("HcalUpgradeMake")

## Load setup python files
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.GlobalTag.globaltag = "STARTUP31X_V1::All"

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.load("Producers.HcalSLHCTrigPrimDigiProducer.RelValQCD_Pt_80_120_CMSSW_3_1_0_STARTUP31X_V1_GEN_SIM_DIGI_RAW_HLTDEBUG_cff")

## Configure the output
process.dump = cms.OutputModule("PoolOutputModule",                              
   outputCommands = cms.untracked.vstring(
        'drop *',
        'keep HBHEDataFramesSorted_*_*_*',
        'keep HFDataFramesSorted_*_*_*',
        'keep HODataFramesSorted_*_*_*'
        ),
   fileName = cms.untracked.string('/tmp/eberry/hcalDigis.root')
)

## Configure paths
process.DUMP    = cms.EndPath (process.dump)

## Schedule the paths
process.schedule = cms.Schedule()
process.schedule.append( process.DUMP )
