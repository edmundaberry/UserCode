import FWCore.ParameterSet.Config as cms

process = cms.Process("HcalUpgradeTPG")

## Load setup python files
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "STARTUP31X_V1::All"

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )
process.load("Producers.HcalSLHCTrigPrimDigiProducer.RelValQCD_Pt_80_120_CMSSW_3_1_0_STARTUP31X_V1_GEN_SIM_DIGI_RAW_HLTDEBUG_cff")

## Load my producer
process.load("Producers.HcalSLHCTrigPrimDigiProducer.hcalSlhcTpDigi_cfi")

## Configure paths
process.TPG = cms.Path(process.simHcalSLHCTriggerPrimitiveDigis)

## Schedule the paths
process.schedule = cms.Schedule()
process.schedule.append( process.TPG )

