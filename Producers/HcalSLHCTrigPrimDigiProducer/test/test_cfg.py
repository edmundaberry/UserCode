import FWCore.ParameterSet.Config as cms

process = cms.Process("HcalUpgradeTest")

## Load setup python files
## process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.GlobalTag.globaltag = "STARTUP31X_V1::All"

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring( 'file:///tmp/eberry/hcalDigis.root' ))

## Load producers
process.load("Producers.HcalSLHCTrigPrimDigiProducer.hcalSlhcTpDigi_cfi")
process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(cms.InputTag('simHcalDigis'),cms.InputTag('simHcalDigis'))

# Load analyzer
process.load("Producers.HcalSLHCTrigPrimDigiProducer.hcalSlhcTpComp_cfi")

## Configure paths
process.PRODUCE_UPGRADE = cms.Path(process.simHcalSLHCTriggerPrimitiveDigis )
process.PRODUCE_DEFAULT = cms.Path(process.simHcalTriggerPrimitiveDigis )
process.ANALYZE = cms.Path (process.hcalSLHCTrigPrimDigiComparer )

## Schedule the paths
process.schedule = cms.Schedule()
process.schedule.append( process.PRODUCE_UPGRADE )
process.schedule.append( process.PRODUCE_DEFAULT )
process.schedule.append( process.ANALYZE  )

