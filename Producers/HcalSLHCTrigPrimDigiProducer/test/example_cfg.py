import FWCore.ParameterSet.Config as cms

process = cms.Process("HcalUpgradeTPG")

## Load setup python files
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )
process.source = cms.Source ("PoolSource",
                     fileNames = cms.untracked.vstring("file:///uscms_data/d2/eberry/SLHC/slhc_pionGun_DIGI.root"))

## Load my producer
process.load("Producers.HcalSLHCTrigPrimDigiProducer.hcalSlhcTpDigi_cfi")
process.simHcalSLHCTriggerPrimitiveDigis.SLHCMode = True

## Configure paths
process.TPG = cms.Path(process.simHcalSLHCTriggerPrimitiveDigis)

## Schedule the paths
process.schedule = cms.Schedule()
process.schedule.append( process.TPG )

