
import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA2")

## Load standard sequences and event content
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.EventContent.EventContent_cff")

## Load timing service
process.PathTimerService = cms.Service( "PathTimerService" )
process.timer = cms.EDProducer( "PathTimerInserter" )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring (
   'file:///uscms/home/eberry/data/output_skimmed.root'
))

## Produce CaloTowers from trigger primitives
process.load("Producers.CaloTowersFromTrigPrimsCreator.calotowersfromtrigprimscreator_cfi")
process.caloTowersFromTrigPrimsCreator.hcalTrigPrimTag = "caloTrigPrimsSkimmer"
process.caloTowersFromTrigPrimsCreator.ecalTrigPrimTag = "caloTrigPrimsSkimmer"

## Declare paths 
process.TOWERS = cms.Path    ( process.caloTowersFromTrigPrimsCreator)
process.TIMER  = cms.EndPath ( process.timer )

## Schedule paths
process.schedule = cms.Schedule()
process.schedule.append(process.TOWERS)
process.schedule.append(process.TIMER )
