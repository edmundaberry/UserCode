
import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA2")

## Load standard sequences and event content
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.EventContent.EventContent_cff")

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2))
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring (
   'file:///uscms/home/eberry/data/output.root'
))

## Produce CaloTowers from trigger primitives
process.load("Producers.CaloTowersFromTrigPrimsCreator.calotowersfromtrigprimscreator_cfi")

## Declare paths 
process.TOWERS = cms.Path ( process.caloTowersFromTrigPrimsCreator)

## Schedule paths
process.schedule = cms.Schedule()
process.schedule.append(process.TOWERS)

