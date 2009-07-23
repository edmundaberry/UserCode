
import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

## Load standard sequences and event content
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring (
   'file:///uscms/home/eberry/data/file.root'
))

## Reco
process.load("RecoLocalCalo.Configuration.RecoLocalCalo_cff")
process.calolocalreco_noEcalPreshower = cms.Sequence(process.ecalLocalRecoSequence_nopreshower+process.hcalLocalRecoSequence)

## Dump everything
process.out = cms.OutputModule("PoolOutputModule",                              
                               outputCommands = cms.untracked.vstring('keep *'),
                               fileName = cms.untracked.string('/uscms/home/eberry/data/output.root'))


## Declare paths
process.UNPACK = cms.Path(process.RawToDigi)
process.CALORECO = cms.Path(process.calolocalreco)
process.OUTPUT = cms.EndPath(process.out)


## Schedule paths
process.schedule = cms.Schedule()
process.schedule.append(process.UNPACK)
process.schedule.append(process.CALORECO)
process.schedule.append(process.OUTPUT)


