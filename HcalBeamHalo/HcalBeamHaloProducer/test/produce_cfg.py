
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("PRODUCE")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring (
                            'rfio:/castor/cern.ch/user/e/eberry/halo_with_digis.root'
))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (200) )

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration/StandardSequences/Geometry_cff")

process.hcalBeamHaloProducer = cms.EDProducer("HcalBeamHaloProducer",
                                              HBHERecHits = cms.InputTag ('hbhereco'),
                                              StandAloneTracks = cms.InputTag("cosmicMuonsEndCapsOnly"),
                                              Width = cms.int32 (5),
                                              MinRecHitEnergy = cms.double(0.5),
                                              MinWindowEnergy = cms.double(5.0) ,
                                              MinWindowCounts = cms.int32(3),
                                              MaxNWindows = cms.int32(3),
                                              Verbose = cms.bool (True) 
)

process.outputSkim = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('produce')),
    fileName = cms.untracked.string("/tmp/eberry/producer_output.root")
    )

process.produce = cms.Path ( process.hcalBeamHaloProducer )
process.output  = cms.EndPath ( process.outputSkim )


