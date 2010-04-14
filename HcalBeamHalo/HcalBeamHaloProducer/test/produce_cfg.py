
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("PRODUCE")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring (
                            'rfio:/castor/cern.ch/user/e/eberry/halo_with_digis.root'
))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (200) )

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR10_P_V4::All'

process.hcalBeamHaloProducer = cms.EDProducer("HcalBeamHaloProducer",
                                              HBHERecHits = cms.InputTag ('hbhereco'),
                                              StandAloneTracks = cms.InputTag("cosmicMuonsEndCapsOnly"),
                                              Width = cms.int32 (5),
                                              MinRecHitEnergy = cms.double(0.5),
                                              MinWindowEnergy = cms.double(5.0) ,
                                              MinWindowCounts = cms.int32(3),
                                              MaxNWindows = cms.int32(3),
                                              Verbose = cms.bool (False),
                                              BlackListCells = cms.vint32 (311401)
)

process.hcalBeamHaloAnalyzer = cms.EDAnalyzer("HcalBeamHaloAnalyzer",
                                              HBHERecHits = cms.InputTag ('hbhereco'),
                                              CSCRecHits = cms.InputTag("csc2DRecHits"),
                                              CSCSegments = cms.InputTag("cscSegments"),
                                              StandAloneTracks = cms.InputTag("cosmicMuonsEndCapsOnly"),
                                              HcalBeamHalo = cms.InputTag ("hcalBeamHaloProducer"),
                                              NtupleFileName = cms.string ("analyzer_output.root"),
                                              MinRecHitEnergy = process.hcalBeamHaloProducer.MinRecHitEnergy
)

process.outputSkim = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('produce')),
    fileName = cms.untracked.string("/tmp/eberry/producer_output.root")
    )

process.produce = cms.Path ( process.hcalBeamHaloProducer )
process.analyze = cms.Path ( process.hcalBeamHaloAnalyzer ) 
# process.output  = cms.EndPath ( process.outputSkim )


