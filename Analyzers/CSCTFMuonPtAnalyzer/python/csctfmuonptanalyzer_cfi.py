import FWCore.ParameterSet.Config as cms

csctfMuonPtAnalyzer = cms.EDAnalyzer('CSCTFMuonPtAnalyzer',
  CSCCorrLCTDigiTag = cms.InputTag("simCscTriggerPrimitiveDigis","MPCSORTED"),
  GenParticlesTag = cms.InputTag("genParticles"),
  AnalyzeGenMuons = cms.bool(True),                                  
  AnalyzeDataMuons = cms.bool(False),                                     
  NPtBins = cms.int32(140),
  MinBinnablePt = cms.double(0),
  MaxBinnablePt = cms.double(140)                             
)
