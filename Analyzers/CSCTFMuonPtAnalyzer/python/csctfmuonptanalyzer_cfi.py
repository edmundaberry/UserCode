import FWCore.ParameterSet.Config as cms

csctfMuonPtAnalyzer = cms.EDAnalyzer('CSCTFMuonPtAnalyzer',
  CSCCorrLCTDigiTag = cms.InputTag("simCscTriggerPrimitiveDigis","MPCSORTED"),
  DTTPDigiTag = cms.InputTag("simDtTriggerPrimitiveDigis"),
  GenParticlesTag = cms.InputTag("genParticles"),
  SimHitsTag = cms.InputTag("g4SimHits","MuonCSCHits"),                                     
  SimTracksTag = cms.InputTag("g4SimHits"),
  AnalyzeGenMuons = cms.bool(True),                                  
  AnalyzeDataMuons = cms.bool(False),
  Verbose = cms.bool (False),                                                                          
  NPtBins = cms.int32(140),
  MinBinnablePt = cms.double(0),
  MaxBinnablePt = cms.double(140),                             
  NEtaBins = cms.int32(16),
  MinBinnableEta = cms.double(0.8),
  MaxBinnableEta = cms.double(2.1),
  FileName = cms.untracked.string("fileName.root"),
  DPhiHistName = cms.untracked.string("dPhiHistName.root"),
  DEtaHistName = cms.untracked.string("dEtaHistName.root"),
  ScalePt = cms.vdouble(-1.0, 0.0, 1.5, 2.0, 2.5,
                         3.0, 3.5, 4.0, 4.5, 5.0,
                         6.0, 7.0, 8.0, 10.0, 12.0,
                         14.0, 16.0, 18.0, 20.0, 25.0,
                         30.0, 35.0, 40.0, 45.0, 50.0,
                         60.0, 70.0, 80.0, 90.0, 100.0,
                         120.0, 140.0, 1000000.0),         
  ScaleEta = cms.vdouble ( -1.0, 0.89, 1.06, 1.14, 1.21, 1.26, 
                            1.32, 1.37, 1.41, 1.46, 1.51, 1.56, 
                            1.60, 1.65, 1.71, 2.40, 1000.0 )                                      
)
