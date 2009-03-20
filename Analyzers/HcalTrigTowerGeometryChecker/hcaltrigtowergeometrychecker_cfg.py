import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")

process.demo = cms.EDAnalyzer('HcalTrigTowerGeometryChecker',
                              verbose              = cms.untracked.bool(True),
                              check_hb_doublecount = cms.untracked.bool(True),
                              check_he_doublecount = cms.untracked.bool(True),
                              check_ho_doublecount = cms.untracked.bool(True),
                              check_hf_doublecount = cms.untracked.bool(True)
)

process.p = cms.Path(process.demo)
