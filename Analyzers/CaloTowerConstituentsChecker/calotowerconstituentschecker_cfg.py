import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")

process.demo = cms.EDAnalyzer('CaloTowerConstituentsChecker',
                              verbose              = cms.untracked.bool(True),
                              check_eb_doublecount = cms.untracked.bool(True),
                              check_ee_doublecount = cms.untracked.bool(True)
)

process.p = cms.Path(process.demo)
