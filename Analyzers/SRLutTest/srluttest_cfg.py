import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0010/A6A47EA6-0D84-DE11-8311-001D09F23F2A.root'
    )
)

process.demo = cms.EDAnalyzer('SRLutTest')

process.p = cms.Path(process.demo)
