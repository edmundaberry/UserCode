import FWCore.ParameterSet.Config as cms

process = cms.Process("TestCaloTowerAna")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.EventContent.EventContent_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'STARTUP_30X::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:///uscms_data/d1/mikeh/CRAFT_on_QCD_NewHcal_RECO_v7.root')
)

process.caloTowerAna = cms.EDAnalyzer('CaloTowerAna')

process.analyze = cms.Path(process.caloTowerAna)

process.schedule = cms.Schedule()
process.schedule.append(process.analyze)
