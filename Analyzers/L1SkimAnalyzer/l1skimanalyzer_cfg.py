import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
        'file:///scratch/eberry/L1SkimRaw_ppex_1.root'
    )
)

process.demo = cms.EDAnalyzer('L1SkimAnalyzer',
                              l1CaloEmCandsTag = cms.untracked.InputTag('gctDigis')
)


process.p = cms.Path(process.demo)
