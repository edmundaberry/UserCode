import FWCore.ParameterSet.Config as cms

process = cms.Process ("AnalyzeMuonPt")

## Load setup python files
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "STARTUP31X_V1::All"

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("Analyzers.CSCTFMuonPtAnalyzer.MyParticleGunMC_cfi")

## Load my analyzer
process.load("Analyzers.CSCTFMuonPtAnalyzer.csctfmuonptanalyzer_cfi")

## Configure paths
process.ANALYZE = cms.Path ( process.csctfMuonPtAnalyzer )
process.csctfMuonPtAnalyzer.NPtBins = 10

## Schedule paths
process.schedule = cms.Schedule()
process.schedule.append( process.ANALYZE )
