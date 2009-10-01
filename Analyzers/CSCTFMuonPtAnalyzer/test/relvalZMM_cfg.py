import FWCore.ParameterSet.Config as cms

process = cms.Process ("AnalyzeMuonPt")

## Load setup python files
process.load("FWCore.MessageService.MessageLogger_cfi")
##process.load("Configuration.StandardSequences.Geometry_cff")
##process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
##process.GlobalTag.globaltag = "STARTUP31X_V4::All"

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
   '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0010/A6A47EA6-0D84-DE11-8311-001D09F23F2A.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/EAE979E9-FF83-DE11-BBD8-0019B9F7312C.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/C86F40F6-FF83-DE11-9922-001617C3B65A.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/C4C40134-0084-DE11-9E81-001D09F27003.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/A2BCFD00-FD83-DE11-94A1-001617E30D12.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/9A645EE7-FC83-DE11-A3C7-001617C3B706.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/98F92016-0084-DE11-BCD3-0019B9F709A4.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/8EA3182C-FD83-DE11-B276-001D09F2424A.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/741DC3C6-FD83-DE11-9C92-001617DC1F70.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/4EAA0122-FD83-DE11-BEA1-001D09F28EA3.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/446B17AA-FD83-DE11-BF21-001D09F2423B.root',
   ## '/store/relval/CMSSW_3_2_4/RelValZMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V4-v1/0009/2A3D9AD8-FD83-DE11-8649-001D09F24D8A.root' 
))

## Load my analyzer
process.load("Analyzers.CSCTFMuonPtAnalyzer.csctfmuonptanalyzer_cfi")
process.csctfMuonPtAnalyzer.NPtBins = 10
process.csctfMuonPtAnalyzer.FileName = cms.string("RelValFile.root")
process.csctfMuonPtAnalyzer.HistName = cms.string("RelValHist.root")

## Configure paths
process.ANALYZE = cms.Path ( process.csctfMuonPtAnalyzer )

## Schedule paths
process.schedule = cms.Schedule()
process.schedule.append( process.ANALYZE )
