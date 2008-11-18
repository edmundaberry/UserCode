import FWCore.ParameterSet.Config as cms

process = cms.Process("unpackGctDigis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'CRAFT_V3P::All'

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/006540F5-32A8-DD11-A0C1-000423D98EC4.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/007F2BE0-51A8-DD11-9B1D-000423D98750.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/00F5DCD0-43A8-DD11-B7FF-000423D98BE8.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/0276CF62-55A8-DD11-9535-000423D98868.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/0467F071-55A8-DD11-831D-000423D991D4.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/04699363-28A8-DD11-9C97-001D09F23174.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/0855BCA1-3FA8-DD11-8408-0019B9F7310E.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/087E67AC-38A8-DD11-A4FC-000423D944FC.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/08B7E2A3-46A8-DD11-9885-000423D99BF2.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/0A208F55-5AA8-DD11-A3F8-0030487A3C9A.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/0AA3AC89-44A8-DD11-840C-0030487A3232.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/0AAEA350-53A8-DD11-937D-000423D6B2D8.root',
        '/store/data/Commissioning08/Calo/RAW/v1/000/068/288/0C3C192F-56A8-DD11-A0BD-000423D9989E.root'
        )
)

## Unpack HLT info from FED
process.load("HLTrigger.HLTanalyzers.HLTopen_cff")

## Unpack GCT info from RAW
process.load( "EventFilter/GctRawToDigi/l1GctHwDigis_cfi" )
process.l1GctHwDigis.inputLabel = cms.InputTag( "source" )
process.DoGCT = cms.Path( process.l1GctHwDigis )

## Analyze
process.analyzer = cms.EDAnalyzer('L1SkimAnalyzer')
process.analyze = cms.Path( process.analyzer )

## Schedule
process.schedule = cms.Schedule(process.DoGCT, process.analyze)
