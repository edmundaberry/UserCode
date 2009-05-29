import FWCore.ParameterSet.Config as cms

process = cms.Process("MyCaloTowerCreator")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/uscms/home/eberry/data/qcdDijets_ReRunL1_withCMSSW_2_2_3_ptBin9_jetThreshold20GeV_job1.root'
    )
)

## Ordinary cff files
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

## My CaloTower producer
process.load("Producers.CaloTowersFromTrigPrimsCreator.calotowersfromtrigprimscreator_cfi")

## My CaloTower comparator 
process.load("Analyzers.CaloTowerComp.calotowercomp_cfi")

## The default CaloTower producer
process.load("RecoJets.JetProducers.CaloTowerSchemeB_cfi")

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('data/final_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_simHcalTriggerPrimitiveDigis_*_*',
        'keep *_simEcalTriggerPrimitiveDigis_*_*',
        'keep CaloTowersSorted_*_*_*'
        )
)

process.calolocalreco_noEcalPreshower = cms.Sequence(process.ecalLocalRecoSequence_nopreshower+process.hcalLocalRecoSequence)

process.Unpack      = cms.Path(process.RawToDigi)
process.Reconstruct = cms.Path(process.calolocalreco_noEcalPreshower)
process.Produce     = cms.Path(process.towerMaker + process.caloTowersFromTrigPrimsCreator)
process.Analyze     = cms.Path(process.caloTowerComparator)
process.SaveAll     = cms.EndPath(process.out)

process.schedule = cms.Schedule()
process.schedule.append(process.Unpack)
process.schedule.append(process.Reconstruct)
process.schedule.append(process.Produce)
process.schedule.append(process.Analyze)
## process.schedule.append(process.SaveAll)

