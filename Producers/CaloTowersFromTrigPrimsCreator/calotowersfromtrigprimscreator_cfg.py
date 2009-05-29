import FWCore.ParameterSet.Config as cms

process = cms.Process("MyCaloTowerCreator")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("Producers.CaloTowersFromTrigPrimsCreator.calotowersfromtrigprimscreator_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/uscms/home/eberry/data/test.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('data/final_output.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_simHcalTriggerPrimitiveDigis_*_*',
        'keep *_simEcalTriggerPrimitiveDigis_*_*',
        'keep CaloTowersSorted_*_*_*'
        )
)

process.unpack = cms.Path(process.RawToDigi)
process.produce = cms.Path(process.caloTowersFromTrigPrimsCreator)
process.finalProcess = cms.EndPath(process.out)

process.schedule = cms.Schedule()
##process.schedule.append(process.unpack)
process.schedule.append(process.produce)
process.schedule.append(process.finalProcess)


