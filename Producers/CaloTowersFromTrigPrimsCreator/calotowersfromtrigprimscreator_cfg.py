import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")

## process.source = cms.Source("EmptySource")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/uscms/home/eberry/data/4406CC7F-3497-DD11-B5F5-001D0967DA99.root'
    )
)

process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

process.unpack = cms.Path(process.RawToDigi)

process.myProducerLabel = cms.EDProducer('CaloTowersFromTrigPrimsCreator')

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('full_output.root'),
    outputCommands = cms.untracked.vstring('keep *')
)

process.produce = cms.Path(process.myProducerLabel)
process.finalProcess = cms.EndPath(process.out)

process.schedule = cms.Schedule()
# process.schedule.append(process.unpack)
process.schedule.append(process.produce)
# process.schedule.append(process.finalProcess)


