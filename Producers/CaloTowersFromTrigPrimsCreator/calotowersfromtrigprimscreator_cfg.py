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
        'file:/uscms/home/eberry/data/test.root'
    )
)

process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
process.load("CalibCalorimetry.EcalTPGTools.ecalTPGScale_cff")

process.unpack = cms.Path(process.RawToDigi)

process.myProducerLabel = cms.EDProducer('CaloTowersFromTrigPrimsCreator')

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('full_output.root'),
    outputCommands = cms.untracked.vstring('keep *')
)

process.produce = cms.Path(process.myProducerLabel)
process.finalProcess = cms.EndPath(process.out)

process.schedule = cms.Schedule()
##process.schedule.append(process.unpack)
process.schedule.append(process.produce)
# process.schedule.append(process.finalProcess)


