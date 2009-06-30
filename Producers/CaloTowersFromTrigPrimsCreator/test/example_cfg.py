
import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.EventContent.EventContent_cff")

process.GlobalTag.globaltag = 'IDEAL_V12::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'dcache:///pnfs/cms/WAX/resilient/eberry/CMSSW_ROOT_Files/ParticleGun/CentralElectronGun/RAW-DIGI-RECO/CentralElectronGun_IDEAL_V12_RAW-DIGI-RECO_job1_20events_100GeV-Focused.root'
))

process.load("Producers.CaloTowersFromTrigPrimsCreator.calotowersfromtrigprimscreator_cfi")

process.PRODUCE = cms.Path( process.caloTowersFromTrigPrimsCreator  )

process.schedule = cms.Schedule()
process.schedule.append(process.PRODUCE)


