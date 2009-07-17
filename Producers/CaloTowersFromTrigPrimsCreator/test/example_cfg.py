
import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

## Load standard sequences and event content
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")

## Configure the source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1))
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring (
   '/store/relval/CMSSW_2_2_9/RelValQCD_Pt_80_120/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_V12_v1/0001/EEB22623-1732-DE11-AC28-000423D944F8.root'
))

## Produce CaloTowers from trigger primitives
process.load("Producers.CaloTowersFromTrigPrimsCreator.calotowersfromtrigprimscreator_cfi")

## Reconstruct jets from new CaloTowers
process.load("RecoJets.JetProducers.iterativeCone5CaloJets_cff")
process.iterativeCone5TPGJets = process.iterativeCone5CaloJets.clone()
process.iterativeCone5TPGJets.src = "caloTowersFromTrigPrimsCreator"
process.iterativeCone5TPGJets.alias = "IC5TPGJet"

## Apply default jet corrections
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff")
process.prefer("L2L3JetCorrectorIC5Calo") 
process.L2L3CorJetIC5TPG = process.L2L3CorJetIC5Calo.clone()
process.L2L3CorJetIC5TPG.src = "iterativeCone5TPGJets"

## Relevant modules to keep
process.out = cms.OutputModule("PoolOutputModule",                              
   outputCommands = cms.untracked.vstring(
        'drop *',
        'keep EcalTriggerPrimitiveDigisSorted_simEcalTriggerPrimitiveDigis_*_*',
        'keep HcalTriggerPrimitiveDigisSorted_simHcalTriggerPrimitiveDigis_*_*',
        'keep CaloTowersSorted_caloTowersFromTrigPrimsCreator_*_*',
        'keep recoCaloJets_L2L3CorJetIC5TPG_*_*',
        'keep recoCaloJets_iterativeCone5TPGJets_*_*'
        ),
   fileName = cms.untracked.string('output.root')
)

## Declare paths
process.TOWERS   = cms.Path    ( process.caloTowersFromTrigPrimsCreator )
process.JETRECO  = cms.Path    ( process.iterativeCone5TPGJets )
process.JETCORR  = cms.Path    ( process.L2L3CorJetIC5TPG )
process.OUTPUT   = cms.EndPath ( process.out )

## Schedule paths
process.schedule = cms.Schedule()
process.schedule.append(process.TOWERS  )
process.schedule.append(process.JETRECO )
process.schedule.append(process.JETCORR )
process.schedule.append(process.OUTPUT  ) 

