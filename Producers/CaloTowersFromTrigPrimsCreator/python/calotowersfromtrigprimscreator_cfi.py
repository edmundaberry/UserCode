import FWCore.ParameterSet.Config as cms

from Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi import *
from RecoJets.Configuration.CaloTowersES_cfi import *

caloTowersFromTrigPrimsCreator = cms.EDProducer('CaloTowersFromTrigPrimsCreator',
                                                
  momHBDepth = cms.untracked.double(0.2),
  momHEDepth = cms.untracked.double(0.4),
  momEBDepth = cms.untracked.double(0.3),
  momEEDepth = cms.untracked.double(0.0),   
  hadThreshold = cms.untracked.double(0.0),                                                
  emThreshold  = cms.untracked.double(0.0),                                                
  hcalTrigPrimTag = cms.untracked.InputTag("simHcalTriggerPrimitiveDigis"),
  ecalTrigPrimTag = cms.untracked.InputTag("simEcalTriggerPrimitiveDigis"),
  defaultCaloTowersTag = cms.untracked.InputTag("towerMaker"),
  useHF = cms.untracked.bool(True),                                               
  verbose = cms.untracked.bool(False)
                                                
)
