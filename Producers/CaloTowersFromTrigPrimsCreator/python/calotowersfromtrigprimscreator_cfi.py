import FWCore.ParameterSet.Config as cms

from Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi import *
from RecoJets.Configuration.CaloTowersES_cfi import *

caloTowersFromTrigPrimsCreator = cms.EDProducer('CaloTowersFromTrigPrimsCreator',                                                
  hadThreshold = cms.untracked.double(0.0),                                                
  emThreshold  = cms.untracked.double(0.0),                                                
  hcalTrigPrimTag = cms.untracked.InputTag("simHcalTriggerPrimitiveDigis"),
  ecalTrigPrimTag = cms.untracked.InputTag("simEcalTriggerPrimitiveDigis"),
  useHF = cms.untracked.bool(True),                                               
  verbose = cms.untracked.bool(False)                                                
)
