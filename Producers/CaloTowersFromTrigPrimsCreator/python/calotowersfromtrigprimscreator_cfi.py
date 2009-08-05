import FWCore.ParameterSet.Config as cms

from Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi import *
from RecoJets.Configuration.CaloTowersES_cfi import *

caloTowersFromTrigPrimsCreator = cms.EDProducer('CaloTowersFromTrigPrimsCreator',                                                
  hadThreshold = cms.double(0.0),                                                
  emThreshold  = cms.double(0.0),                                                
  hcalTrigPrimTag = cms.InputTag("simHcalTriggerPrimitiveDigis"),
  ecalTrigPrimTag = cms.InputTag("simEcalTriggerPrimitiveDigis"),
  useHF   = cms.bool(True),                                               
  verbose = cms.bool(False)                                                
)
