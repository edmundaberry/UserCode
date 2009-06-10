import FWCore.ParameterSet.Config as cms

from Geometry.CaloEventSetup.CaloTowerConstituents_cfi import *
from Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi import *
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import *
from CalibCalorimetry.EcalTPGTools.ecalTPGScale_cff import *
from RecoJets.Configuration.CaloTowersES_cfi import *

caloTowersFromTrigPrimsCreator = cms.EDProducer('CaloTowersFromTrigPrimsCreator',
                                                
  momHBDepth = cms.untracked.double(0.2),
  momHEDepth = cms.untracked.double(0.4),
  momEBDepth = cms.untracked.double(0.3),
  momEEDepth = cms.untracked.double(0.0),   
  hcalTrigPrimTag = cms.untracked.InputTag("hcalDigis"),
  ecalTrigPrimTag = cms.untracked.InputTag("ecalDigis","EcalTriggerPrimitives"),
  defaultCaloTowersTag = cms.untracked.InputTag("towerMaker"),
  verbose = cms.untracked.bool(False)
                                                
)
