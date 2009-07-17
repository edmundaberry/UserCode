import FWCore.ParameterSet.Config as cms

caloTowersFromTrigPrimsAnalyzer = cms.EDAnalyzer('CaloTowersFromTrigPrimsAnalyzer',
                                                 caloTowerTag = cms.untracked.InputTag("caloTowersFromTrigPrimsCreator"),
                                                 hcalTrigPrimTag = cms.untracked.InputTag("simHcalTriggerPrimitiveDigis"),
                                                 ecalTrigPrimTag = cms.untracked.InputTag("simEcalTriggerPrimitiveDigis"),
                                                 outputFileName = cms.untracked.string("CaloTowersFromTrigPrimsAnalyzer.root"),
                                                 verbose = cms.untracked.bool(False)
                                                 )
