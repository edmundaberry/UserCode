import FWCore.ParameterSet.Config as cms

caloTowerComparator = cms.EDAnalyzer('CaloTowerComp',
                                     defaultCaloTowerTag = cms.untracked.InputTag("towerMaker"),
                                     createdCaloTowerTag = cms.untracked.InputTag("caloTowersFromTrigPrimsCreator"),
                                     defaultCaloTowerFileName = cms.untracked.string("data/CaloTowerInfo_Default.root"),
                                     createdCaloTowerFileName = cms.untracked.string("data/CaloTowerInfo_Created.root")
                                     )
