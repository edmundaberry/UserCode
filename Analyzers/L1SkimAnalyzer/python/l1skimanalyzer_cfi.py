import FWCore.ParameterSet.Config as cms

l1SkimAnalyzer = cms.EDAnalyzer('L1SkimAnalyzer',
                 
                                l1CaloEmCandsTag = cms.untracked.InputTag("gctDigis"),
                                l1GctJetCandsCenJetsTag = cms.untracked.InputTag("gctDigis","cenJets"),
                                l1GctJetCandsTauJetsTag = cms.untracked.InputTag("gctDigis","tauJets"),
                                l1GctJetCandsForJetsTag = cms.untracked.InputTag("gctDigis","forJets"),    
                                l1GctEtHadsTag = cms.untracked.InputTag("gctDigis"),
                                l1DecisionWordTag = cms.untracked.InputTag("newGtDigis"),
                                hltRecoCaloJetCandsTag = cms.untracked.InputTag("hltIterativeCone5CaloJets"),
                                outSuffix = cms.untracked.string("_L1EmulatorOnMC")

)
