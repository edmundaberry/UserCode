import FWCore.ParameterSet.Config as cms

l1SkimAnalyzerData = cms.EDAnalyzer('L1SkimAnalyzer',
                                
                                l1CaloEmCandsTag = cms.untracked.InputTag("simRctDigis"),
                                l1GctJetCandsCenJetsTag = cms.untracked.InputTag("simGctDigis","cenJets"),
                                l1GctJetCandsTauJetsTag = cms.untracked.InputTag("simGctDigis","tauJets"),
                                l1GctJetCandsForJetsTag = cms.untracked.InputTag("simGctDigis","forJets"),    
                                l1GctEtHadsTag = cms.untracked.InputTag("simGctDigis"),
                                l1DecisionWordTag = cms.untracked.InputTag("simGtDigis"),
                                hltRecoCaloJetCandsTag = cms.untracked.InputTag("hltIterativeCone5CaloJets"),
                                outSuffix = cms.untracked.string("_L1EmulatorOnData")
                                
)
