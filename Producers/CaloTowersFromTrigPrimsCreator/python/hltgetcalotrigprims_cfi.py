import FWCore.ParameterSet.Config as cms

hltGetCaloTrigPrims = cms.EDAnalyzer('HLTGetCaloTrigPrims',
                                     hcalTrigPrimTag = cms.InputTag("simHcalTriggerPrimitiveDigis"),
                                     ecalTrigPrimTag = cms.InputTag("simEcalTriggerPrimitiveDigis")
)
