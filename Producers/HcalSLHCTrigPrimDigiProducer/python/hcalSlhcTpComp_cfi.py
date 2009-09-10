import FWCore.ParameterSet.Config as cms

hcalSLHCTrigPrimDigiComparer = cms.EDAnalyzer('HcalSLHCTrigPrimDigiComparer',
  defaultTrigPrimTag = cms.InputTag("simHcalTriggerPrimitiveDigis"),
  upgradeTrigPrimTag = cms.InputTag("simHcalSLHCTriggerPrimitiveDigis")
)


