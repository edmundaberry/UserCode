import FWCore.ParameterSet.Config as cms


from CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi import *
HcalTPGCoderULUT = cms.ESProducer("HcalTPGCoderULUT",
    read_Ascii_LUTs = cms.bool(False),
    read_XML_LUTs = cms.bool(False),
    LUTGenerationMode = cms.bool(True),
    DumpL1TriggerObjects = cms.bool(False),
    TagName = cms.string("LUTFromHCALDb"),
    AlgoName = cms.string("LUTFromHCALDb"),
    inputLUTs = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/inputLUTcoder_physics.dat'),
    filename = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/RecHit-TPG-calib.dat')
)

simHcalSLHCTriggerPrimitiveDigis = cms.EDProducer('HcalSLHCTrigPrimDigiProducer',
  hbheDigis = cms.InputTag("simHcalDigis"),
  hfDigis   = cms.InputTag("simHcalDigis")                                               
)
