import FWCore.ParameterSet.Config as cms

from CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi import *
HcalTPGCoderULUT = cms.ESProducer("HcalTPGCoderULUT",
    read_Ascii_LUTs = cms.bool(True),
##    read_Ascii_LUTs = cms.bool(False),
##    read_XML_LUTs = cms.bool(False),
    read_XML_LUTs = cms.bool(True),
    LUTGenerationMode = cms.bool(True),
    DumpL1TriggerObjects = cms.bool(False),
    TagName = cms.string("LUTFromHCALDb"),
    AlgoName = cms.string("LUTFromHCALDb"),
##    inputLUTs = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/inputLUTcoder_SLHC.dat'),
    inputLUTs = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/inputLUTcoder_physics.dat'),
    filename = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/RecHit-TPG-calib.dat')
)

simHcalSLHCTriggerPrimitiveDigis = cms.EDProducer('HcalSLHCTrigPrimDigiProducer',
  hbheDigis = cms.InputTag("simHcalDigis"),
  hfDigis   = cms.InputTag("simHcalDigis"),
  latency = cms.int32(1),
  weights = cms.vdouble(1.0, 1.0), ##hardware algo 
  peakFilter = cms.bool(True),
  firstTPSample = cms.int32(2),
  TPSize = cms.int32(4),
  FG_threshold = cms.uint32(12),   ## threshold for setting fine grain bit
  minIsoDepth = cms.int32(0),
  maxIsoDepth = cms.int32(20),
  excludeDepth5 = cms.bool(True),
  CompressionLSB = cms.double(1.0), ## LSB for HBHEDataFrame ADC -> IntegerCaloSamples
  SLHCMode = cms.bool(True)         ## When converting from HBHEDataFrame ADC -> IntegerCaloSamples,
                                    ## use the normal coder (False) ? Or use the LSB (True) ?                                          
)
