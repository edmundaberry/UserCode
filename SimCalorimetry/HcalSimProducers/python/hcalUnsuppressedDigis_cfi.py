import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HcalSimProducers.hcalSimParameters_cfi import *
simHcalUnsuppressedDigisMyOverlay = cms.EDProducer("HcalDigiProducer",
    hcalSimParameters,
    doNoise = cms.bool(True),
    doHPDNoise = cms.bool(False),
    doTimeSlew = cms.bool(True),
    doHFWindow = cms.bool(True),
    hitsProducer = cms.string('g4SimHits'),

    takeNoiseFromCRUZETData = cms.bool(True),
    startingEventNumber = cms.int32(0),
                                                    
    hbFile = cms.string("/uscms/home/eberry/data/HcalFinalOutput_HB_data.root"),
    heFile = cms.string("/uscms/home/eberry/data/HcalFinalOutput_HE_data.root"),
    hoFile = cms.string("/uscms/home/eberry/data/HcalFinalOutput_HO_data.root"),
    hfFile = cms.string("/uscms/home/eberry/data/HcalFinalOutput_HF_data.root")

)



