import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")
generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(140),
        MinPt = cms.double(10),
        PartID = cms.vint32(13),
        MaxEta = cms.double(1.8),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(0.9),
        MinPhi = cms.double(-3.14159265359) ## in radians

    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('single mu pt 100'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
