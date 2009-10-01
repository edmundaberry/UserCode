import FWCore.ParameterSet.Config as cms

source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
        'file:///uscms_data/d2/eberry/CSCMuonGun/GEN-SIM-DIGI-L1/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_jobs501-750_0-140GeV.root',
        'file:///uscms_data/d2/eberry/CSCMuonGun/GEN-SIM-DIGI-L1/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_jobs751-1000_0-140GeV.root',
        'file:///uscms_data/d2/eberry/CSCMuonGun/GEN-SIM-DIGI-L1/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_jobs1001-1250_0-140GeV.root',
        'file:///uscms_data/d2/eberry/CSCMuonGun/GEN-SIM-DIGI-L1/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_jobs1251-1500_0-140GeV.root',
),
duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

