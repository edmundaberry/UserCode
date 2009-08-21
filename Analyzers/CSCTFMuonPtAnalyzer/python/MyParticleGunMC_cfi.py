import FWCore.ParameterSet.Config as cms

source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
        "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job100_100events0-140GeV.root",
        "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job101_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job102_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job103_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job104_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job105_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job106_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job107_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job108_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job109_100events0-140GeV.root",
         "rfio:/castor/cern.ch/user/e/eberry/MC_SAMPLES/PARTICLE_GUN/CSC_MUON_GUN/CMSSW_3_2_4/CSCMuonGun_STARTUP31X_V4_GEN-SIM-DIGI-L1_job10_100events0-140GeV.root"
),
duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

