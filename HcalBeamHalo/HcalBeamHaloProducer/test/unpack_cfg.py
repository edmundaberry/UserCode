import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("UNPACK")

process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring (
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_1.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_2.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_3.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_4.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_5.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_6.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_7.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_8.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_9.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_10.root',
        'rfio:/castor/cern.ch/user/r/rcr/ExpressPhysics/Skims/Run130445/ExpressPhysics__Commissioning10-Express-v3__FEVT____run130445_11.root'
      )
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration/StandardSequences/MagneticField_cff")
process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load("RecoMuon.Configuration.RecoMuonCosmics_cff")
process.GlobalTag.globaltag = 'GR10_P_V4::All'

process.outputSkim = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('analyze')),
    fileName = cms.untracked.string("/tmp/eberry/halo_with_digis.root")
    )

process.muonCSCDigis.UnpackStatusDigis = True

## process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32 (1000) )

process.analyze = cms.Path ( process.muonCSCDigis * process.STAmuontrackingforcosmicsEnsCapsOnly )
process.output  = cms.EndPath ( process.outputSkim )
