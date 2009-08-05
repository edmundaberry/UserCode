import FWCore.ParameterSet.Config as cms

process = cms.Process("HLT2")

# Message logging
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

# Load source
process.load("Producers.CaloTowersFromTrigPrimsCreator.RelValQCD_LowPt_cfi")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

# Load geometry/detector info
process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Timing service
process.PathTimerService = cms.Service( "PathTimerService" )
process.hltTimer = cms.EDProducer( "PathTimerInserter" )

# Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V12::All'

# Setup HLT for running my mini menu
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Producers.CaloTowersFromTrigPrimsCreator.MyMiniHLTMenu_cff")
process.schedule = process.MyHLTSchedule 

# Setup output
process.load("Configuration.EventContent.EventContent_cff")

## process.hltPoolOutput = cms.OutputModule("PoolOutputModule",
##                                          outputCommands = cms.untracked.vstring( 'drop *', 'keep HLTPerformanceInfo_*_*_*' ),
##                                          fileName = cms.untracked.string('file:HLT.root')
## )

process.HLTOutput = cms.EndPath(process.hltTimer)

# Add output to the schedule
process.schedule.append(process.HLTOutput)
