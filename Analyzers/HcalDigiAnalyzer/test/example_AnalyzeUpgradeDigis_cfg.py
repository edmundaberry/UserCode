import FWCore.ParameterSet.Config as cms

process = cms.Process("MyHcalDigiAnalyzer")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")

##-----------------------------------------------
## Input number of events 
##-----------------------------------------------

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

##-----------------------------------------------
## Source goes here.  
##-----------------------------------------------

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),    
    fileNames = cms.untracked.vstring(
	'file:/tigress-hsm/tpestun/CMSSW_2_2_5/src/myfile.root'
    )
)

##-----------------------------------------------
## Unpack the digis
##-----------------------------------------------

process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.unpack = cms.Path(process.RawToDigi)

##-----------------------------------------------
## Declare analyzer objects for each DataFrame type
##-----------------------------------------------

process.hbDigiAnalyzer          = cms.EDAnalyzer("HBHEDigiAnalyzer", 
	hcalDigiTag             = cms.untracked.InputTag ("simHcalUnsuppressedDigis"),
	caloTowerTag            = cms.untracked.InputTag ("towerMaker"),
	hcalTrigPrimTag         = cms.untracked.InputTag ("simHcalTriggerPrimitiveDigis"),
	outPath                 = cms.untracked.string   ("/home/eberry/data/"),
	outSuffix               = cms.untracked.string   ("_testHardcode"),
	subdetName              = cms.untracked.string   ("HB"),
	useEventList            = cms.untracked.bool     (False)
)			        
				    
process.heDigiAnalyzer          = cms.EDAnalyzer("HBHEDigiAnalyzer",		
	hcalDigiTag             = cms.untracked.InputTag ("simHcalUnsuppressedDigis"),
	caloTowerTag            = cms.untracked.InputTag ("towerMaker"),
	hcalTrigPrimTag         = cms.untracked.InputTag ("simHcalTriggerPrimitiveDigis"),
	outPath                 = cms.untracked.string   ("/home/eberry/data/"),
	outSuffix               = cms.untracked.string   ("_testHardcode"),
	subdetName              = cms.untracked.string   ("HE"),
	useEventList            = cms.untracked.bool     (False)
)			        
			        
process.hoDigiAnalyzer          = cms.EDAnalyzer("HODigiAnalyzer",
	hcalDigiTag             = cms.untracked.InputTag ("simHcalUnsuppressedDigis"),
	caloTowerTag            = cms.untracked.InputTag ("towerMaker"),
	hcalTrigPrimTag         = cms.untracked.InputTag ("simHcalTriggerPrimitiveDigis"),
	outPath                 = cms.untracked.string   ("/home/eberry/data/"),
	outSuffix               = cms.untracked.string   ("_testHardcode"),
	subdetName              = cms.untracked.string   ("HO"),
	useEventList            = cms.untracked.bool     (False)
)			        
			        
process.hfDigiAnalyzer          = cms.EDAnalyzer("HFDigiAnalyzer",
	hcalDigiTag             = cms.untracked.InputTag ("simHcalUnsuppressedDigis"),
	caloTowerTag            = cms.untracked.InputTag ("towerMaker"),
	hcalTrigPrimTag         = cms.untracked.InputTag ("simHcalTriggerPrimitiveDigis"),
	outPath                 = cms.untracked.string   ("/home/eberry/data/"),
	outSuffix               = cms.untracked.string   ("_testHardcode"),
	subdetName              = cms.untracked.string   ("HF"),
        useEventList            = cms.untracked.bool     (False),
	doRecHitPulseCorrection = cms.untracked.bool     (False)
)


process.analyze = cms.Path(
		  process.hbDigiAnalyzer
		* process.heDigiAnalyzer
	        * process.hoDigiAnalyzer
		* process.hfDigiAnalyzer
		)

process.schedule = cms.Schedule()
process.schedule.append(process.analyze)

