#------------------------------------------------------------
# import libraries
#------------------------------------------------------------

import FWCore.ParameterSet.Config as cms
import sys, os

#------------------------------------------------------------
# Define raw input with default response
#------------------------------------------------------------

def get_input_with_default(default, prompt=">>> "):
    result = raw_input(str(prompt) + ' ['+str(default)+']:\n')
    if result == "": result = default
    return default

#------------------------------------------------------------
# Let the user say which run to look at 
#------------------------------------------------------------

RUNNUMBER = str ( raw_input ( "Enter Runnumber:\n" ) )

#------------------------------------------------------------
# Let the user say where to store the data
# TODO: Give a default area
#------------------------------------------------------------

folder = get_input_with_default ( "/tmp/eberry/LocalRuns/" , "Where to store data?"  )

#------------------------------------------------------------
# If the data storage area does not exist, try to create it 
#------------------------------------------------------------

if ( not os.path.exists ( folder ) ) :
    folder_cmd = "mkdir -p "+folder
    os.system ( folder_cmd ) 

#------------------------------------------------------------
# If the file isn't already on local disk, move it there
#------------------------------------------------------------

file_name  = "USC_" + str(RUNNUMBER) + ".root"
input_file = folder+"/"+file_name

if ( not os.path.exists ( input_file ) ) :
    file_input_remote_location = "/bigspool/usc/"+file_name
    scp_cmd = "scp cmshcal03.cern.ch:///" + file_input_remote_location + " " + input_file
    print (scp_cmd) 
    os.system ( scp_cmd ) 

#------------------------------------------------------------
# If the file/folder still aren't accessible, bail.
#------------------------------------------------------------

if ( not os.path.exists ( input_file ) ) :
    print "Cannot copy information on run " + str(RUNNUMBER) + " from cmshcal03... Bailing."
    sys.exit()

#------------------------------------------------------------
# Now the CMSSW configuration...
#------------------------------------------------------------

process = cms.Process('UNPACK')

#------------------------------------------------------------
# Configure the source
#------------------------------------------------------------

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("HcalTBSource",
    fileNames= cms.untracked.vstring(
        "file://"+input_file
    ),
    streams = cms.untracked.vstring('HCAL_Trigger', 'HCAL_SlowData','HCAL_DCC700')
)

#------------------------------------------------------------
# Configure the unpacker -> convert local run format
# to CMSSW EDM format
#------------------------------------------------------------

process.tbunpacker = cms.EDFilter("HcalTBObjectUnpacker",
    HcalTriggerFED       = cms.untracked.int32(1),
    HcalSlowDataFED      = cms.untracked.int32(3),
    HcalTDCFED           = cms.untracked.int32(-1),
    HcalSourcePosFED     = cms.untracked.int32(-1),
    IncludeUnmatchedHits = cms.untracked.bool(False)
)

#------------------------------------------------------------
# DB information... unclear why/if this is needed
#------------------------------------------------------------

process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)

#------------------------------------------------------------
# Load standard configuration files
#------------------------------------------------------------

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")
process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR10_P_V2::All'

#------------------------------------------------------------
# Use different tags, if needed/desired
#------------------------------------------------------------

#################### Load and Prefer New Resonse Corrections ###############

# process.es_pool = cms.ESSource("PoolDBESSource",
#                                process.CondDBSetup,
#                                timetype = cms.string('runnumber'),
#                                toGet = cms.VPSet(
#     cms.PSet(
#     record = cms.string("HcalGainsRcd"),
#     tag = cms.string("HcalGains_v2.08_hlt_TEST")
#     )),
#                                connect = cms.string('frontier://cmsfrontier.cern.ch:8000/FrontierProd/CMS_COND_31X_HCAL'),
#                                authenticationMethod = cms.untracked.uint32(0),
#                                )
#
# process.es_prefer_es_pool = cms.ESPrefer( "PoolDBESSource", "es_pool" )
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

#------------------------------------------------------------
# Unpack from RAW to HCAL digis
#------------------------------------------------------------

process.hcalDigis = cms.EDProducer("HcalRawToDigi",
    FilterDataQuality = cms.bool(True),
    HcalFirstFED = cms.untracked.int32(700),
    InputLabel = cms.InputTag("source"),
    UnpackCalib = cms.untracked.bool(True),
    FEDs = cms.untracked.vint32(
        700, 701, 702, 703, 704,
        705, 706, 707, 708, 709,
        710, 711, 712, 713, 714,
        715, 716, 717, 718, 719,
        720, 721, 722, 723, 724,
        725, 726, 727, 728, 729,
        730, 731, 732),
    streams = cms.untracked.vstring(
          'HCAL_Trigger','HCAL_SlowData','HCAL_QADCTDC'
    ),
    lastSample = cms.int32(9),
    firstSample = cms.int32(0),
    ComplainEmptyData = cms.untracked.bool(True)
)

#------------------------------------------------------------
# Configure the analyzer
#------------------------------------------------------------

process.analysis = cms.EDAnalyzer('HcalDigiAnalyzer',
                                  outPath   = cms.untracked.string ( folder + "/" ),
                                  outSuffix = cms.untracked.string ( "_" + str(RUNNUMBER ) )
)

#------------------------------------------------------------
# Define paths: DIGI, RECO, ANALYSIS
#------------------------------------------------------------

process.DIGI = cms.Path(process.hcalDigis ) 
process.RECO = cms.Path(process.hbheprereco  * process.horeco * process.hfreco )
process.ANALYSIS = cms.EndPath ( process.analysis ) 

#------------------------------------------------------------
# Schedule paths: DIGI, RECO, ANALYSIS
#------------------------------------------------------------

process.schedule = cms.Schedule()
process.schedule.append ( process.DIGI     ) 
process.schedule.append ( process.RECO     ) 
process.schedule.append ( process.ANALYSIS ) 

#------------------------------------------------------------
# DONE!
#------------------------------------------------------------
