import FWCore.ParameterSet.Config as cms

## Load the full HLT menu
from HLTrigger.Configuration.HLT_2E30_cff import *

## Load my CaloTower producer
from Producers.CaloTowersFromTrigPrimsCreator.calotowersfromtrigprimscreator_cfi import *

## Load my mini analyzer to get the HCAL/ECAL trigger primitives
from Producers.CaloTowersFromTrigPrimsCreator.hltgetcalotrigprims_cfi import *

## Reconstruct jets from my CaloTowers
hltIterativeCone5TPGJets = hltIterativeCone5CaloJets.clone()
hltIterativeCone5TPGJets.src = "caloTowersFromTrigPrimsCreator"
hltIterativeCone5TPGJets.alias = "IC5TPGJet"

## Apply default jet corrections to my jets
hltL2L3CorJetIC5TPG = hltMCJetCorJetIcone5.clone()
hltL2L3CorJetIC5TPG.src = "hltIterativeCone5TPGJets"

## Declare my jet reco sequence
HLTDoTPGJetRecoSequence = cms.Sequence( hltIterativeCone5TPGJets + hltL2L3CorJetIC5TPG )

## Declare the default HLT calo jet sequence
MyDoHLTCaloJetsSequence = cms.Sequence (HLTDoCaloSequence +             ## Calorimeter RECO sequence (includes CaloTowersCreator)
                                        HLTDoJetRecoSequence)           ## IterativeCone5 jet reconstruction and correction


## Declare my HLT TPG jet sequence
MyDoHLTTPGJetsSequence = cms.Sequence (
    ## hltGetCaloTrigPrims +
    ## caloTowersFromTrigPrimsCreator +
    caloTowersFromTrigPrimsCreator ## + ## My CaloTower producer
    ## HLTDoTPGJetRecoSequence
)         ## IterativeCone5 jet reconstruction and correction

MyHLTBeginSequence = cms.Sequence ( HLTBeginSequence + hltGetCaloTrigPrims )

## Declare trigger paths
## MyHLTriggerFirstPath   = cms.Path ( HLTBeginSequence + hltGetRaw + hltPreFirstPath + hltBoolFirstPath + HLTEndSequence )
MyHLTriggerTPGJetPath  = cms.Path ( HLTBeginSequence + MyDoHLTTPGJetsSequence  + HLTEndSequence )
## MyHLTriggerCaloJetPath = cms.Path ( HLTBeginSequence + MyDoHLTCaloJetsSequence + HLTEndSequence )

MyHLTSchedule = cms.Schedule()
## MyHLTSchedule.append( MyHLTriggerFirstPath   )
MyHLTSchedule.append( MyHLTriggerTPGJetPath  )
## MyHLTSchedule.append( MyHLTriggerCaloJetPath )


