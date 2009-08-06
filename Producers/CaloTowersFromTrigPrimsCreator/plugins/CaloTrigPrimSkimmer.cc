#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTrigPrimSkimmer.h"

//-----------------------------------------------------
// Constructor
//-----------------------------------------------------

CaloTrigPrimSkimmer::CaloTrigPrimSkimmer(const edm::ParameterSet& iConfig):
  m_hcalTrigPrimTag( iConfig.getParameter<edm::InputTag>("hcalTrigPrimTag" )),
  m_ecalTrigPrimTag( iConfig.getParameter<edm::InputTag>("ecalTrigPrimTag" )) {
	       
  produces<HcalTrigPrimDigiCollection>("");
  produces<EcalTrigPrimDigiCollection>("");
  
}

//-----------------------------------------------------
// Destructor
//-----------------------------------------------------

CaloTrigPrimSkimmer::~CaloTrigPrimSkimmer(){}

//-----------------------------------------------------
// Event work function
//-----------------------------------------------------

void CaloTrigPrimSkimmer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
 
  //-----------------------------------------------------
  // Declare output 
  //-----------------------------------------------------

  std::auto_ptr<HcalTrigPrimDigiCollection> SkimmedHcalTrigPrimDigis (new HcalTrigPrimDigiCollection());
  std::auto_ptr<EcalTrigPrimDigiCollection> SkimmedEcalTrigPrimDigis (new EcalTrigPrimDigiCollection());
  
  //-----------------------------------------------------
  // Get trigger primitives from the event
  //-----------------------------------------------------
  
  edm::Handle<HcalTrigPrimDigiCollection> InitialHcalTrigPrimDigis;
  edm::Handle<EcalTrigPrimDigiCollection> InitialEcalTrigPrimDigis;
  
  iEvent.getByLabel(m_hcalTrigPrimTag, InitialHcalTrigPrimDigis );
  iEvent.getByLabel(m_ecalTrigPrimTag, InitialEcalTrigPrimDigis );
  
  //-----------------------------------------------------
  // Perform the skimming
  //-----------------------------------------------------

  skim ( *InitialHcalTrigPrimDigis, *SkimmedHcalTrigPrimDigis );
  skim ( *InitialEcalTrigPrimDigis, *SkimmedEcalTrigPrimDigis );
  
  //-----------------------------------------------------
  // Put the skimmed collections into the event
  //-----------------------------------------------------
  
  iEvent.put(SkimmedHcalTrigPrimDigis);
  iEvent.put(SkimmedEcalTrigPrimDigis);
  
}

//-----------------------------------------------------
// Hcal skimming function
//-----------------------------------------------------

void CaloTrigPrimSkimmer::skim ( const HcalTrigPrimDigiCollection & InitialHcalCol,
				 HcalTrigPrimDigiCollection & SkimmedHcalCol ){

  HcalTrigPrimDigiCollection::const_iterator initial_hcal_tp     = InitialHcalCol.begin();
  HcalTrigPrimDigiCollection::const_iterator initial_hcal_tp_end = InitialHcalCol.end();

  for (; initial_hcal_tp != initial_hcal_tp_end; ++initial_hcal_tp ){
    
    if ( initial_hcal_tp -> SOI_fineGrain()    != 0 ||
	 initial_hcal_tp -> SOI_compressedEt() != 0 )  SkimmedHcalCol.push_back (*initial_hcal_tp);
    
    
  }
  
}

//-----------------------------------------------------
// Ecal skimming function
//-----------------------------------------------------

void CaloTrigPrimSkimmer::skim ( const EcalTrigPrimDigiCollection & InitialEcalCol,
				 EcalTrigPrimDigiCollection & SkimmedEcalCol ){
  
  EcalTrigPrimDigiCollection::const_iterator initial_ecal_tp     = InitialEcalCol.begin();
  EcalTrigPrimDigiCollection::const_iterator initial_ecal_tp_end = InitialEcalCol.end();
  
  for (; initial_ecal_tp != initial_ecal_tp_end; ++initial_ecal_tp ){
    
    if ( initial_ecal_tp -> fineGrain()    != 0 ||
  	 initial_ecal_tp -> compressedEt() != 0 )  SkimmedEcalCol.push_back ( *initial_ecal_tp);
  
  }

}

//-----------------------------------------------------
// Begin/end job routines
//-----------------------------------------------------

void CaloTrigPrimSkimmer::beginJob(const edm::EventSetup&){}

void CaloTrigPrimSkimmer::endJob() {}

//-----------------------------------------------------
// Define the module
//-----------------------------------------------------

DEFINE_FWK_MODULE(CaloTrigPrimSkimmer);
