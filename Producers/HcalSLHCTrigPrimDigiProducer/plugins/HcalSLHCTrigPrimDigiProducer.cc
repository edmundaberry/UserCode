#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiProducer.h"
#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiCollection.h"

//------------------------------------------------------
// Constructor, 
//   get tags for data collections
//   and instruction bools here.
//------------------------------------------------------

HcalSLHCTrigPrimDigiProducer::HcalSLHCTrigPrimDigiProducer(const edm::ParameterSet& iConfig):
  m_hcalSLHCTriggerPrimitiveAlgo ( new HcalSLHCTriggerPrimitiveAlgo())
{ produces<HcalSLHCTrigPrimDigiCollection>("HcalSLHCTrigPrimDigiCollection"); }

//------------------------------------------------------
// Destructor
//------------------------------------------------------

HcalSLHCTrigPrimDigiProducer::~HcalSLHCTrigPrimDigiProducer(){
  if (m_hcalSLHCTriggerPrimitiveAlgo != 0) delete m_hcalSLHCTriggerPrimitiveAlgo;
}

//------------------------------------------------------
// Main production function
//------------------------------------------------------

void HcalSLHCTrigPrimDigiProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
    
  //------------------------------------------------------
  // Create an empty collection
  //------------------------------------------------------

  std::auto_ptr<HcalSLHCTrigPrimDigiCollection> collection(new HcalSLHCTrigPrimDigiCollection());
  
  //------------------------------------------------------
  // Add the final CaloTowerCollection to the event
  //------------------------------------------------------
  
  iEvent.put(collection,"HcalSLHCTrigPrimDigiCollection");

}

void HcalSLHCTrigPrimDigiProducer::beginJob(){}

void HcalSLHCTrigPrimDigiProducer::endJob() {}

DEFINE_FWK_MODULE(HcalSLHCTrigPrimDigiProducer);
