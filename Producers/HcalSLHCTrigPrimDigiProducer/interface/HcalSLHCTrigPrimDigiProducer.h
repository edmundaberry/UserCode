#ifndef PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLHCTRIGPRIMDIGIPRODUCER_H
#define PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLHCTRIGPRIMDIGIPRODUCER_H

//------------------------------------------------------
// Include files
//------------------------------------------------------

#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Algorithm
#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTriggerPrimitiveAlgo.h"

class HcalSLHCTrigPrimDigiProducer : public edm::EDProducer {
public:
  explicit HcalSLHCTrigPrimDigiProducer(const edm::ParameterSet&);
  ~HcalSLHCTrigPrimDigiProducer();
  
private:

  //------------------------------------------------------
  // Analaysis functions
  //------------------------------------------------------

  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;  

  //------------------------------------------------------
  // InputTags
  //------------------------------------------------------
  
  

  //------------------------------------------------------
  // Algorithm
  //------------------------------------------------------
 
  HcalSLHCTriggerPrimitiveAlgo * m_hcalSLHCTriggerPrimitiveAlgo;
  
};

#endif
