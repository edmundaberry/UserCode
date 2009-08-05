#ifndef PRODUCERS_CALOTOWERSFROMTRIGPRIMSCREATOR_CALOTOWERSFROMTRIGPRIMSCREATOR_H
#define PRODUCERS_CALOTOWERSFROMTRIGPRIMSCREATOR_CALOTOWERSFROMTRIGPRIMSCREATOR_H

//------------------------------------------------------
// Include files
//------------------------------------------------------

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

// Algorithm
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAlgo.h"

class CaloTowersFromTrigPrimsCreator : public edm::EDProducer {
public:
  explicit CaloTowersFromTrigPrimsCreator(const edm::ParameterSet&);
  ~CaloTowersFromTrigPrimsCreator();
  
private:

  //------------------------------------------------------
  // Analysis functions
  //------------------------------------------------------

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  //------------------------------------------------------
  // InputTags 
  //------------------------------------------------------
  
  edm::InputTag m_hcalTrigPrimTag;
  edm::InputTag m_ecalTrigPrimTag;
  
  //------------------------------------------------------
  // Algorithm
  //------------------------------------------------------

  CaloTowersFromTrigPrimsAlgo * m_caloTowersFromTrigPrimsAlgo;

};

#endif
