#ifndef CALOTRIGPRIMSKIMMER_H
#define CALOTRIGPRIMSKIMMER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

class CaloTrigPrimSkimmer : public edm::EDProducer {
public:
  explicit CaloTrigPrimSkimmer(const edm::ParameterSet&);
  ~CaloTrigPrimSkimmer();
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void skim (const HcalTrigPrimDigiCollection & InitialHcalCol,
	     HcalTrigPrimDigiCollection & SkimmedHcalCol);
  
  void skim (const EcalTrigPrimDigiCollection & InitialEcalCol,
	     EcalTrigPrimDigiCollection & SkimmedEcalCol);
  
  edm::InputTag m_hcalTrigPrimTag;
  edm::InputTag m_ecalTrigPrimTag;
      
};

#endif
