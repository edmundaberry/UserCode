#ifndef PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLHCTRIGPRIMDIGICOMPARER_H
#define PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLHCTRIGPRIMDIGICOMPARER_H

//------------------------------------------------------
// Include files
//------------------------------------------------------

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HcalSLHCTrigPrimDigiComparer : public edm::EDAnalyzer {
   public:
      explicit HcalSLHCTrigPrimDigiComparer(const edm::ParameterSet&);
      ~HcalSLHCTrigPrimDigiComparer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      edm::InputTag m_upgradeDigisTag;
      edm::InputTag m_defaultDigisTag;

};

#endif
