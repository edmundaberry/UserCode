// -*- C++ -*-
//
// Package:    SRLutTest
// Class:      SRLutTest
// 
/**\class SRLutTest SRLutTest.cc Analyzers/SRLutTest/src/SRLutTest.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Edmund A BERRY
//         Created:  Thu Sep 17 12:05:48 CEST 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class decleration
//

class SRLutTest : public edm::EDAnalyzer {
   public:
  explicit SRLutTest(const edm::ParameterSet&);
  ~SRLutTest();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SRLutTest::SRLutTest(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
}


SRLutTest::~SRLutTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SRLutTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   
   edm::Handle< CSCCorrelatedLCTDigiCollection > CorrLCTs;
   
   iEvent.getByLabel (edm::InputTag("simCscTriggerPrimitiveDigis","MPCSORTED"), CorrLCTs);
   
   
   CSCCorrelatedLCTDigiCollection::DigiRangeIterator chamber     = CorrLCTs.product() -> begin();
   CSCCorrelatedLCTDigiCollection::DigiRangeIterator chamber_end = CorrLCTs.product() -> end();
   
   for (; chamber != chamber_end; ++chamber){
     
     CSCDetId * cscDetId = &((*chamber).first);
     
     int station      = cscDetId -> station      ();
     int chamber      = cscDetId -> chamber      ();
     int triggerSector= cscDetId -> triggerSector();
     
     bool interesting_chamber = (chamber > 9 && station == 4 && triggerSector == 2);
     
     if (interesting_chamber) {
       std::cout << "CSCDetId is : " << *cscDetId << std::endl;
       std::cout << "  chamber = " << chamber << std::endl;       
     }
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
SRLutTest::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SRLutTest::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SRLutTest);
