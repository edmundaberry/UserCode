#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"

#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiProducer.h"
#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiCollection.h"

//------------------------------------------------------
// Constructor, 
//   get tags for data collections
//   and instruction bools here.
//------------------------------------------------------

HcalSLHCTrigPrimDigiProducer::HcalSLHCTrigPrimDigiProducer(const edm::ParameterSet& iConfig):
  m_hbheDigisTag ( iConfig.getParameter<edm::InputTag>("hbheDigis")),
  m_hfDigisTag   ( iConfig.getParameter<edm::InputTag>("hfDigis"  )),
  m_hcalSLHCTriggerPrimitiveAlgo ( new HcalSLHCTriggerPrimitiveAlgo(iConfig.getParameter<bool>("peakFilter"),
								    iConfig.getParameter<std::vector<double> >("weights"),
								    iConfig.getParameter<int>("latency"),
								    iConfig.getParameter<uint32_t>("FG_threshold"),
								    iConfig.getParameter<uint32_t>("ZS_threshold"),
								    iConfig.getParameter<int>("firstTPSample"),
								    iConfig.getParameter<int>("TPSize"),
								    iConfig.getParameter<int>("minIsoDepth"),
								    iConfig.getParameter<int>("maxIsoDepth")))				   
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
  // Get the edm::Handles from the Event
  //------------------------------------------------------

  edm::Handle<HBHEDigiCollection> hbheDigis;
  edm::Handle<HFDigiCollection>   hfDigis;

  bool gotHBHE = iEvent.getByLabel ( m_hbheDigisTag , hbheDigis );
  bool gotHF   = iEvent.getByLabel ( m_hfDigisTag   , hfDigis   );
  
  if (!gotHBHE) { 
    edm::LogWarning("HcalSLHCTrigPrimDigiProducer") << "Cannot get " << m_hbheDigisTag; 
    return; 
  }
  
  if (!gotHF  ) { 
    edm::LogWarning("HcalSLHCTrigPrimDigiProducer") << "Cannot get " << m_hfDigisTag  ; 
    return; 
  }

  //------------------------------------------------------
  // Get the coders from EventSetup
  //------------------------------------------------------

  edm::ESHandle<HcalTPGCoder> inputCoder;
  edm::ESHandle<CaloTPGTranscoder> outTranscoder;

  iSetup.get<HcalTPGRecord>().get(inputCoder);  
  iSetup.get<CaloTPGRecord>().get(outTranscoder);
  
  outTranscoder->setup(iSetup,CaloTPGTranscoder::HcalTPG);

  //------------------------------------------------------
  // Create an empty collection
  //------------------------------------------------------

  std::auto_ptr<HcalSLHCTrigPrimDigiCollection> result (new HcalSLHCTrigPrimDigiCollection());
    
  //------------------------------------------------------
  // Run the algorithm
  //------------------------------------------------------

  m_hcalSLHCTriggerPrimitiveAlgo -> run(inputCoder.product(),
					outTranscoder->getHcalCompressor().get(),
					*hbheDigis,  *hfDigis, *result);

  //------------------------------------------------------
  // Add the final CaloTowerCollection to the event
  //------------------------------------------------------
  
  iEvent.put(result,"HcalSLHCTrigPrimDigiCollection");

  outTranscoder->releaseSetup();

}

void HcalSLHCTrigPrimDigiProducer::beginJob(){}

void HcalSLHCTrigPrimDigiProducer::endJob() {}

DEFINE_FWK_MODULE(HcalSLHCTrigPrimDigiProducer);
