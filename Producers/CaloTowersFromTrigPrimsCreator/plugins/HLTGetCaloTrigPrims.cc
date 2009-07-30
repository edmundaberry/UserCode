#include <memory>

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HLTGetCaloTrigPrims : public edm::EDAnalyzer {
public:
  explicit HLTGetCaloTrigPrims(const edm::ParameterSet&);
  ~HLTGetCaloTrigPrims();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag m_hcalTrigPrimTag;
  edm::InputTag m_ecalTrigPrimTag;
  
};

HLTGetCaloTrigPrims::HLTGetCaloTrigPrims(const edm::ParameterSet& iConfig):
  m_hcalTrigPrimTag ( iConfig.getParameter<edm::InputTag>("hcalTrigPrimTag")),
  m_ecalTrigPrimTag ( iConfig.getParameter<edm::InputTag>("ecalTrigPrimTag"))
{}


HLTGetCaloTrigPrims::~HLTGetCaloTrigPrims(){}

void HLTGetCaloTrigPrims::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<HcalTrigPrimDigiCollection> HCALTrigPrimDigis;
  edm::Handle<EcalTrigPrimDigiCollection> ECALTrigPrimDigis;

  iEvent.getByLabel(m_hcalTrigPrimTag, HCALTrigPrimDigis );
  iEvent.getByLabel(m_ecalTrigPrimTag, ECALTrigPrimDigis );

  int i = 0; 
  i += (HCALTrigPrimDigis -> size());
  i += (ECALTrigPrimDigis -> size());

}

void HLTGetCaloTrigPrims::beginJob(const edm::EventSetup&) {}


void HLTGetCaloTrigPrims::endJob() {}

DEFINE_FWK_MODULE(HLTGetCaloTrigPrims);
