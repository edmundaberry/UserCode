#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiComparer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

HcalSLHCTrigPrimDigiComparer::HcalSLHCTrigPrimDigiComparer(const edm::ParameterSet& iConfig){}

HcalSLHCTrigPrimDigiComparer::~HcalSLHCTrigPrimDigiComparer(){}


void HcalSLHCTrigPrimDigiComparer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

}


void HcalSLHCTrigPrimDigiComparer::beginJob(){}

void HcalSLHCTrigPrimDigiComparer::endJob() {}

DEFINE_FWK_MODULE(HcalSLHCTrigPrimDigiComparer);
