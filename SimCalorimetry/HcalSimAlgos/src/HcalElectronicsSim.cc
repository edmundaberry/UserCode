#include "SimCalorimetry/HcalSimAlgos/interface/HcalElectronicsSim.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalAmplifier.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalCoderFactory.h"
#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"
#include "DataFormats/HcalDigi/interface/HFDataFrame.h"
#include "DataFormats/HcalDigi/interface/ZDCDataFrame.h"
#include "CLHEP/Random/RandFlat.h"

HcalElectronicsSim::HcalElectronicsSim(HcalAmplifier * amplifier, const HcalCoderFactory * coderFactory)
  : theAmplifier(amplifier),
    theCoderFactory(coderFactory),
    theRandFlat(0),
    theStartingCapId(0)
{}


HcalElectronicsSim::~HcalElectronicsSim()
{
  delete theRandFlat;
}


void HcalElectronicsSim::setRandomEngine(CLHEP::HepRandomEngine & engine)
{
  theRandFlat = new CLHEP::RandFlat(engine);
}


template<class Digi> 
void HcalElectronicsSim::convert(CaloSamples & frame, Digi & result) {
  
  result.setSize(frame.size());
  theAmplifier->amplify(frame);
  theCoderFactory->coder(frame.id())->fC2adc(frame, result, theStartingCapId);
}

//------------------------------------------------------
// My convert function

template<class Digi> 
void HcalElectronicsSim::convertCRUZET(CaloSamples & frame, Digi & result){
  
  result.setSize(frame.size());
  theAmplifier->amplifyCRUZET(frame);  
  theCoderFactory->coder(frame.id())->fC2adc(frame, result, theStartingCapId);
}

//------------------------------------------------------

void HcalElectronicsSim::analogToDigital(CaloSamples & lf, HBHEDataFrame & result) {
  convert<HBHEDataFrame>(lf, result);
}


void HcalElectronicsSim::analogToDigital(CaloSamples & lf, HODataFrame & result) {
  convert<HODataFrame>(lf, result);
}


void HcalElectronicsSim::analogToDigital(CaloSamples & lf, HFDataFrame & result) {
  convert<HFDataFrame>(lf, result);
}

void HcalElectronicsSim::analogToDigital(CaloSamples & lf, ZDCDataFrame & result) {
  convert<ZDCDataFrame>(lf, result);
}

//------------------------------------------------------
// My analogToDigital functions

void HcalElectronicsSim::analogToDigitalCRUZET(CaloSamples & lf, HBHEDataFrame & result) {
  convertCRUZET<HBHEDataFrame>(lf, result);
}


void HcalElectronicsSim::analogToDigitalCRUZET(CaloSamples & lf, HODataFrame & result) {
  convertCRUZET<HODataFrame>(lf, result);
}


void HcalElectronicsSim::analogToDigitalCRUZET(CaloSamples & lf, HFDataFrame & result) {
  convertCRUZET<HFDataFrame>(lf, result);
}

void HcalElectronicsSim::analogToDigitalCRUZET(CaloSamples & lf, ZDCDataFrame & result) {
  convertCRUZET<ZDCDataFrame>(lf, result);
}

//------------------------------------------------------


void HcalElectronicsSim::newEvent() {
  // pick a new starting Capacitor ID
  theStartingCapId = theRandFlat->fireInt(4);
  theAmplifier->setStartingCapId(theStartingCapId);
}

