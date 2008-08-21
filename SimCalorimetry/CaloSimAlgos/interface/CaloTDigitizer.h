#ifndef CaloSimAlgos_CaloTDigitizer_h
#define CaloSimAlgos_CaloTDigitizer_h

/** Turns hits into digis.  Assumes that 
    there's an ElectroncsSim class with the
    interface analogToDigital(const CaloSamples &, Digi &);

*/
#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVNoiseHitGenerator.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CalibFormats/HcalObjects/interface/HcalText2DetIdConverter.h"
#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalCalibDetId.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"


#include <cassert>

template<class Traits>
class CaloTDigitizer
{
public:
  /// these are the types that need to be defined in the Traits
  /// class.  The ElectronicsSim needs to have an interface
  /// that you'll see in the run() method
  typedef typename Traits::ElectronicsSim ElectronicsSim;
  typedef typename Traits::Digi Digi;
  typedef typename Traits::DigiCollection DigiCollection;

  CaloTDigitizer(CaloHitResponse * hitResponse, ElectronicsSim * electronicsSim, bool addNoise)
  :  theHitResponse(hitResponse),
     theNoiseHitGenerator(0),
     theElectronicsSim(electronicsSim),
     theDetIds(0),
     addNoise_(addNoise)
  {
  }


  /// doesn't delete the pointers passed in
  ~CaloTDigitizer() {}

  /// tell the digitizer which cells exist
  void setDetIds(const std::vector<DetId> & detIds) {theDetIds = detIds;}

  void setNoiseHitGenerator(CaloVNoiseHitGenerator * generator) 
  {
    theNoiseHitGenerator = generator;
  }

  /// turns hits into digis
  void run(MixCollection<PCaloHit> & input, DigiCollection & output) {
    assert(theDetIds.size() != 0);

    theHitResponse->run(input);

    if(theNoiseHitGenerator != 0) addNoiseHits();

    theElectronicsSim->newEvent();

    // reserve space for how many digis we expect
    int nDigisExpected = addNoise_ ? theDetIds.size() : theHitResponse->nSignals();
    output.reserve(nDigisExpected);

    // make a raw digi for evey cell
    for(std::vector<DetId>::const_iterator idItr = theDetIds.begin();
        idItr != theDetIds.end(); ++idItr)
    {
       Digi digi(*idItr);
       CaloSamples * analogSignal = theHitResponse->findSignal(*idItr);
       bool needToDeleteSignal = false;
       // don't bother digitizing if no signal and no noise
       if(analogSignal == 0 && addNoise_) {
         // I guess we need to make a blank signal for this cell.
         // Don't bother storing it anywhere.
         analogSignal = new CaloSamples(theHitResponse->makeBlankSignal(*idItr));
         needToDeleteSignal = true;
       }
       if(analogSignal != 0) { 
         theElectronicsSim->analogToDigital(*analogSignal , digi);
         output.push_back(digi);
         if(needToDeleteSignal) delete analogSignal;
      }
    }

    // free up some memory
    theHitResponse->clear();
  }

  //------------------------------------------------------
  // My added function
  //------------------------------------------------------

  void runCRUZET(MixCollection<PCaloHit> & input, DigiCollection & output){

    assert(theDetIds.size() != 0);

    theHitResponse->run(input); 

    if(theNoiseHitGenerator != 0) addNoiseHits();

    theElectronicsSim->newEvent();

    int subdetId;
    HcalSubdetector hcalsubdet;

    // reserve space for how many digis we expect
    int nDigisExpected = addNoise_ ? theDetIds.size() : theHitResponse->nSignals();
    output.reserve(nDigisExpected);

    // make a raw digi for evey cell
    for(std::vector<DetId>::const_iterator idItr = theDetIds.begin();
        idItr != theDetIds.end(); ++idItr){
      
      hcalsubdet = HcalSubdetector((*idItr).subdetId());
      subdetId = (int) hcalsubdet;
      
      Digi digi(*idItr);
      CaloSamples * analogSignal = theHitResponse->findSignal(*idItr);
      bool needToDeleteSignal = false;
       // don't bother digitizing if no signal and no noise
      if(analogSignal == 0 && addNoise_) {
	// I guess we need to make a blank signal for this cell.
	// Don't bother storing it anywhere.
	analogSignal = new CaloSamples(theHitResponse->makeBlankSignal(*idItr));
	needToDeleteSignal = true;
      }
      if(analogSignal != 0) {
	if (subdetId == 1 ||
	    subdetId == 3 ||
	    subdetId == 4   )
	  theElectronicsSim->analogToDigitalCRUZET(*analogSignal , digi);
	else theElectronicsSim->analogToDigital(*analogSignal , digi);
	output.push_back(digi);
	if(needToDeleteSignal) delete analogSignal;
      }
    }
    
    // free up some memory
    theHitResponse->clear();

  }

  void addNoiseHits()
  {
    std::vector<PCaloHit> noiseHits;
    theNoiseHitGenerator->getNoiseHits(noiseHits);
    for(std::vector<PCaloHit>::const_iterator hitItr = noiseHits.begin(),
        hitEnd = noiseHits.end(); hitItr != hitEnd; ++hitItr)
    {
      theHitResponse->add(*hitItr);
    }
  }


private:
  CaloHitResponse * theHitResponse;
  CaloVNoiseHitGenerator * theNoiseHitGenerator;
  ElectronicsSim * theElectronicsSim;
  std::vector<DetId> theDetIds;
  bool addNoise_;
};

#endif

