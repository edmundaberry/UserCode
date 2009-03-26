#ifndef HcalSimAlgos_HcalAmplifier_h
#define HcalSimAlgos_HcalAmplifier_h

#include "CalibFormats/HcalObjects/interface/HcalCalibrationWidths.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVSimParameterMap.h"
#include "CLHEP/Random/RandGaussQ.h"

class HcalDbService;

class HcalAmplifier {
public:
  HcalAmplifier(const CaloVSimParameterMap * parameters, bool addNoise);
  virtual ~HcalAmplifier(){ delete theRandGaussQ; }
  
  // the Producer will probably update this every event
  void setDbService(const HcalDbService * service) {
    theDbService = service;
  }
  
  //------------------------------------------------------
  // My added functions:

  void setFileNames (std::string hbFile, std::string heFile, std::string hofile, std::string hfFile);
  void setAmplifierEvent(int iEvent);
  void setId(HcalDetId id);
  void readTree();
  virtual void amplifyCRUZET(CaloSamples & linearFrame) const;

  //------------------------------------------------------

  void setRandomEngine(CLHEP::HepRandomEngine & engine);
  
  virtual void amplify(CaloSamples & linearFrame) const;

  void makeNoise (const HcalCalibrationWidths& width, int fFrames, double* fGauss, double* fNoise) const;


  void setStartingCapId(int capId) {theStartingCapId = capId;}

private:

  const HcalDbService * theDbService;
  CLHEP::RandGaussQ * theRandGaussQ;
  const CaloVSimParameterMap * theParameterMap;

  unsigned theStartingCapId;
  bool addNoise_;
  
  //------------------------------------------------------
  // My added parameters

  double noiseArray_[5][200][100][5][10];
  double haveDataForThisChannel_[5][83][73][4];
  bool takeNoiseFromCRUZETData_;
  int amplifierEvent_;
  std::string hbFile_, heFile_, hoFile_, hfFile_;

  //------------------------------------------------------x
  
};

#endif
