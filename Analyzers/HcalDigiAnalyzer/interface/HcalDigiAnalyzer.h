#ifndef Analyzers_HcalDigiAnalyzer_HcalDigiAnalyzer_h
#define Analyzers_HcalDigiAnalyzer_HcalDigiAnalyzer_h

// System includes
#include <string> 

// HCAL includes
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

// Framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// My includes
#include "Analyzers/HcalDigiAnalyzer/interface/DigiTree.h"
#include "Analyzers/HcalDigiAnalyzer/interface/FillDigiTree.h"

class HcalDigiAnalyzer : public edm::EDAnalyzer {
public:
  explicit HcalDigiAnalyzer(const edm::ParameterSet&);
  ~HcalDigiAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  template <class DigiCollection, class RecoCollection> 
  void analyzeDigiCollection ( const HcalDbService & conditions,  const DigiCollection & digis , const RecoCollection & recos );
  
  //-----------------------------------------------------
  // Root tree objects
  //-----------------------------------------------------

  DigiTree     m_digiTree;
  FillDigiTree m_fillDigi;

  //-----------------------------------------------------
  // Out file info
  //-----------------------------------------------------

  std::string m_outPath;
  std::string m_outSuffix;
  std::string m_rootFile;  



};

#endif
