#ifndef PRODUCERS_CALOTOWERSFROMTRIGPRIMSCREATOR_H
#define PRODUCERS_CALOTOWERSFROMTRIGPRIMSCREATOR_H

//------------------------------------------------------
// Include files
//------------------------------------------------------

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Algorithm
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAlgo.h"

class CaloTowersFromTrigPrimsCreator : public edm::EDProducer {
public:
  explicit CaloTowersFromTrigPrimsCreator(const edm::ParameterSet&);
  ~CaloTowersFromTrigPrimsCreator();
  
private:

  //------------------------------------------------------
  // Analysis functions
  //------------------------------------------------------

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
    
  //------------------------------------------------------
  // Instruction bools
  //------------------------------------------------------
    
  bool m_useHF;
  bool m_verbose;

  //------------------------------------------------------
  // CaloTower thresholds
  //------------------------------------------------------

  double m_hadThreshold;
  double m_emThreshold ;

  //------------------------------------------------------
  // Calorimeter depth information
  //------------------------------------------------------

  double m_momHBDepth;
  double m_momHEDepth;
  double m_momEBDepth;
  double m_momEEDepth;

  //------------------------------------------------------
  // InputTags 
  //------------------------------------------------------
  
  edm::InputTag m_hcalTrigPrimTag;
  edm::InputTag m_ecalTrigPrimTag;
  edm::InputTag m_defaultCaloTowerTag;

  //------------------------------------------------------
  // Geometry objects
  //------------------------------------------------------
  
  edm::ESHandle<EcalElectronicsMapping>       m_ecalElectronicsMapping;
  edm::ESHandle<CaloGeometry>                 m_geometry;
  edm::ESHandle<CaloSubdetectorGeometry>      m_ecalBarrelGeometry;
  edm::ESHandle<CaloSubdetectorGeometry>      m_ecalEndcapGeometry;
  HcalTrigTowerGeometry                       m_trigTowerGeometry;

  //------------------------------------------------------
  // Constituent mappings
  //------------------------------------------------------
  
  edm::ESHandle<EcalTrigTowerConstituentsMap> m_ecalTrigTowerConstituentsMap;
  edm::ESHandle<CaloTowerConstituentsMap>     m_caloTowerConstituentsMap;

  //------------------------------------------------------
  // Energy scales
  //------------------------------------------------------

  edm::ESHandle<L1CaloEcalScale>   m_l1CaloEcalScale;
  edm::ESHandle<L1CaloHcalScale>   m_l1CaloHcalScale;

  //------------------------------------------------------
  // Algorithm
  //------------------------------------------------------

  CaloTowersFromTrigPrimsAlgo m_caloTowersFromTrigPrimsAlgo;

};

#endif
