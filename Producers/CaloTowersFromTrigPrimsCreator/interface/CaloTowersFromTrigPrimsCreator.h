#ifndef PRODUCERS_CALOTOWERSFROMTRIGPRIMSCREATOR_H
#define PRODUCERS_CALOTOWERSFROMTRIGPRIMSCREATOR_H

//------------------------------------------------------
// Include files
//------------------------------------------------------

// Framework files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// Geometry
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Energy scale info
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloEcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

// Algorithm
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAlgo.h"

// class CaloTowerConstituentsMap;

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
  edm::InputTag m_hoTrigPrimTag;
  edm::InputTag m_ecalTrigPrimTag;
  edm::InputTag m_defaultCaloTowersTag;
  edm::InputTag m_ebDigiTag;
  edm::InputTag m_eeDigiTag;

  //------------------------------------------------------
  // Geometry objects
  //------------------------------------------------------

  edm::ESHandle<CaloGeometry>                 m_geometry;
  edm::ESHandle<CaloSubdetectorGeometry>      m_ecalBarrelGeometry;
  edm::ESHandle<CaloSubdetectorGeometry>      m_ecalEndcapGeometry;
  HcalTrigTowerGeometry                  m_trigTowerGeometry;

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
