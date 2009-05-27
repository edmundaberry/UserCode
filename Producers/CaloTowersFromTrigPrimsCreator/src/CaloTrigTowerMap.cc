// System files
#include <iostream>
#include <map>
#include <utility>
#include <vector>

// Main header file
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTrigTowerMap.h"

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Framework
#include "FWCore/Framework/interface/ESHandle.h"

using namespace edm;
using namespace std;

CaloTrigTowerMap::CaloTrigTowerMap(){}

CaloTrigTowerMap::CaloTrigTowerMap(const CaloTrigTowerMap& copy){
  m_geometry                 = copy.m_geometry;
  m_caloTowerConstituentsMap = copy.m_caloTowerConstituentsMap;
}

CaloTrigTowerMap::~CaloTrigTowerMap(){}

void CaloTrigTowerMap::setGeometry( const CaloGeometry *geometry, 
					   const CaloTowerConstituentsMap *caloTowerConstituentsMap,
					   const EcalTrigTowerConstituentsMap *ecalTrigTowerConstituentsMap
					   ){
  m_caloTowerConstituentsMap     = caloTowerConstituentsMap;
  m_ecalTrigTowerConstituentsMap = ecalTrigTowerConstituentsMap;
  m_geometry                     = geometry;
}

CaloTrigTowerMap & CaloTrigTowerMap::operator=(const CaloTrigTowerMap &copy){

  m_geometry                     = copy.m_geometry;
  m_ecalTrigTowerConstituentsMap = copy.m_ecalTrigTowerConstituentsMap;
  m_caloTowerConstituentsMap     = copy.m_caloTowerConstituentsMap;
  

  return *this;
}

vector<CaloTowerDetId> CaloTrigTowerMap::getCaloTowers(HcalTrigTowerDetId hcalTrigTowerDetId){
  
  vector <CaloTowerDetId> retval;

  vector<HcalDetId> HcalDetIds = m_hcalTrigTowerGeometry.detIds(hcalTrigTowerDetId);
  
  vector<HcalDetId>::iterator hcalDetId = HcalDetIds.begin();
  
  for (; hcalDetId != HcalDetIds.end(); hcalDetId++){
    CaloTowerDetId caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(*hcalDetId);
    retval.push_back (caloTowerDetId);
  }

  sort   (retval.begin(), retval.end());  
  unique (retval.begin(), retval.end());
  
  return retval;
  
}

vector<CaloTowerDetId> CaloTrigTowerMap::getCaloTowers(EcalTrigTowerDetId ecalTrigTowerDetId){
  
  vector <CaloTowerDetId> retval;

  vector<DetId> DetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(ecalTrigTowerDetId);
  
  vector<DetId>::iterator detId = DetIds.begin();

  for (; detId != DetIds.end(); detId++){

    CaloTowerDetId caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(*detId);
    
    if (!caloTowerDetId.null()) retval.push_back(caloTowerDetId);
    
  }

  sort   (retval.begin(), retval.end());  
  unique (retval.begin(), retval.end());

  return retval;

}

