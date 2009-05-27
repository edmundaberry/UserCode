#ifndef CALOTRIGTOWERMAP_H
#define CALOTRIGTOWERMAP_H

// System files
#include <memory>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// Geometry
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Framework
#include "FWCore/Framework/interface/ESHandle.h"

// Detector ids
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloIdCompiler.h"

using namespace std;
using namespace edm;

class CaloTrigTowerMap {
  
 public:

  //----------------------------------------------
  // Constructors, destructors, and operators
  //----------------------------------------------

  CaloTrigTowerMap();
  
  CaloTrigTowerMap(const CaloTrigTowerMap&);

  virtual ~CaloTrigTowerMap();

  CaloTrigTowerMap & operator=(const CaloTrigTowerMap&);

 public:

  //----------------------------------------------
  // User-accessible maps
  //----------------------------------------------

  void setGeometry( const CaloGeometry *geometry, 
		    const CaloTowerConstituentsMap *caloTowerConstituentsMap,
		    const EcalTrigTowerConstituentsMap *ecalTrigTowerConstituentsMap );
		    
  vector<CaloTowerDetId> getCaloTowers(HcalTrigTowerDetId hcalTrigTowerDetId);
  vector<CaloTowerDetId> getCaloTowers(EcalTrigTowerDetId ecalTrigTowerDetId);

 private:

  //----------------------------------------------
  // Calorimeter trigger geometries
  //----------------------------------------------

  const CaloGeometry                 *m_geometry;
  const CaloTowerConstituentsMap     *m_caloTowerConstituentsMap;
  const EcalTrigTowerConstituentsMap *m_ecalTrigTowerConstituentsMap;
  HcalTrigTowerGeometry               m_hcalTrigTowerGeometry;

};

#endif 
