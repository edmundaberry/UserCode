#ifndef CALOTRIGTOWERMAPBUILDER_H
#define CALOTRIGTOWERMAPBUILDER_H

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
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Framework
#include "FWCore/Framework/interface/ESHandle.h"

// Detector ids
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloIdCompiler.h"

using namespace std;
using namespace edm;

class CaloTrigTowerMapBuilder {
  
 public:

  //----------------------------------------------
  // Constructors, destructors, and operators
  //----------------------------------------------

  CaloTrigTowerMapBuilder();
  
  CaloTrigTowerMapBuilder(const CaloTrigTowerMapBuilder&);

  virtual ~CaloTrigTowerMapBuilder();

  CaloTrigTowerMapBuilder & operator=(const CaloTrigTowerMapBuilder&);

 public:

  //----------------------------------------------
  // User-accessible maps
  //----------------------------------------------

  typedef multimap<CaloTowerDetId, HcalTrigTowerDetId , less< CaloTowerDetId > > CaloToTrigTowerMap;
  typedef multimap<CaloTowerDetId, HcalTrigTowerDetId , less< CaloTowerDetId > >::value_type CaloToTrigTowerMap_Item;
  CaloToTrigTowerMap caloToTrigTowerMap;
    
  typedef multimap<HcalTrigTowerDetId , CaloTowerDetId, less< HcalTrigTowerDetId > > TrigToCaloTowerMap;
  typedef multimap<HcalTrigTowerDetId , CaloTowerDetId, less< HcalTrigTowerDetId > >::value_type TrigToCaloTowerMap_Item;
  TrigToCaloTowerMap trigToCaloTowerMap;
  
  void buildMap();  
  void setGeometry( const CaloGeometry *geometry, const CaloTowerConstituentsMap *caloTowerConstituentsMap);

  CaloToTrigTowerMap getCaloToTrigTowerMap() {return caloToTrigTowerMap;}
  TrigToCaloTowerMap getTrigToCaloTowerMap() {return trigToCaloTowerMap;}

 private:

  //----------------------------------------------
  // Internal functions
  //----------------------------------------------

  void getHcalTrigTowerDetIds_fromSubdet(std::vector<DetId> cellIds, std::vector<HcalTrigTowerDetId> & towerIds);
  void getEcalTrigTowerDetIds_fromSubdet(std::vector<DetId> cellIds, std::vector<EcalTrigTowerDetId> & towerIds);
  void getCaloTowerDetIds_fromSubdet    (std::vector<DetId> cellIds, std::vector<CaloTowerDetId>     & towerIds);
  
  //----------------------------------------------
  // HcalDetId's
  //----------------------------------------------

  vector<DetId> m_hbCells;
  vector<DetId> m_heCells;
  vector<DetId> m_hoCells;
  vector<DetId> m_hfCells;

  //----------------------------------------------
  // Mapping object for other det ids
  //----------------------------------------------

  CaloIdCompiler caloIdCompiler;

  //----------------------------------------------
  // Calorimeter trigger geometries
  //----------------------------------------------

  HcalTrigTowerGeometry m_hcalTrigTowerGeometry;
  //EcalTrigTowerGeometry m_ecalTrigTowerGeometry;

  //----------------------------------------------
  // User-set data
  //----------------------------------------------

  ESHandle<CaloGeometry>     m_geometry;
  const CaloTowerConstituentsMap * m_caloTowerConstituentsMap;
  
};

#endif 
