#ifndef CALOIDCOMPILER_H
#define CALOIDCOMPILER_H

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

using namespace std;
using namespace edm;

class CaloIdCompiler {

 public:
  CaloIdCompiler();
  ~CaloIdCompiler();

  //----------------------------------------------
  // User-accessible functions
  //----------------------------------------------

  void setGeometry( const CaloGeometry *geometry, const CaloTowerConstituentsMap *caloTowerConstituentsMap);

  vector<HcalTrigTowerDetId> getAllHcalTrigTowerDetIds();
  vector<EcalTrigTowerDetId> getAllEcalTrigTowerDetIds();
  vector<CaloTowerDetId>     getAllCaloTowerDetIds    ();

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
  // Calorimeter trigger geometries
  //----------------------------------------------

  HcalTrigTowerGeometry m_hcalTrigTowerGeometry;
  //EcalTrigTowerGeometry m_ecalTrigTowerGeometry;

  //----------------------------------------------
  // User-set data
  //----------------------------------------------

  ESHandle<CaloGeometry>           m_geometry;
  const CaloTowerConstituentsMap * m_caloTowerConstituentsMap;

  //----------------------------------------------
  // Internal mapping tools
  //----------------------------------------------

  typedef multimap<HcalTrigTowerDetId, int, less<HcalTrigTowerDetId> > HcalTrigTowerDetIdCountMap;
  typedef multimap<HcalTrigTowerDetId, int, less<HcalTrigTowerDetId> >::value_type HcalTrigTowerDetIdCount;
  HcalTrigTowerDetIdCountMap m_hcalTrigTowerDetIdCountMap;
  
  typedef multimap<CaloTowerDetId, int, less<CaloTowerDetId> > CaloTowerDetIdCountMap;
  typedef multimap<CaloTowerDetId, int, less<CaloTowerDetId> >::value_type CaloTowerDetIdCount;
  CaloTowerDetIdCountMap m_caloTowerDetIdCountMap;



};

#endif
