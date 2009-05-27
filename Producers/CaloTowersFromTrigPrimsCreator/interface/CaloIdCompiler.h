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
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
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

  void setGeometry( const CaloGeometry *geometry, 
		    const CaloTowerConstituentsMap *caloTowerConstituentsMap,
		    const EcalTrigTowerConstituentsMap *ecalTrigTowerConstituentsMap  );

  vector<HcalTrigTowerDetId> getAllHcalTrigTowerDetIds();
  vector<EcalTrigTowerDetId> getAllEcalTrigTowerDetIds();
  vector<CaloTowerDetId>     getAllCaloTowerDetIds    ();

 private:

  //----------------------------------------------
  // Internal functions
  //----------------------------------------------

  template < class EcalDetId > 
  void getEcalTrigTowerDetIds_fromSubdet(vector<DetId> cellIds, vector<EcalTrigTowerDetId> &towerIds);  
  void getHcalTrigTowerDetIds_fromSubdet(vector<DetId> cellIds, vector<HcalTrigTowerDetId> &towerIds);
  void getCaloTowerDetIds_fromSubdet    (vector<DetId> cellIds, vector<CaloTowerDetId>     &towerIds);
  
  //----------------------------------------------
  // HcalDetId's
  //----------------------------------------------

  vector<DetId> m_hbCells;
  vector<DetId> m_heCells;
  vector<DetId> m_hoCells;
  vector<DetId> m_hfCells;

  vector<DetId> m_ebCells;
  vector<DetId> m_eeCells;
  vector<DetId> m_esCells;

  //----------------------------------------------
  // Calorimeter trigger geometries
  //----------------------------------------------

  HcalTrigTowerGeometry m_hcalTrigTowerGeometry;
  //EcalTrigTowerGeometry m_ecalTrigTowerGeometry;

  //----------------------------------------------
  // User-set data
  //----------------------------------------------

  const CaloGeometry                 * m_geometry;
  const EcalTrigTowerConstituentsMap * m_ecalTrigTowerConstituentsMap;
  const CaloTowerConstituentsMap     * m_caloTowerConstituentsMap;

  //----------------------------------------------
  // Internal mapping tools
  //----------------------------------------------

  typedef multimap<HcalTrigTowerDetId, int, less<HcalTrigTowerDetId> > HcalTrigTowerDetIdCountMap;
  typedef multimap<HcalTrigTowerDetId, int, less<HcalTrigTowerDetId> >::value_type HcalTrigTowerDetIdCount;
  HcalTrigTowerDetIdCountMap m_hcalTrigTowerDetIdCountMap;

  typedef multimap<EcalTrigTowerDetId, int, less<EcalTrigTowerDetId> > EcalTrigTowerDetIdCountMap;
  typedef multimap<EcalTrigTowerDetId, int, less<EcalTrigTowerDetId> >::value_type EcalTrigTowerDetIdCount;
  EcalTrigTowerDetIdCountMap m_ecalTrigTowerDetIdCountMap;
  
  typedef multimap<CaloTowerDetId, int, less<CaloTowerDetId> > CaloTowerDetIdCountMap;
  typedef multimap<CaloTowerDetId, int, less<CaloTowerDetId> >::value_type CaloTowerDetIdCount;
  CaloTowerDetIdCountMap m_caloTowerDetIdCountMap;



};

#endif
