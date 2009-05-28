#ifndef CALOIDCOMPILER_H
#define CALOIDCOMPILER_H

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// Geometry
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

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

  std::vector<HcalTrigTowerDetId> getAllHcalTrigTowerDetIds();
  std::vector<EcalTrigTowerDetId> getAllEcalTrigTowerDetIds();
  std::vector<CaloTowerDetId>     getAllCaloTowerDetIds    ();

 private:

  //----------------------------------------------
  // Internal functions
  //----------------------------------------------

  template < class EcalDetId > 
  void getEcalTrigTowerDetIds_fromSubdet(std::vector<DetId> cellIds, std::vector<EcalTrigTowerDetId> &towerIds);  
  void getHcalTrigTowerDetIds_fromSubdet(std::vector<DetId> cellIds, std::vector<HcalTrigTowerDetId> &towerIds);
  void getCaloTowerDetIds_fromSubdet    (std::vector<DetId> cellIds, std::vector<CaloTowerDetId>     &towerIds);
  
  //----------------------------------------------
  // HcalDetId's
  //----------------------------------------------

  std::vector<DetId> m_hbCells;
  std::vector<DetId> m_heCells;
  std::vector<DetId> m_hoCells;
  std::vector<DetId> m_hfCells;

  std::vector<DetId> m_ebCells;
  std::vector<DetId> m_eeCells;
  std::vector<DetId> m_esCells;

  //----------------------------------------------
  // Calorimeter trigger geometries
  //----------------------------------------------

  HcalTrigTowerGeometry m_hcalTrigTowerGeometry;

  //----------------------------------------------
  // User-set data
  //----------------------------------------------

  const CaloGeometry                 * m_geometry;
  const EcalTrigTowerConstituentsMap * m_ecalTrigTowerConstituentsMap;
  const CaloTowerConstituentsMap     * m_caloTowerConstituentsMap;

  //----------------------------------------------
  // Internal mapping tools
  //----------------------------------------------

  typedef std::multimap<HcalTrigTowerDetId, int, std::less<HcalTrigTowerDetId> > HcalTrigTowerDetIdCountMap;
  typedef std::multimap<HcalTrigTowerDetId, int, std::less<HcalTrigTowerDetId> >::value_type HcalTrigTowerDetIdCount;
  HcalTrigTowerDetIdCountMap m_hcalTrigTowerDetIdCountMap;

  typedef std::multimap<EcalTrigTowerDetId, int, std::less<EcalTrigTowerDetId> > EcalTrigTowerDetIdCountMap;
  typedef std::multimap<EcalTrigTowerDetId, int, std::less<EcalTrigTowerDetId> >::value_type EcalTrigTowerDetIdCount;
  EcalTrigTowerDetIdCountMap m_ecalTrigTowerDetIdCountMap;
  
  typedef std::multimap<CaloTowerDetId, int, std::less<CaloTowerDetId> > CaloTowerDetIdCountMap;
  typedef std::multimap<CaloTowerDetId, int, std::less<CaloTowerDetId> >::value_type CaloTowerDetIdCount;
  CaloTowerDetIdCountMap m_caloTowerDetIdCountMap;



};

#endif
