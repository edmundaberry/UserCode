#ifndef CALOTRIGTOWERMAP_H
#define CALOTRIGTOWERMAP_H

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// Geometry
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

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
  
  void setGeometry( const CaloGeometry                 *geometry, 
		    const CaloTowerConstituentsMap     *caloTowerConstituentsMap,
		    const EcalTrigTowerConstituentsMap *ecalTrigTowerConstituentsMap );
  
  std::vector<CaloTowerDetId> getCaloTowers(HcalTrigTowerDetId hcalTrigTowerDetId);
  std::vector<CaloTowerDetId> getCaloTowers(EcalTrigTowerDetId ecalTrigTowerDetId);

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
