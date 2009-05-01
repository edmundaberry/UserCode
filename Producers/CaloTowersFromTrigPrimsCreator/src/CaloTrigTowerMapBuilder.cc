// System files
#include <iostream>
#include <map>
#include <utility>
#include <vector>

// Main header file
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTrigTowerMapBuilder.h"

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

CaloTrigTowerMapBuilder::CaloTrigTowerMapBuilder(){}

CaloTrigTowerMapBuilder::CaloTrigTowerMapBuilder(const CaloTrigTowerMapBuilder& copy){
  m_geometry                 = copy.m_geometry;
  m_caloTowerConstituentsMap = copy.m_caloTowerConstituentsMap;
  
  m_hbCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalBarrel );
  m_heCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalEndcap );
  m_hoCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalOuter  );
  m_hoCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalForward);
}

CaloTrigTowerMapBuilder::~CaloTrigTowerMapBuilder(){}

void CaloTrigTowerMapBuilder::setGeometry( const CaloGeometry *geometry, const CaloTowerConstituentsMap *caloTowerConstituentsMap){

  m_caloTowerConstituentsMap = caloTowerConstituentsMap;
  m_geometry = geometry;

  m_hbCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalBarrel );
  m_heCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalEndcap );
  m_hoCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalOuter  );
  m_hoCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalForward);

}

CaloTrigTowerMapBuilder & CaloTrigTowerMapBuilder::operator=(const CaloTrigTowerMapBuilder &copy){

  m_geometry                 = copy.m_geometry;
  m_caloTowerConstituentsMap = copy.m_caloTowerConstituentsMap;

  m_hbCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalBarrel );
  m_heCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalEndcap );
  m_hoCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalOuter  );
  m_hoCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalForward);

  return *this;

}

void CaloTrigTowerMapBuilder::getEcalTrigTowerDetIds_fromSubdet(std::vector<DetId> cellIds, 
								 std::vector<EcalTrigTowerDetId> & towerIds){

  cout << "Not implemented yet!" << endl;

}

void CaloTrigTowerMapBuilder::buildMap(){

  caloIdCompiler.setGeometry(m_geometry.product(), m_caloTowerConstituentsMap);

  vector<HcalTrigTowerDetId> hcalTrigTowerDetIds = caloIdCompiler.getAllHcalTrigTowerDetIds();
  vector<CaloTowerDetId>     caloTowerDetIds     = caloIdCompiler.getAllCaloTowerDetIds    ();

  vector<HcalTrigTowerDetId>::iterator hcalTrigTowerDetId_iter = hcalTrigTowerDetIds.begin();

  for(; hcalTrigTowerDetId_iter != hcalTrigTowerDetIds.end(); ++hcalTrigTowerDetId_iter){
    
    vector<HcalDetId> HcalDetIds = m_hcalTrigTowerGeometry.detIds(*hcalTrigTowerDetId_iter);
    
    vector<HcalDetId>::iterator hcalDetId = HcalDetIds.begin();
    
    for (; hcalDetId != HcalDetIds.end(); hcalDetId++){
      
      CaloTowerDetId caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(*hcalDetId);
      
      if (caloToTrigTowerMap.count(caloTowerDetId) == 0) {
	
	caloToTrigTowerMap.insert(CaloToTrigTowerMap_Item(caloTowerDetId, *hcalTrigTowerDetId_iter));
	trigToCaloTowerMap.insert(TrigToCaloTowerMap_Item(*hcalTrigTowerDetId_iter, caloTowerDetId));
	
      }
      
      else {
	
	CaloToTrigTowerMap::const_iterator caloToTrigTowerMap_iter = caloToTrigTowerMap.begin();
	
	bool pairAlreadyMapped = false;
	
	for (; caloToTrigTowerMap_iter !=  caloToTrigTowerMap.end(); caloToTrigTowerMap_iter++){
	  pairAlreadyMapped |= (caloToTrigTowerMap_iter -> first  ==  caloTowerDetId &&
				caloToTrigTowerMap_iter -> second == *hcalTrigTowerDetId_iter );
	}
	
	if (!pairAlreadyMapped){
	  caloToTrigTowerMap.insert(CaloToTrigTowerMap_Item(caloTowerDetId, *hcalTrigTowerDetId_iter));
	  trigToCaloTowerMap.insert(TrigToCaloTowerMap_Item(*hcalTrigTowerDetId_iter, caloTowerDetId));
	}
	
      }
    }
  }
  
}
