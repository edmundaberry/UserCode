// Header file
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTrigTowerMap.h"

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
  
  std::cout << "    Setting CaloTowerConstituentsMap " << std::endl;
  
  m_caloTowerConstituentsMap     = caloTowerConstituentsMap;

  std::cout << "    Setting EcalTrigTowerConstituentsMap " << std::endl;

  m_ecalTrigTowerConstituentsMap = ecalTrigTowerConstituentsMap;

  std::cout << "    Setting CaloGeometry " << std::endl;

  m_geometry                     = geometry;
}

CaloTrigTowerMap & CaloTrigTowerMap::operator=(const CaloTrigTowerMap &copy){

  m_geometry                     = copy.m_geometry;
  m_ecalTrigTowerConstituentsMap = copy.m_ecalTrigTowerConstituentsMap;
  m_caloTowerConstituentsMap     = copy.m_caloTowerConstituentsMap;  

  return *this;
}

std::vector<CaloTowerDetId> CaloTrigTowerMap::getCaloTowers(HcalTrigTowerDetId hcalTrigTowerDetId){
  
  std::cout << "      In CaloTrigTowerMap::getCaloTowers " << std::endl;

  std::vector <CaloTowerDetId> retval;

  std::cout << "      Getting HcalDetId's for this trigger tower: "
	    << hcalTrigTowerDetId
	    << std::endl;

  std::vector<HcalDetId> HcalDetIds = m_hcalTrigTowerGeometry.detIds(hcalTrigTowerDetId);
  
  std::cout << "      Got HcalDetId's for this trigger tower: "
	    << hcalTrigTowerDetId
	    << std::endl;

  std::vector<HcalDetId>::iterator hcalDetId = HcalDetIds.begin();

  std::cout << "      Got THESE HcalDetId's for this trigger tower: " 
	    << hcalTrigTowerDetId
	    << std::endl;
  
  for (; hcalDetId != HcalDetIds.end(); hcalDetId++){
    std::cout << "        " << *hcalDetId 
	      << std::endl;
    CaloTowerDetId caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(*hcalDetId);
    retval.push_back (caloTowerDetId);
  }

  std::sort   (retval.begin(), retval.end());  
  std::unique (retval.begin(), retval.end());
  
  return retval;
  
}

std::vector<CaloTowerDetId> CaloTrigTowerMap::getCaloTowers(EcalTrigTowerDetId ecalTrigTowerDetId){
  
  std::vector <CaloTowerDetId> retval;

  std::vector<DetId> DetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(ecalTrigTowerDetId);
  
  std::vector<DetId>::iterator detId = DetIds.begin();

  for (; detId != DetIds.end(); detId++){

    CaloTowerDetId caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(*detId);
    
    if (!caloTowerDetId.null()) retval.push_back(caloTowerDetId);
    
  }

  std::sort   (retval.begin(), retval.end());  
  std::unique (retval.begin(), retval.end());

  return retval;

}

