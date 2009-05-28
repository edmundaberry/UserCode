#include "Tools/CaloIdCompiler/interface/CaloIdCompiler.h"

CaloIdCompiler::CaloIdCompiler(){}

CaloIdCompiler::~CaloIdCompiler(){}

void CaloIdCompiler::setGeometry( const CaloGeometry *geometry, 
				  const CaloTowerConstituentsMap *caloTowerConstituentsMap,
				  const EcalTrigTowerConstituentsMap * ecalTrigTowerConstituentsMap ){

  m_ecalTrigTowerConstituentsMap = ecalTrigTowerConstituentsMap;
  m_caloTowerConstituentsMap = caloTowerConstituentsMap;
  m_geometry = geometry;

  m_hbCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalBarrel );
  m_heCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalEndcap );
  m_hoCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalOuter  );
  m_hoCells = m_geometry -> getValidDetIds(DetId::Hcal, HcalForward);  

  m_ebCells = m_geometry -> getValidDetIds(DetId::Ecal, EcalBarrel);
  m_eeCells = m_geometry -> getValidDetIds(DetId::Ecal, EcalEndcap);
  m_esCells = m_geometry -> getValidDetIds(DetId::Ecal, EcalPreshower);


}

void CaloIdCompiler::getHcalTrigTowerDetIds_fromSubdet(std::vector<DetId> cellIds, 
						       std::vector<HcalTrigTowerDetId> & towerIds){
  
  int entries;
  
  std::vector<DetId>::iterator cellId_iter = cellIds.begin();
  
  for (; cellId_iter != cellIds.end(); ++cellId_iter){
    HcalDetId hcalDetId = HcalDetId(*cellId_iter);
    std::vector<HcalTrigTowerDetId> temp = m_hcalTrigTowerGeometry.towerIds(*cellId_iter);
    std::vector<HcalTrigTowerDetId>::iterator temp_iter = temp.begin();
    for (; temp_iter != temp.end(); ++temp_iter){
      
      entries = m_hcalTrigTowerDetIdCountMap.count(*temp_iter);
      
      if (entries == 0) towerIds.push_back(*temp_iter);
      m_hcalTrigTowerDetIdCountMap.insert(HcalTrigTowerDetIdCount(*temp_iter,entries+1));
	
    }
  }
}

template < class EcalDetId > 
void CaloIdCompiler::getEcalTrigTowerDetIds_fromSubdet (std::vector<DetId> cellIds, 
							std::vector<EcalTrigTowerDetId> & towerIds){
  
  
  
  int entries;
  std::vector<DetId>::iterator cellId_iter = cellIds.begin();
  for (; cellId_iter != cellIds.end(); ++cellId_iter){
    
    EcalDetId ecalDetId = EcalDetId(*cellId_iter);
    EcalTrigTowerDetId temp = m_ecalTrigTowerConstituentsMap -> towerOf(ecalDetId);
  
    entries = m_ecalTrigTowerDetIdCountMap.count(temp);

    if (entries == 0) towerIds.push_back(temp);

    m_ecalTrigTowerDetIdCountMap.insert(EcalTrigTowerDetIdCount(temp,entries+1));
  }
}


void CaloIdCompiler::getCaloTowerDetIds_fromSubdet(std::vector<DetId> cellIds, 
								 std::vector<CaloTowerDetId> & towerIds){

  int entries;

  std::vector<DetId>::iterator cellId_iter = cellIds.begin();
  
  for (; cellId_iter != cellIds.end(); ++cellId_iter){

    HcalDetId hcalDetId = HcalDetId(*cellId_iter);
    CaloTowerDetId temp = m_caloTowerConstituentsMap -> towerOf(*cellId_iter);

    entries = m_caloTowerDetIdCountMap.count(temp);

    if (entries == 0) towerIds.push_back(temp);
    
    m_caloTowerDetIdCountMap.insert(CaloTowerDetIdCount(temp,entries+1));
  }
}

std::vector<HcalTrigTowerDetId> CaloIdCompiler::getAllHcalTrigTowerDetIds(){
  
  std::vector <HcalTrigTowerDetId> retval;
  
  getHcalTrigTowerDetIds_fromSubdet(m_hbCells, retval);
  getHcalTrigTowerDetIds_fromSubdet(m_heCells, retval);
  getHcalTrigTowerDetIds_fromSubdet(m_hoCells, retval);
  getHcalTrigTowerDetIds_fromSubdet(m_hfCells, retval);
  
  std::sort   (retval.begin(), retval.end());  
  std::unique (retval.begin(), retval.end());
  
  return retval;
}

std::vector<EcalTrigTowerDetId> CaloIdCompiler::getAllEcalTrigTowerDetIds(){

  std::vector <EcalTrigTowerDetId> retval;

  getEcalTrigTowerDetIds_fromSubdet<EBDetId>(m_ebCells, retval);
  getEcalTrigTowerDetIds_fromSubdet<EEDetId>(m_eeCells, retval);
  getEcalTrigTowerDetIds_fromSubdet<ESDetId>(m_esCells, retval);
  
  return retval;
}

std::vector<CaloTowerDetId>     CaloIdCompiler::getAllCaloTowerDetIds(){

  std::vector <CaloTowerDetId> retval;

  getCaloTowerDetIds_fromSubdet(m_hbCells, retval);
  getCaloTowerDetIds_fromSubdet(m_heCells, retval);
  getCaloTowerDetIds_fromSubdet(m_hoCells, retval);
  getCaloTowerDetIds_fromSubdet(m_hfCells, retval);

  std::vector<CaloTowerDetId>::iterator iter = retval.begin();

  std::sort   (retval.begin(), retval.end());  
  std::unique (retval.begin(), retval.end());
  
  return retval;
}
