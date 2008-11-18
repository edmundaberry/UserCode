#include "Analyzers/L1SkimAnalyzer/interface/L1RegionLookup.h"

void L1RegionLookup::SetArrays(){
  
  int temp_minIeta_for[12] = {29,32,35,38,-999,-999,-999,-999,29,32,35,38};
  int temp_maxIeta_for[12] = {31,34,37,41,-999,-999,-999,-999,31,34,37,41};
  
  int temp_minIeta_cen[15] = {1,5,9 ,13,17,21,25,-999,1,5,9 ,13,17,21,25};
  int temp_maxIeta_cen[15] = {4,8,12,16,20,24,28,-999,4,8,12,16,20,24,28};

  for (int i = 0; i < 12; i++){
    m_minHcalIetaAbs_for[i] = temp_minIeta_for[i];
    m_maxHcalIetaAbs_for[i] = temp_maxIeta_for[i];      
  }

  for (int i = 0; i < 15; i++){
    m_minHcalIetaAbs_cen[i] = temp_minIeta_cen[i];
    m_maxHcalIetaAbs_cen[i] = temp_maxIeta_cen[i];
  }
    
}

L1RegionLookup::L1RegionLookup(){
  
  SetArrays();
  
}

L1RegionLookup::~L1RegionLookup(){}

//-----------------------------------------------
// For each L1 Gct region, return an HCAL tower
// min ieta and max ieta
//-----------------------------------------------

int L1RegionLookup::getMinHcalIeta(unsigned int etaIndex, bool isForward){
  
  int min_ietaAbs, max_ietaAbs;
  int min_ieta;

  if (isForward){

    min_ietaAbs = m_minHcalIetaAbs_for[etaIndex];
    max_ietaAbs = m_maxHcalIetaAbs_for[etaIndex];
    
  }
  
  if (!isForward){

    min_ietaAbs = m_minHcalIetaAbs_cen[etaIndex];
    max_ietaAbs = m_maxHcalIetaAbs_cen[etaIndex];
    
  }

  if (etaIndex > 7) min_ieta = max_ietaAbs * (-1); // ieta < 0
  else if (etaIndex < 7) min_ieta = min_ietaAbs;
  else min_ieta = -999;
  
  return min_ieta;

}

int L1RegionLookup::getMaxHcalIeta(unsigned int etaIndex, bool isForward){
  
  int min_ietaAbs, max_ietaAbs;
  int max_ieta;

  if (isForward){

    min_ietaAbs = m_minHcalIetaAbs_for[etaIndex];
    max_ietaAbs = m_maxHcalIetaAbs_for[etaIndex];
    
  }
  
  if (!isForward){

    min_ietaAbs = m_minHcalIetaAbs_cen[etaIndex];
    max_ietaAbs = m_maxHcalIetaAbs_cen[etaIndex];
    
  }

  if (etaIndex > 7) max_ieta = min_ietaAbs * (-1); // ieta < 0
  else if (etaIndex < 7) max_ieta = max_ietaAbs;
  else max_ieta = -999;
  
  return max_ieta;

}
