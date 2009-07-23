#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAlgo.h"
#include <limits> 

//----------------------------------------------------
// Constructor & Destructor
//----------------------------------------------------

CaloTowersFromTrigPrimsAlgo::CaloTowersFromTrigPrimsAlgo(double hadThreshold,
							 double emThreshold,
							 bool useHF,
							 bool verbose):
  m_hadThreshold(hadThreshold),
  m_emThreshold(emThreshold),
  m_useHF(useHF),
  m_verbose(verbose),
  MAX_ECALTT_IETA(28),
  MAX_ECALTT_IPHI(72),
  MAX_HCALTT_IETA(32),
  MAX_HCALTT_IPHI(72),
  MIN_HFTT_IETA(29),
  MAX_CT_IETA(CaloTowerDetId::kMaxIEta),
  MAX_CT_IPHI(CaloTowerDetId::kMaxIPhi),
  m_nDetIdsInCTFromEcal( MAX_CT_IETA * 2 + 1, std::vector < int > ( MAX_CT_IPHI + 1 , -999 )),
  m_nDetIdsInCTFromHcal( MAX_CT_IETA * 2 + 1, std::vector < int > ( MAX_CT_IPHI + 1 , -999 )),
  m_ecalTTMeanEta( MAX_ECALTT_IETA * 2 + 1, std::vector < double > ( MAX_ECALTT_IPHI + 1 , -999 )),
  m_hcalTTMeanEta( MAX_HCALTT_IETA * 2 + 1, std::vector < double > ( MAX_HCALTT_IPHI + 1 , -999 )),
  m_ecalTTMap( MAX_ECALTT_IETA * 2 + 1, std::vector< std::vector <CaloTowerDetId> >( MAX_ECALTT_IPHI + 1, std::vector<CaloTowerDetId>(0))),
  m_hcalTTMap( MAX_HCALTT_IETA * 2 + 1, std::vector< std::vector <CaloTowerDetId> >( MAX_HCALTT_IPHI + 1, std::vector<CaloTowerDetId>(0)))
  
{}

CaloTowersFromTrigPrimsAlgo::~CaloTowersFromTrigPrimsAlgo(){}

//----------------------------------------------------
// Loop through the TPG collections and assign energies
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::process(const EcalTrigPrimDigiCollection& EcalTrigPrimDigis){
					
  EcalTrigPrimDigiCollection::const_iterator trigPrimDigi = EcalTrigPrimDigis.begin();
  EcalTrigPrimDigiCollection::const_iterator trigPrimDigi_end = EcalTrigPrimDigis.end();

  double assignedEnergy = 0.0;    

  for(; trigPrimDigi != trigPrimDigi_end; ++trigPrimDigi) 
    assignedEnergy += assignEnergy(&(*trigPrimDigi));
  
  m_tpgEnergy += assignedEnergy;
  
}

void CaloTowersFromTrigPrimsAlgo::process(const HcalTrigPrimDigiCollection& HcalTrigPrimDigis){					
  
  HcalTrigPrimDigiCollection::const_iterator trigPrimDigi = HcalTrigPrimDigis.begin();
  HcalTrigPrimDigiCollection::const_iterator trigPrimDigi_end = HcalTrigPrimDigis.end();
  
  double assignedEnergy = 0.0;

  for(; trigPrimDigi != trigPrimDigi_end; ++trigPrimDigi)
    assignedEnergy += assignEnergy(&(*trigPrimDigi));

  m_tpgEnergy += assignedEnergy;
}

//----------------------------------------------------
// 1. Extract energy from trigger primitives, 
// 2. Divide energy among mapped CaloTowers
//----------------------------------------------------

template <typename TriggerPrimitiveDigi>
double CaloTowersFromTrigPrimsAlgo::assignEnergy(const TriggerPrimitiveDigi* trigPrimDigi){

  DetId::Detector det = (*trigPrimDigi).id().det();

  using std::numeric_limits; 

  //----------------------------------------------------
  // Get the energy from this trig prim
  //
  // -- If the total energy is zero, don't bother continuing
  // -- If the energy is infinitely big (map error), bomb.
  //----------------------------------------------------
  
  double totalEnergy = getTrigTowerEnergy(trigPrimDigi);
  if (totalEnergy == 0.0) return totalEnergy;
  
  assert ( totalEnergy != numeric_limits<double>::infinity() );
  assert ( totalEnergy == totalEnergy );
  
  //----------------------------------------------------
  // Get the trigger tower DetId & map it to CaloTowers
  //----------------------------------------------------
  
  typedef typename TriggerPrimitiveDigi::key_type TrigTowerDetId;
  
  TrigTowerDetId trigTowerDetId = trigPrimDigi -> id();

  std::vector<CaloTowerDetId> CaloTowerDetIds = getCaloTowers( trigTowerDetId );
  std::vector<CaloTowerDetId>::iterator caloTowerDetId = CaloTowerDetIds.begin();
  std::vector<CaloTowerDetId>::iterator caloTowerDetId_end = CaloTowerDetIds.end();

  double towerEnergy = totalEnergy / (double) CaloTowerDetIds.size();
    
  //----------------------------------------------------
  // Loop over the mapped CaloTowers and distribute
  // the energy to the appropriate MetaTowers
  //----------------------------------------------------

  for (; caloTowerDetId != caloTowerDetId_end; ++caloTowerDetId){
    
    //----------------------------------------------------
    // Find the MetaTower
    //----------------------------------------------------

    MetaTower &metaTower = find(*caloTowerDetId);
    
    //----------------------------------------------------
    // Now distribute the EM and hadronic energies.
    // E_outer is zeroed.
    //----------------------------------------------------
    
    metaTower.E_outer = 0;
    metaTower.E += towerEnergy;
    
    //----------------------------------------------------
    // First distribute hadronic energy
    //----------------------------------------------------
    
    if ( det == DetId::Hcal){     

      if (trigTowerDetId.ietaAbs() >= MIN_HFTT_IETA){	

	if ( m_useHF ){	
	  // Arbitrarily assume the HF is 50/50 EM vs Hadronic energy
	  metaTower.E_had += towerEnergy / 2.0;
	  metaTower.E_em  += towerEnergy / 2.0;
	}
	
      }
      
      else metaTower.E_had += towerEnergy;

    }
    
    //----------------------------------------------------
    // Now distribute the EM energy
    //----------------------------------------------------
    
    else if ( det == DetId::Ecal){
      metaTower.E_em  += towerEnergy;
    }
    
    else edm::LogWarning("CaloTowersFromTrigPrimsCreator") << "This trigger tower detector ID doesn't make sense -- " << trigTowerDetId << std::endl;
    
  }

  return totalEnergy;
} 

//----------------------------------------------------
// 1. Loop through the map of MetaTowers
// 2. Covert MetaTowers to CaloTowers
// 3. Push CaloTowers with energy and constituents
//    into the Collection
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::finish(CaloTowerCollection& result){					

  MetaTowerMap::const_iterator mapItr = m_metaTowerMap.begin();
  MetaTowerMap::const_iterator mapItr_end = m_metaTowerMap.end();

  for(; mapItr != mapItr_end; ++mapItr){

    // Don't bother with zero-energy MetaTowers
    double metaTowerEnergy = (mapItr->second).E_em + (mapItr->second).E_had;    
    if (metaTowerEnergy == 0.0) continue;
    
    // Convert to CaloTower
    CaloTower ct=convert(mapItr->first,mapItr->second);
    
    // Count this CaloTower's energy towards the events' total
    double caloTowerEnergy = ct.emEnergy() + ct.hadEnergy();
    m_cctEnergy += caloTowerEnergy;
    
    // Warn if you're going to lose a tower with energy because there are no constituents
    if (ct.constituentsSize() == 0 && caloTowerEnergy != 0.0)
      edm::LogError("CaloTowersFromTrigPrimsAlgo") << "This CaloTower: " << (mapItr->first) << " "
						   << "has no constituents and will be dropped. " 
						   << caloTowerEnergy << " GeV will be lost!";
    
    // Push the CaloTowers with energy and constituents into the collection
    if (ct.constituentsSize() > 0 && caloTowerEnergy > 0.0 ) 	
      result.push_back(ct);
  }   
  
  m_metaTowerMap.clear(); 
}

//----------------------------------------------------
// Convert MetaTower to CaloTower
//----------------------------------------------------

CaloTower CaloTowersFromTrigPrimsAlgo::convert(const CaloTowerDetId& id, const MetaTower& mt){

  //----------------------------------------------------
  // Transfer information from meta towers  
  //----------------------------------------------------

  double E=mt.E;
  double E_em=mt.E_em;
  double E_had=mt.E_had;
  double E_outer=mt.E_outer;
  
  //----------------------------------------------------  
  // Get location information
  //----------------------------------------------------
  
  GlobalPoint emPoint, hadPoint;  
  CaloTower::PolarLorentzVector towerP4;
  
  GlobalPoint p= m_caloTowerGeometry->getGeometry(id)->getPosition();
  double pf=1.0/cosh(p.eta());
  if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*pf, p.eta(), p.phi(), 0);  // simple momentum assignment, same position
  emPoint  = p;    
  hadPoint = p;  

  //----------------------------------------------------
  // Create the CaloTower
  //----------------------------------------------------
  
  CaloTower retval(id, E_em, E_had, E_outer, -1, -1, towerP4, emPoint, hadPoint);
  
  //----------------------------------------------------
  // Set timing 
  //----------------------------------------------------

  float  ecalTime = 0;
  float  hcalTime = 0;

  retval.setEcalTime(compactTime(ecalTime));
  retval.setHcalTime(compactTime(hcalTime));
  
  //----------------------------------------------------
  // Add the constituents
  //----------------------------------------------------
  /*
    std::vector<std::pair<DetId,double> > metaContains=mt.metaConstituents;
    std::vector<std::pair<DetId,double> >::iterator metaContained = metaContains.begin();
    std::vector<std::pair<DetId,double> >::iterator metaContained_end = metaContains.end();
    
    std::vector<DetId> contains;
    
    for (; metaContained != metaContained_end; ++metaContained) 
    contains.push_back(metaContained->first);
  */
  
  std::vector<DetId> contains = m_caloTowerConstituentsMap -> constituentsOf(id);

  retval.addConstituents( contains );
  
  //----------------------------------------------------
  // Return the finished CaloTower
  //----------------------------------------------------

  return retval;

} 

//----------------------------------------------------
// Decompress the energy from an HCAL TPG
//----------------------------------------------------

double CaloTowersFromTrigPrimsAlgo::getMeanEta(const HcalTrigTowerDetId& id){
  
  int ieta = id.ieta();
  int iphi = id.iphi();

  double arrayVal = m_hcalTTMeanEta[ieta + MAX_HCALTT_IETA][iphi];
  
  if ( arrayVal != -999.0 ) return arrayVal;

  double etaMin, etaMax, etaMean;

  m_hcalTrigTowerGeometry.towerEtaBounds( id.ieta(), etaMin, etaMax);

  etaMean = 0.5 * ( etaMin + etaMax );

  m_hcalTTMeanEta[ieta + MAX_HCALTT_IETA][iphi] = etaMean;
  
  return etaMean;
  
}

double CaloTowersFromTrigPrimsAlgo::getTrigTowerET(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi){

  unsigned short hcalRctInput = hcalTrigPrimDigi -> SOI_compressedEt() * 2 + hcalTrigPrimDigi -> SOI_fineGrain();
  hcalRctInput /= 2;
  
  unsigned short ietaAbs = hcalTrigPrimDigi -> id().ietaAbs();
  short sign             = hcalTrigPrimDigi -> id().zside();

  double et = (double) m_l1CaloHcalScale -> et (hcalRctInput, ietaAbs, sign);

  return et;

}

double CaloTowersFromTrigPrimsAlgo::getTrigTowerEnergy(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi){

  double et = getTrigTowerET( hcalTrigPrimDigi );

  if ( et < m_hadThreshold ) return 0.0;

  double etaMean = getMeanEta(hcalTrigPrimDigi -> id());
  
  double energy = et * TMath::CosH( etaMean );
    
  return energy;
  
}

//----------------------------------------------------
// Decompress the energy from an ECAL TPG
//----------------------------------------------------

double CaloTowersFromTrigPrimsAlgo::getMeanEta(const EcalTrigTowerDetId& id){
  
  double arrayVal = m_ecalTTMeanEta[id.ieta() + MAX_ECALTT_IETA][id.iphi()];
  
  if (arrayVal != -999.0) return arrayVal;

  double meanTheta = 0.0;
  
  std::vector<DetId> EcalDetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(id);  

  if ( EcalDetIds.size() == 0 ){
    EcalTrigTowerDetId partnerTower = getEcalPseudoTowerPartner ( id );
    EcalDetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(partnerTower);
  }

  assert ( EcalDetIds.size() != 0 );

  std::vector<DetId>::iterator ecalDetId = EcalDetIds.begin();
  std::vector<DetId>::iterator ecalDetId_end = EcalDetIds.end();
  
  for (; ecalDetId != ecalDetId_end; ++ecalDetId){       

    if ((*ecalDetId).subdetId() == EcalBarrel) {
      meanTheta += (double) m_ecalBarrelGeometry -> getGeometry( *ecalDetId ) -> getPosition().theta();
    }

    if ((*ecalDetId).subdetId() == EcalEndcap) 
      meanTheta += (double) m_ecalEndcapGeometry -> getGeometry( *ecalDetId ) -> getPosition().theta();
    
  }
  
  meanTheta /= (double) EcalDetIds.size();
  
  double meanEta = (-1.0) * TMath::Log( TMath::Tan ( (meanTheta / 2.0) ) );

  m_ecalTTMeanEta[id.ieta() + MAX_ECALTT_IETA][id.iphi()] = meanEta;

  return meanEta;
  
}

double CaloTowersFromTrigPrimsAlgo::getTrigTowerET( const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi){
  
  unsigned short ecalRctInput = ecalTrigPrimDigi -> compressedEt() * 2 + ecalTrigPrimDigi -> fineGrain();
  ecalRctInput /= 2;

  unsigned short ietaAbs = ecalTrigPrimDigi -> id().ietaAbs();
  short sign             = ecalTrigPrimDigi -> id().zside();

  double et = (double) m_l1CaloEcalScale -> et (ecalRctInput, ietaAbs, sign);

  return et;

}

double CaloTowersFromTrigPrimsAlgo::getTrigTowerEnergy(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi){  
  
  double et      = getTrigTowerET( ecalTrigPrimDigi );

  if (et < m_emThreshold) return 0.0;

  double etaMean = getMeanEta    ( ecalTrigPrimDigi -> id() );

  double energy = et * TMath::CosH( etaMean );

  return energy;

}

//----------------------------------------------------
// How many ECAL/HCAL DetId's are there in this CaloTower?
// We need this to determine how to divide the 
// CaloTower's energy among its constituents
//----------------------------------------------------

int CaloTowersFromTrigPrimsAlgo::getNDetIdsFromThisCalo (const CaloTowerDetId& id, const DetId::Detector& trigTowerDet){

  int ieta = id.ieta();
  int iphi = id.iphi();
  int temp_ieta = ieta + MAX_CT_IETA;

  if      (trigTowerDet == DetId::Ecal && m_nDetIdsInCTFromEcal[temp_ieta][iphi] != -999)
    return m_nDetIdsInCTFromEcal[temp_ieta][iphi];
  
  else if (trigTowerDet == DetId::Hcal && m_nDetIdsInCTFromHcal[temp_ieta][iphi] != -999)
    return m_nDetIdsInCTFromHcal[temp_ieta][iphi];

  else {

    int retval = 0;

    std::vector<DetId> DetIds = m_caloTowerConstituentsMap -> constituentsOf(id);
    std::vector<DetId>::iterator detId = DetIds.begin();
    std::vector<DetId>::iterator detId_end = DetIds.end();     

    for (; detId != detId_end; ++detId)
      if ((*detId).det() == trigTowerDet ) ++retval;
    
    if (trigTowerDet == DetId::Ecal) m_nDetIdsInCTFromEcal[temp_ieta][iphi] = retval;
    if (trigTowerDet == DetId::Hcal) m_nDetIdsInCTFromHcal[temp_ieta][iphi] = retval;

    return retval;

  }

}

//----------------------------------------------------
// MetaTower initializer and find function
//----------------------------------------------------

CaloTowersFromTrigPrimsAlgo::MetaTower::MetaTower() : 
  E(0),E_em(0),E_had(0),E_outer(0), emSumTimeTimesE(0), hadSumTimeTimesE(0), emSumEForTime(0), hadSumEForTime(0) {}


CaloTowersFromTrigPrimsAlgo::MetaTower & CaloTowersFromTrigPrimsAlgo::find(const CaloTowerDetId & detId) {
  MetaTowerMap::iterator itr = m_metaTowerMap.find(detId);
  if(itr == m_metaTowerMap.end()) {

    // need to build a new tower
    MetaTower t;

    // store it in the map
    m_metaTowerMap.insert(std::pair<CaloTowerDetId, CaloTowersFromTrigPrimsAlgo::MetaTower>(detId, t));
    itr = m_metaTowerMap.find(detId);
  }
  return itr->second;
}

//----------------------------------------------------
// Holdover from CaloTowersCreator
//----------------------------------------------------

int CaloTowersFromTrigPrimsAlgo::compactTime(float time) {

  const float timeUnit = 0.01; // discretization (ns)

  if (time>  300.0) return  30000;
  if (time< -300.0) return -30000;

  return int(time/timeUnit + 0.5);

}

//----------------------------------------------------
// For 'pseudo' towers on the two inner ECAL endcap 
// rings, get the tower adjacent in phi
//----------------------------------------------------

EcalTrigTowerDetId CaloTowersFromTrigPrimsAlgo::getEcalPseudoTowerPartner (const EcalTrigTowerDetId& id){

  EcalSubdetector subdet = id.subDet();
  int ieta               = id.ieta();
  int ietaAbs            = id.ietaAbs();
  int iphi               = id.iphi();
  int zside              = ieta / ietaAbs;

  if ( ietaAbs != 28 && ietaAbs != 27 ){
    edm::LogError("CaloTowersFromTrigPrimsAlgo") << " Trying to find the partner of an ECAL tower not on the inner ring!";
    return id;
  }

  bool iphi_is_even = ( iphi == (iphi/2)*2);
  
  if ( iphi_is_even) --iphi;
  if (!iphi_is_even) ++iphi;
  
  EcalTrigTowerDetId retval = EcalTrigTowerDetId(zside,subdet,ietaAbs,iphi);

  return retval;

}

//----------------------------------------------------
// Map from TriggerTowers to CaloTowers
//----------------------------------------------------

std::vector<CaloTowerDetId> CaloTowersFromTrigPrimsAlgo::getCaloTowers(const HcalTrigTowerDetId&  hcalTrigTowerDetId){

  std::vector<CaloTowerDetId> storedDetIds;
  storedDetIds = m_hcalTTMap [hcalTrigTowerDetId.ieta() + MAX_HCALTT_IETA][hcalTrigTowerDetId.iphi()];
  
  if ( !storedDetIds.empty() ) return storedDetIds;
  
  std::vector <uint32_t> rawIds;
  std::vector<CaloTowerDetId> retval;
  
  std::vector<HcalDetId> HcalDetIds = m_hcalTrigTowerGeometry.detIds(hcalTrigTowerDetId);
  std::vector<HcalDetId>::iterator hcalDetId = HcalDetIds.begin();
  std::vector<HcalDetId>::iterator hcalDetId_end = HcalDetIds.end();

  assert ( HcalDetIds.size() != 0 );
  
  for (; hcalDetId != hcalDetId_end; ++hcalDetId){
    CaloTowerDetId caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(*hcalDetId);    
    rawIds.push_back (caloTowerDetId.rawId());
  }

  std::sort (rawIds.begin(), rawIds.end());  
  rawIds.erase(std::unique(rawIds.begin(), rawIds.end()), rawIds.end());
  
  std::vector<uint32_t>::iterator rawId = rawIds.begin();
  std::vector<uint32_t>::iterator rawId_end = rawIds.end();

  for (; rawId != rawId_end; ++rawId)
    retval.push_back(CaloTowerDetId(*rawId));

  m_hcalTTMap [hcalTrigTowerDetId.ieta() + MAX_HCALTT_IETA][hcalTrigTowerDetId.iphi()] = retval;
  
  return retval;
  
}

std::vector<CaloTowerDetId> CaloTowersFromTrigPrimsAlgo::getCaloTowers(const EcalTrigTowerDetId& ecalTrigTowerDetId){
  
  std::vector<CaloTowerDetId> storedDetIds;
  storedDetIds = m_ecalTTMap [ecalTrigTowerDetId.ieta() + MAX_ECALTT_IETA][ecalTrigTowerDetId.iphi()];

  if ( !storedDetIds.empty() ) return storedDetIds;  

  std::vector <uint32_t> rawIds;
  std::vector <CaloTowerDetId> retval;

  std::vector<DetId> EcalDetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(ecalTrigTowerDetId);

  if ( EcalDetIds.size() == 0 ){
    EcalTrigTowerDetId partnerTower = getEcalPseudoTowerPartner(ecalTrigTowerDetId);
    EcalDetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(partnerTower);
  }

  assert ( EcalDetIds.size() != 0 );

  std::vector<DetId>::iterator ecalDetId = EcalDetIds.begin();
  std::vector<DetId>::iterator ecalDetId_end = EcalDetIds.end();

  for (; ecalDetId != ecalDetId_end; ++ecalDetId){
    
    CaloTowerDetId caloTowerDetId(0);

    if ((*ecalDetId).subdetId() == EcalBarrel){     
      
      EBDetId ebDetId = EBDetId(*ecalDetId);
      
      EcalTrigTowerDetId tempETTID = m_ecalTrigTowerConstituentsMap -> towerOf (ebDetId);
      caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(ebDetId);
      rawIds.push_back(caloTowerDetId.rawId());

    }

    else if ((*ecalDetId).subdetId() == EcalEndcap){
      EEDetId eeDetId = EEDetId(*ecalDetId);
      caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(eeDetId);      
      rawIds.push_back(caloTowerDetId.rawId());
    }
  }

  std::sort (rawIds.begin(), rawIds.end());  
  rawIds.erase(std::unique(rawIds.begin(), rawIds.end()), rawIds.end());
  
  std::vector<uint32_t>::iterator rawId = rawIds.begin();
  std::vector<uint32_t>::iterator rawId_end = rawIds.end();
  
  for (; rawId != rawId_end; ++rawId)
    retval.push_back(CaloTowerDetId(*rawId));

  m_ecalTTMap [ecalTrigTowerDetId.ieta() + MAX_ECALTT_IETA][ecalTrigTowerDetId.iphi()] = retval;

  return retval;

}

//----------------------------------------------------
// Set TPG L1 Calo scales
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::setL1CaloScales(const L1CaloEcalScale *ecalScale, 
						  const L1CaloHcalScale *hcalScale){
  
  
  m_l1CaloEcalScale = ecalScale;
  m_l1CaloHcalScale = hcalScale;

}

//----------------------------------------------------
// Set geometry and mapping tools
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::setGeometry( const CaloGeometry *geometry, 
					       const CaloTowerConstituentsMap *caloTowerConstituentsMap,
					       const EcalTrigTowerConstituentsMap *ecalTrigTowerConstituentsMap,
					       const CaloSubdetectorGeometry *ecalBarrelGeometry,
					       const CaloSubdetectorGeometry *ecalEndcapGeometry    ){
  
  m_caloTowerConstituentsMap     = caloTowerConstituentsMap;
  m_ecalTrigTowerConstituentsMap = ecalTrigTowerConstituentsMap;
  m_ecalBarrelGeometry           = ecalBarrelGeometry;
  m_ecalEndcapGeometry           = ecalEndcapGeometry;
  m_caloTowerGeometry            = geometry -> getSubdetectorGeometry(DetId::Calo, CaloTowerDetId::SubdetId);

}


