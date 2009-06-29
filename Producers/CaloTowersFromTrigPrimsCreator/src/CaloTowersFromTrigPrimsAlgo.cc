#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAlgo.h"

//----------------------------------------------------
// Constructor
//----------------------------------------------------

CaloTowersFromTrigPrimsAlgo::CaloTowersFromTrigPrimsAlgo():
  m_verbose(false),
  m_geometry(0),
  m_caloTowerConstituentsMap(0),
  m_ecalTrigTowerConstituentsMap(0),
  m_ecalBarrelGeometry(0),
  m_ecalEndcapGeometry(0),
  m_momHBDepth(-999.0),
  m_momHEDepth(-999.0),
  m_momEBDepth(-999.0),
  m_momEEDepth(-999.0),
  m_hadThreshold(-1.0),
  m_emThreshold(-1.0)
{}

//----------------------------------------------------
// Destructor
//----------------------------------------------------

CaloTowersFromTrigPrimsAlgo::~CaloTowersFromTrigPrimsAlgo(){}

//----------------------------------------------------
// Set TPG L1 Calo scales
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::setL1CaloScales(const L1CaloEcalScale *ecalScale, 
						  const L1CaloHcalScale *hcalScale){
  
  
  m_l1CaloEcalScale = ecalScale;
  m_l1CaloHcalScale = hcalScale;

}

//----------------------------------------------------
// Get input from user
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::setDetectorDepths(double momHBDepth, double momHEDepth, 
						    double momEBDepth, double momEEDepth){
  
  m_momHBDepth = momHBDepth;
  m_momHEDepth = momHEDepth;
  m_momEBDepth = momEBDepth;
  m_momEEDepth = momEEDepth;

}

//----------------------------------------------------
// 1. Extract energy from trigger primitives, 
// 2. Divide energy among mapped CaloTowers
//----------------------------------------------------

template <typename TriggerPrimitiveDigi>
void CaloTowersFromTrigPrimsAlgo::assignEnergy(const TriggerPrimitiveDigi* trigPrimDigi){

  //----------------------------------------------------
  // Get the energy from this trig prim
  // If the total energy is zero, don't bother 
  // distributing it.
  //----------------------------------------------------
  
  double totalEnergy = getTrigTowerEnergy(trigPrimDigi);
  if (totalEnergy == 0.0) return;
    
  //----------------------------------------------------
  // Get the trigger tower DetId & map it to CaloTowers
  //----------------------------------------------------
  
  typedef typename TriggerPrimitiveDigi::key_type TrigTowerDetId;
  
  TrigTowerDetId trigTowerDetId = trigPrimDigi -> id();

  std::vector<CaloTowerDetId> CaloTowerDetIds = getCaloTowers( trigTowerDetId );
  std::vector<CaloTowerDetId>::iterator caloTowerDetId = CaloTowerDetIds.begin();

  double towerEnergy = totalEnergy / (double) CaloTowerDetIds.size();
  
  //----------------------------------------------------
  // Try to compare to an existing CaloTower
  //----------------------------------------------------
  
  caloTowerDetId = CaloTowerDetIds.begin();
  
  //----------------------------------------------------
  // Loop over the mapped CaloTowers and distribute
  // the energy to the appropriate MetaTowers
  //----------------------------------------------------
  
  for (; caloTowerDetId != CaloTowerDetIds.end(); caloTowerDetId++){
    
    //----------------------------------------------------
    // Find the MetaTower
    //----------------------------------------------------

    MetaTower &metaTower = find(*caloTowerDetId);
    
    std::vector<DetId> DetIds = m_caloTowerConstituentsMap -> constituentsOf(*caloTowerDetId);
    std::vector<DetId>::iterator detId = DetIds.begin();
  
    //----------------------------------------------------
    // We need to know how many DetId's from this
    // calorimeter are in this CaloTower, and we need to
    // distribute the energy evenly among them.
    //
    // Then we can give each CaloTower a list of 
    // constituent DetId's and their energies.
    // 
    // This is for consistency with the generic CaloTower
    // creation method, which used this information to 
    // determine shower positioning.
    //----------------------------------------------------
    
    int nDetIdsFromThisCalo = 0;
    
    for (; detId != DetIds.end(); detId++)      
      if ((*detId).det() == trigTowerDetId.det()) nDetIdsFromThisCalo++;
    
    detId = DetIds.begin();
    
    for (; detId != DetIds.end(); detId++){
      if ((*detId).det() == trigTowerDetId.det() && (*detId).subdetId() != HcalOuter) {
	double constituent_energy = towerEnergy / (double) nDetIdsFromThisCalo;
	std::pair<DetId,double> mc(*detId, constituent_energy);
	metaTower.metaConstituents.push_back(mc);
      }
    }
    
    //----------------------------------------------------
    // Now distribute the EM and hadronic energies.
    // E_outer is zeroed for now.
    //----------------------------------------------------
    
    metaTower.E_outer = 0;
    metaTower.E += towerEnergy;
    
    //----------------------------------------------------
    // First distribute hadronic energy
    //----------------------------------------------------
    
    if (trigTowerDetId.det() == DetId::Hcal){
      
      // Assume the HF is 50/50 EM and Hadronic energy
      
      if ( trigTowerDetId.ieta() > m_hcalTopology.firstHFRing() ) {
      
      	metaTower.E_had += towerEnergy / 2.0;
      	metaTower.E_em  += towerEnergy / 2.0;      
      }
      
      else metaTower.E_had += towerEnergy;
      // metaTower.E_had += towerEnergy;
    }
    
    //----------------------------------------------------
    // Now distribute the EM energy
    //----------------------------------------------------
    
    else if (trigTowerDetId.det() == DetId::Ecal){
      metaTower.E_em  += towerEnergy;
    }
    
    else edm::LogWarning("CaloTowersFromTrigPrimsCreator") << "This trigger tower detector ID doesn't make sense -- " << trigTowerDetId << std::endl;
    
  }
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
  m_geometry                     = geometry;
  m_ecalBarrelGeometry           = ecalBarrelGeometry;
  m_ecalEndcapGeometry           = ecalEndcapGeometry;
  m_caloTowerGeometry            = geometry -> getSubdetectorGeometry(DetId::Calo, CaloTowerDetId::SubdetId);

}

//----------------------------------------------------
// Loop through the TPG collections and assign energies
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::process(const EcalTrigPrimDigiCollection& EcalTrigPrimDigis){
					
  EcalTrigPrimDigiCollection::const_iterator trigPrimDigi = EcalTrigPrimDigis.begin();

  for(; trigPrimDigi != EcalTrigPrimDigis.end(); ++trigPrimDigi) 
    assignEnergy(&(*trigPrimDigi));

}


void CaloTowersFromTrigPrimsAlgo::process(const HcalTrigPrimDigiCollection& HcalTrigPrimDigis){					
  
  HcalTrigPrimDigiCollection::const_iterator trigPrimDigi = HcalTrigPrimDigis.begin();
  
  for(; trigPrimDigi != HcalTrigPrimDigis.end(); ++trigPrimDigi)
    assignEnergy(&(*trigPrimDigi));
    
}


//----------------------------------------------------
// 1. Loop through the map of MetaTowers
// 2. Covert MetaTowers to CaloTowers
// 3. Push CaloTowers into the Collection
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::finish(CaloTowerCollection& result){					

  for(MetaTowerMap::const_iterator mapItr = m_metaTowerMap.begin(); mapItr != m_metaTowerMap.end(); ++mapItr) {

    if ( (mapItr->second).metaConstituents.size()<1) continue;
    
    CaloTower ct=convert(mapItr->first,mapItr->second);

    if (ct.constituentsSize() > 0 && 
	( ct.emEnergy() != 0.0 || ct.hadEnergy() != 0.0 )) {      
      result.push_back(ct);
    }    
  }

  m_metaTowerMap.clear(); 
}

//----------------------------------------------------
// Decompress the energy from an HCAL TPG
//----------------------------------------------------

double CaloTowersFromTrigPrimsAlgo::getMeanEta(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi){

  double etaMin, etaMax, etaMean;

  m_hcalTrigTowerGeometry.towerEtaBounds( hcalTrigPrimDigi -> id().ieta(), etaMin, etaMax);

  etaMean = 0.5 * ( etaMin + etaMax );
  
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

  double etaMean = getMeanEta(hcalTrigPrimDigi);
  
  double energy = et * TMath::CosH( etaMean );
    
  return energy;
  
}

//----------------------------------------------------
// Decompress the energy from an ECAL TPG
//----------------------------------------------------

double CaloTowersFromTrigPrimsAlgo::getMeanEta(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi){
  
  double meanTheta = 0.0;
  
  std::vector<DetId> EcalDetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(ecalTrigPrimDigi -> id());
  std::vector<DetId>::iterator ecalDetId = EcalDetIds.begin();
  
  for (; ecalDetId != EcalDetIds.end(); ecalDetId++){
    
    if ((*ecalDetId).subdetId() == EcalBarrel) 
      meanTheta += (double) m_ecalBarrelGeometry -> getGeometry( *ecalDetId ) -> getPosition().theta();
    if ((*ecalDetId).subdetId() == EcalEndcap) 
      meanTheta += (double) m_ecalEndcapGeometry -> getGeometry( *ecalDetId ) -> getPosition().theta();
    
  }
  
  if (EcalDetIds.size() != 0 ) meanTheta /= (double) EcalDetIds.size();
  if (EcalDetIds.size() == 0 ) {
    return -999.0;
  }
  
  double meanEta = (-1.0) * TMath::Log( TMath::Tan ( (meanTheta / 2.0) ) );
  
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
  double etaMean = getMeanEta    ( ecalTrigPrimDigi );
  if (etaMean == -999.0 ){    
    return 0.0;
  }

  double energy = et * TMath::CosH( etaMean );
  
  return energy;

}

//----------------------------------------------------
// Determine the position of an EM shower
//----------------------------------------------------

GlobalPoint CaloTowersFromTrigPrimsAlgo::emShwrPos(std::vector<std::pair<DetId,double> >& metaContains, 
						   float fracDepth, double emE) {
  
  if (emE<=0) return GlobalPoint(0,0,0);
  
  double emX = 0.0;
  double emY = 0.0;
  double emZ = 0.0;

  double eSum = 0;

  std::vector<std::pair<DetId,double> >::iterator mc_it = metaContains.begin();
  for (; mc_it!=metaContains.end(); ++mc_it) {
    if (mc_it->first.det() != DetId::Ecal) continue;
    GlobalPoint p = emCrystalShwrPos(mc_it->first, fracDepth);
    double e = mc_it->second;

    if (e>0) {
      emX += p.x() * e;
      emY += p.y() * e;
      emZ += p.z() * e;
      eSum += e;
    }

  }

   return GlobalPoint(emX/eSum, emY/eSum, emZ/eSum);
}

GlobalPoint CaloTowersFromTrigPrimsAlgo::emCrystalShwrPos(DetId detId, float fracDepth) {
   const CaloCellGeometry* cellGeometry = m_geometry -> getGeometry(detId);
   GlobalPoint point = cellGeometry->getPosition();  // face of the cell

   if      (fracDepth<0) fracDepth=0;
   else if (fracDepth>1) fracDepth=1;

   if (fracDepth>0.0) {
     CaloCellGeometry::CornersVec cv = cellGeometry->getCorners();
     GlobalPoint backPoint = GlobalPoint( 0.25*( cv[4].x() + cv[5].x() + cv[6].x() + cv[7].x() ),
                                          0.25*( cv[4].y() + cv[5].y() + cv[6].y() + cv[7].y() ),
                                          0.25*( cv[4].z() + cv[5].z() + cv[6].z() + cv[7].z() ) );
     point += fracDepth * (backPoint-point);
   }

   return point;
}

//----------------------------------------------------
// Determine the position of a hadronic shower
//----------------------------------------------------

GlobalPoint CaloTowersFromTrigPrimsAlgo::hadShwrPos(std::vector<std::pair<DetId,double> >& metaContains,
                                               float fracDepth, double hadE) {
  
  // this is based on available RecHits, can lead to different actual depths if
  // hits in multi-depth towers are not all there
  if (hadE<=0) return GlobalPoint(0,0,0);

  double hadX = 0.0;
  double hadY = 0.0;
  double hadZ = 0.0;

  int nConst = 0;

  std::vector<std::pair<DetId,double> >::iterator mc_it = metaContains.begin();
  for (; mc_it!=metaContains.end(); ++mc_it) {
    if (mc_it->first.det() != DetId::Hcal) continue;
    // do not use HO for deirection calculations for now
    if (HcalDetId(mc_it->first).subdet() == HcalOuter) continue;
    ++nConst;

    GlobalPoint p = hadSegmentShwrPos(mc_it->first, fracDepth);

    // longitudinal segmentation: do not weight by energy,
    // get the geometrical position
    hadX += p.x();
    hadY += p.y();
    hadZ += p.z();
  }
  
  return GlobalPoint(hadX/nConst, hadY/nConst, hadZ/nConst);
}

GlobalPoint CaloTowersFromTrigPrimsAlgo::hadSegmentShwrPos(DetId detId, float fracDepth) {
  const CaloCellGeometry* cellGeometry = m_geometry -> getGeometry(detId);
  GlobalPoint point = cellGeometry->getPosition();  // face of the cell
  
  if      (fracDepth<0) fracDepth=0;
  else if (fracDepth>1) fracDepth=1;
  
  if (fracDepth>0.0) {
    CaloCellGeometry::CornersVec cv = cellGeometry->getCorners();
    GlobalPoint backPoint = GlobalPoint( 0.25*( cv[4].x() + cv[5].x() + cv[6].x() + cv[7].x() ),
					 0.25*( cv[4].y() + cv[5].y() + cv[6].y() + cv[7].y() ),
					 0.25*( cv[4].z() + cv[5].z() + cv[6].z() + cv[7].z() ) );
    point += fracDepth * (backPoint-point);
  }
  
  return point;
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
  // Depth information.
  //----------------------------------------------------
  
  double momHadDepth;
  double momEmDepth;

  if (id.ietaAbs()<=17) {
    momHadDepth = m_momHBDepth;
    momEmDepth  = m_momEBDepth;
  }
  else {
    momHadDepth = m_momHEDepth;
    momEmDepth  = m_momEEDepth;
  }
  
  //----------------------------------------------------
  // Timing information zeroed out
  //----------------------------------------------------  
 
  float  ecalTime = -9999.0;
  float  hcalTime = -9999.0;
  
  std::vector<std::pair<DetId,double> > metaContains=mt.metaConstituents;
  
  //----------------------------------------------------  
  // Get location information
  //----------------------------------------------------
  
  GlobalPoint emPoint, hadPoint;  
  CaloTower::PolarLorentzVector towerP4;
  
  /*
  // For ECAL, HB/HE/HO
  if (id.ietaAbs()<=29) {
    if (E_em>0) {
      emPoint   = emShwrPos(metaContains, momEmDepth, E_em);
      double emPf = 1.0/cosh(emPoint.eta());
      towerP4 += CaloTower::PolarLorentzVector(E_em*emPf, emPoint.eta(), emPoint.phi(), 0); 
    }
    if (E_had>0) {
      //double E_had_tot = (theHOIsUsed && id.ietaAbs()<16)? E_had+E_outer : E_had;
      hadPoint  = hadShwrPos(metaContains, momHadDepth, E_had);
      double hadPf = 1.0/cosh(hadPoint.eta());
      towerP4 += CaloTower::PolarLorentzVector(E_had*hadPf, hadPoint.eta(), hadPoint.phi(), 0); 
    }
  }
  */

  //else {  // forward detector: use the CaloTower position 
  GlobalPoint p= m_caloTowerGeometry->getGeometry(id)->getPosition();
  double pf=1.0/cosh(p.eta());
  if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*pf, p.eta(), p.phi(), 0);  // simple momentum assignment, same position
  emPoint  = p;    
  hadPoint = p;  
    //}


  //----------------------------------------------------
  // Create the CaloTower
  //----------------------------------------------------
  
  CaloTower retval(id, E_em, E_had, E_outer, -1, -1, towerP4, emPoint, hadPoint);
  
  //----------------------------------------------------
  // Set timing (zero this out)
  //----------------------------------------------------

  retval.setEcalTime(compactTime(ecalTime));
  retval.setHcalTime(compactTime(hcalTime));
  
  //----------------------------------------------------
  // Add the constituents
  //----------------------------------------------------
  
  std::vector<DetId> contains;
  for (std::vector<std::pair<DetId,double> >::iterator i=metaContains.begin(); i!=metaContains.end(); ++i) 
    contains.push_back(i->first);
  
  retval.addConstituents(contains);

  //----------------------------------------------------
  // Return the finished CaloTower
  //----------------------------------------------------

  return retval;

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
// Map from TriggerTowers to CaloTowers
// Formerly performed by its own class.
// Seems to run faster within this Algorithm class.
//----------------------------------------------------

std::vector<CaloTowerDetId> CaloTowersFromTrigPrimsAlgo::getCaloTowers(HcalTrigTowerDetId hcalTrigTowerDetId){
  
  std::vector <uint32_t> rawIds;
  std::vector<CaloTowerDetId> retval;
  
  std::vector<HcalDetId> HcalDetIds = m_hcalTrigTowerGeometry.detIds(hcalTrigTowerDetId);
  std::vector<HcalDetId>::iterator hcalDetId = HcalDetIds.begin();
  
  for (; hcalDetId != HcalDetIds.end(); hcalDetId++){
    CaloTowerDetId caloTowerDetId = m_caloTowerConstituentsMap -> towerOf(*hcalDetId);    
    rawIds.push_back (caloTowerDetId.rawId());
  }

  std::sort (rawIds.begin(), rawIds.end());  
  rawIds.erase(std::unique(rawIds.begin(), rawIds.end()), rawIds.end());
  
  std::vector<uint32_t>::iterator rawIds_iter = rawIds.begin();
  
  for (; rawIds_iter != rawIds.end(); rawIds_iter++)
    retval.push_back(CaloTowerDetId(*rawIds_iter));

  return retval;
  
}

std::vector<CaloTowerDetId> CaloTowersFromTrigPrimsAlgo::getCaloTowers(EcalTrigTowerDetId ecalTrigTowerDetId){
  
  std::vector <uint32_t> rawIds;
  std::vector <CaloTowerDetId> retval;

  std::vector<DetId> EcalDetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(ecalTrigTowerDetId);
  std::vector<DetId>::iterator ecalDetId = EcalDetIds.begin();

  for (; ecalDetId != EcalDetIds.end(); ecalDetId++){
    
    CaloTowerDetId caloTowerDetId(0);

    if ((*ecalDetId).subdetId() == EcalBarrel){
      EBDetId ebDetId = EBDetId(*ecalDetId);
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
  
  std::vector<uint32_t>::iterator rawIds_iter = rawIds.begin();
  
  for (; rawIds_iter != rawIds.end(); rawIds_iter++)
    retval.push_back(CaloTowerDetId(*rawIds_iter));

  return retval;

}

