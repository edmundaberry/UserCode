#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAlgo.h"

//----------------------------------------------------
// Constructor
//----------------------------------------------------

CaloTowersFromTrigPrimsAlgo::CaloTowersFromTrigPrimsAlgo():
  m_geometry(0),
  m_caloTowerConstituentsMap(0),
  m_ecalTrigTowerConstituentsMap(0),
  m_ecalBarrelGeometry(0),
  m_ecalEndcapGeometry(0),
  theMomHBDepth(0.2),
  theMomHEDepth(0.4),
  theMomEBDepth(0.3),
  theMomEEDepth(0.0)
{}

//----------------------------------------------------
// Destructor
//----------------------------------------------------

CaloTowersFromTrigPrimsAlgo::~CaloTowersFromTrigPrimsAlgo(){}

//----------------------------------------------------
// Pass off the TPG handles to assign hits
//
// This can probably be templated
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::process(const HcalTrigPrimDigiCollection& HcalTrigPrimDigis) { 

  for(HcalTrigPrimDigiCollection::const_iterator hcalTrigPrimDigi = HcalTrigPrimDigis.begin();
      hcalTrigPrimDigi != HcalTrigPrimDigis.end(); ++hcalTrigPrimDigi)    
    assignHit(&(*hcalTrigPrimDigi));

}

void CaloTowersFromTrigPrimsAlgo::process(const HOTrigPrimDigiCollection& HOTrigPrimDigis) { 

  for(HOTrigPrimDigiCollection::const_iterator hoTrigPrimDigi = HOTrigPrimDigis.begin();
      hoTrigPrimDigi != HOTrigPrimDigis.end(); ++hoTrigPrimDigi)    
    assignHit(&(*hoTrigPrimDigi));

}

void CaloTowersFromTrigPrimsAlgo::process(const EcalTrigPrimDigiCollection& EcalTrigPrimDigis) { 

  for(EcalTrigPrimDigiCollection::const_iterator ecalTrigPrimDigi = EcalTrigPrimDigis.begin();
      ecalTrigPrimDigi != EcalTrigPrimDigis.end(); ++ecalTrigPrimDigi)
    assignHit(&(*ecalTrigPrimDigi));
  
}

//----------------------------------------------------
// Assign hits for the HO.
// HO trigger primitives contain only single yes/no
// bits.  Add the threshold energy.
//----------------------------------------------------

void CaloTowersFromTrigPrimsAlgo::assignHit(const HOTriggerPrimitiveDigi * hoTrigPrimDigi){
}

void CaloTowersFromTrigPrimsAlgo::assignHit(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi) {

  const HcalTrigTowerDetId hcalTrigTowerDetId = hcalTrigPrimDigi -> id();
  
  vector<CaloTowerDetId> CaloTowerDetIds = m_caloTrigTowerMap -> getCaloTowers( hcalTrigTowerDetId );
  vector<CaloTowerDetId>::iterator caloTowerDetId = CaloTowerDetIds.begin();

  double totalEnergy = getTrigTowerEnergy(hcalTrigPrimDigi);
  double towerEnergy = totalEnergy / (double) CaloTowerDetIds.size();

  for (; caloTowerDetId != CaloTowerDetIds.end(); caloTowerDetId++){

    MetaTower &metaTower = find(*caloTowerDetId);
    
    std::vector<DetId> DetIds = m_caloTowerConstituentsMap -> constituentsOf(*caloTowerDetId);
    std::vector<DetId>::iterator detId = DetIds.begin();

    int nHcalIds = 0;
    
    for (; detId != DetIds.end(); detId++)      
      if ((*detId).det() == DetId::Hcal && (*detId).subdetId() != HcalOuter) nHcalIds++;
    
    detId = DetIds.begin();

    for (; detId != DetIds.end(); detId++){
      if ((*detId).det() == DetId::Hcal && (*detId).subdetId() != HcalOuter) {
	double constituent_energy = towerEnergy / (double) nHcalIds;
	std::pair<DetId,double> mc(*detId, constituent_energy);
	metaTower.metaConstituents.push_back(mc);
      }
    }

    metaTower.E += towerEnergy;

    if ( hcalTrigTowerDetId.ieta() > m_hcalTopology.firstHFRing() ) {
      metaTower.E_had += towerEnergy / 2.0;
      metaTower.E_em  += towerEnergy / 2.0;      
    }
    
    else {
      metaTower.E_had += towerEnergy;
      metaTower.E_em  += 0.0;
    }
  }
  
}

void CaloTowersFromTrigPrimsAlgo::assignHit(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi) {

  EcalTrigTowerDetId ecalTrigTowerDetId = ecalTrigPrimDigi -> id();
  
  vector<CaloTowerDetId> CaloTowerDetIds = m_caloTrigTowerMap -> getCaloTowers( ecalTrigTowerDetId );
  vector<CaloTowerDetId>::iterator caloTowerDetId = CaloTowerDetIds.begin();

  double totalEnergy = getTrigTowerEnergy( ecalTrigPrimDigi );
  double towerEnergy = totalEnergy / (double) CaloTowerDetIds.size();

  for (; caloTowerDetId != CaloTowerDetIds.end(); caloTowerDetId++){

    MetaTower &metaTower = find(*caloTowerDetId);
    
    std::vector<DetId> DetIds = m_caloTowerConstituentsMap -> constituentsOf(*caloTowerDetId);
    std::vector<DetId>::iterator detId = DetIds.begin();

    int nEcalIds = 0;
    
    for (; detId != DetIds.end(); detId++)      
      if ((*detId).det() == DetId::Ecal) nEcalIds++;
    
    detId = DetIds.begin();

    for (; detId != DetIds.end(); detId++){
      if ((*detId).det() == DetId::Ecal) {
	double constituent_energy = towerEnergy / (double) nEcalIds;
	std::pair<DetId,double> mc(*detId, constituent_energy);
	metaTower.metaConstituents.push_back(mc);
      }
    }

    metaTower.E    += towerEnergy;
    metaTower.E_em += towerEnergy;

  }
}

void CaloTowersFromTrigPrimsAlgo::finish(CaloTowerCollection& result) {

  for(MetaTowerMap::const_iterator mapItr = theTowerMap.begin(); mapItr != theTowerMap.end(); ++mapItr) {
    
    if ( (mapItr->second).metaConstituents.size()<1) continue;
    
    CaloTower ct=convert(mapItr->first,mapItr->second);
    
    if (ct.constituentsSize()>0) result.push_back(ct);
  }

  theTowerMap.clear(); 
}

double CaloTowersFromTrigPrimsAlgo::getTrigTowerEnergy(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi){

  double et = m_caloTPGTranscoder -> hcaletValue( hcalTrigPrimDigi -> id().ieta(),
						  hcalTrigPrimDigi -> SOI_compressedEt());

  double etaMin, etaMax, etaMean;

  m_hcalTrigTowerGeometry.towerEtaBounds( hcalTrigPrimDigi -> id().ieta(), etaMin, etaMax);

  etaMean = 0.5 * ( etaMin + etaMax );

  double energy = et * TMath::CosH( etaMean );
  
  return energy;
  
}

double CaloTowersFromTrigPrimsAlgo::getTrigTowerEnergy(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi){
  
  double et = m_ecalTPGScale.getTPGInGeV(*ecalTrigPrimDigi);
  double meanEta = 0.0;

  vector<DetId> EcalDetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(ecalTrigPrimDigi -> id());
  vector<DetId>::iterator ecalDetId = EcalDetIds.begin();

  for (; ecalDetId != EcalDetIds.end(); ecalDetId++){

    if ((*ecalDetId).subdetId() == EcalBarrel) 
      meanEta += (double) m_ecalBarrelGeometry -> getGeometry( *ecalDetId ) -> getPosition().eta();
    if ((*ecalDetId).subdetId() == EcalEndcap) 
      meanEta += (double) m_ecalEndcapGeometry -> getGeometry( *ecalDetId ) -> getPosition().eta();

  }

  if (EcalDetIds.size() != 0 ) meanEta /= (double) EcalDetIds.size();
  if (EcalDetIds.size() == 0 ) return 0.0;

  double energy = et * TMath::CosH( meanEta );
  
  return energy;

}

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


CaloTowersFromTrigPrimsAlgo::MetaTower::MetaTower() : 
  E(0),E_em(0),E_had(0),E_outer(0), emSumTimeTimesE(0), hadSumTimeTimesE(0), emSumEForTime(0), hadSumEForTime(0) {}


CaloTowersFromTrigPrimsAlgo::MetaTower & CaloTowersFromTrigPrimsAlgo::find(const CaloTowerDetId & detId) {
  MetaTowerMap::iterator itr = theTowerMap.find(detId);
  if(itr == theTowerMap.end()) {

    // need to build a new tower
    MetaTower t;

    // store it in the map
    theTowerMap.insert(std::pair<CaloTowerDetId, CaloTowersFromTrigPrimsAlgo::MetaTower>(detId, t));
    itr = theTowerMap.find(detId);
  }
  return itr->second;
}

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

CaloTower CaloTowersFromTrigPrimsAlgo::convert(const CaloTowerDetId& id, const MetaTower& mt) {

  //----------------------------------------------------
  // Transfer information from meta towers  
  //----------------------------------------------------

  double E=mt.E;
  double E_em=mt.E_em;
  double E_had=mt.E_had;
  double E_outer=mt.E_outer;
  
  //----------------------------------------------------
  // Depth information.  Come back to this.
  //----------------------------------------------------
  
  if (id.ietaAbs()<=17) {
    theMomHadDepth = theMomHBDepth;
    theMomEmDepth  = theMomEBDepth;
  }
  else {
    theMomHadDepth = theMomHEDepth;
    theMomEmDepth  = theMomEEDepth;
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
  
  // For ECAL, HB/HE/HO
  if (id.ietaAbs()<=29) {
    if (E_em>0) {
      emPoint   = emShwrPos(metaContains, theMomEmDepth, E_em);
      double emPf = 1.0/cosh(emPoint.eta());
      towerP4 += CaloTower::PolarLorentzVector(E_em*emPf, emPoint.eta(), emPoint.phi(), 0); 
    }
    if (E_had>0) {
      //double E_had_tot = (theHOIsUsed && id.ietaAbs()<16)? E_had+E_outer : E_had;
      hadPoint  = hadShwrPos(metaContains, theMomHadDepth, E_had);
      double hadPf = 1.0/cosh(hadPoint.eta());
      towerP4 += CaloTower::PolarLorentzVector(E_had*hadPf, hadPoint.eta(), hadPoint.phi(), 0); 
    }
  }
  
  else {  // forward detector: use the CaloTower position 
    GlobalPoint p= m_caloTowerGeometry->getGeometry(id)->getPosition();
    double pf=1.0/cosh(p.eta());
    if (E>0) towerP4 = CaloTower::PolarLorentzVector(E*pf, p.eta(), p.phi(), 0);  // simple momentum assignment, same position
    emPoint  = p;    
    hadPoint = p;  
  }


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

int CaloTowersFromTrigPrimsAlgo::compactTime(float time) {

  const float timeUnit = 0.01; // discretization (ns)

  if (time>  300.0) return  30000;
  if (time< -300.0) return -30000;

  return int(time/timeUnit + 0.5);

}

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
