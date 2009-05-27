#ifndef CALOTOWERSFROMTRIGPRIMSALGO_H
#define CALOTOWERSFROMTRIGPRIMSALGO_H

#include <TMath.h>

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"

// Mapping function
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTrigTowerMap.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

// TPG transcoder info
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibCalorimetry/EcalTPGTools/interface/EcalTPGScale.h"

class CaloTowersFromTrigPrimsAlgo {

 public:

  

  CaloTowersFromTrigPrimsAlgo();
  ~CaloTowersFromTrigPrimsAlgo();

  void setGeometry     ( const CaloGeometry *geometry, 
		         const CaloTowerConstituentsMap *caloTowerConstituentsMap,
		         const EcalTrigTowerConstituentsMap *ecalTrigTowerConstituentsMap,
			 const CaloSubdetectorGeometry *ecalBarrelGeometry,
			 const CaloSubdetectorGeometry *ecalEndcapGeometry );
  		       
  void setMap          ( CaloTrigTowerMap * map   ) { m_caloTrigTowerMap  = map;   }
  void setCoder        ( const CaloTPGTranscoder * coder ) { m_caloTPGTranscoder = coder; }
  void setEcalTPGScale ( EcalTPGScale              scale ) { m_ecalTPGScale      = scale; }

  void process(const HcalTrigPrimDigiCollection& HcalTrigPrimDigis);
  void process(const HOTrigPrimDigiCollection&   HOTrigPrimDigis  );
  void process(const EcalTrigPrimDigiCollection& EcalTrigPrimDigis);
  void finish (CaloTowerCollection& FinalCaloTowerCollection);

 private:
  
  EcalTPGScale              m_ecalTPGScale;
  HcalTrigTowerGeometry     m_hcalTrigTowerGeometry;
  CaloTrigTowerMap        * m_caloTrigTowerMap;
  const CaloTPGTranscoder * m_caloTPGTranscoder;  

  struct MetaTower {
    MetaTower();
    double E, E_em, E_had, E_outer;
    std::vector< std::pair<DetId, double> > metaConstituents;
    double emSumTimeTimesE, hadSumTimeTimesE, emSumEForTime, hadSumEForTime;
  };

  typedef std::map<CaloTowerDetId, MetaTower> MetaTowerMap;
  MetaTowerMap theTowerMap;

  const CaloGeometry                 *m_geometry;
  const CaloTowerConstituentsMap     *m_caloTowerConstituentsMap;
  const EcalTrigTowerConstituentsMap *m_ecalTrigTowerConstituentsMap;
  const CaloSubdetectorGeometry      *m_ecalBarrelGeometry;
  const CaloSubdetectorGeometry      *m_ecalEndcapGeometry;
  const CaloSubdetectorGeometry      *m_caloTowerGeometry;
  const HcalTopology                  m_hcalTopology;

  MetaTower & find(const CaloTowerDetId &id);


  void assignHit(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);
  void assignHit(const HOTriggerPrimitiveDigi   * hoTrigPrimDigi  );
  void assignHit(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);

  CaloTower convert(const CaloTowerDetId& id, const MetaTower& mt);

  int compactTime(float time);

  double getTrigTowerEnergy(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);
  double getTrigTowerEnergy(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);
  
  GlobalPoint emCrystalShwrPos (DetId detId, float fracDepth); 
  GlobalPoint hadSegmentShwrPos(DetId detId, float fracDepth);

  GlobalPoint hadShwrPos(std::vector<std::pair<DetId,double> >& metaContains, float fracDepth, double hadE  );
  GlobalPoint emShwrPos (std::vector<std::pair<DetId,double> >& metaContains, float fracDepth, double totEmE);


  double theMomHadDepth;
  double theMomEmDepth;

  double theMomHBDepth;
  double theMomHEDepth;   
  double theMomEBDepth;
  double theMomEEDepth;
  
  

};

#endif
