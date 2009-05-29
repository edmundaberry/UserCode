#ifndef CALOTOWERSFROMTRIGPRIMSALGO_H
#define CALOTOWERSFROMTRIGPRIMSALGO_H

// ROOT functions
#include <TMath.h>

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"

// Mapping/Geometry info
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

// Energy scaling
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibCalorimetry/EcalTPGTools/interface/EcalTPGScale.h"

// Framework
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class CaloTowersFromTrigPrimsAlgo {

 public:

  //------------------------------------------------------
  // Constructor & Destructor
  //------------------------------------------------------

  CaloTowersFromTrigPrimsAlgo();
  ~CaloTowersFromTrigPrimsAlgo();

  //------------------------------------------------------
  // Preliminary setup functions
  //------------------------------------------------------

  void setGeometry       ( const CaloGeometry *geometry, 
			   const CaloTowerConstituentsMap *caloTowerConstituentsMap,
			   const EcalTrigTowerConstituentsMap *ecalTrigTowerConstituentsMap,
			   const CaloSubdetectorGeometry *ecalBarrelGeometry,
			   const CaloSubdetectorGeometry *ecalEndcapGeometry );
  
  void setHcalTPGCoder   ( const CaloTPGTranscoder * coder ) { m_caloTPGTranscoder = coder; }
  void setEcalTPGScale   ( EcalTPGScale              scale ) { m_ecalTPGScale      = scale; }
  void setDetectorDepths ( double momHBDepth, double momHEDepth,
			   double momEBDepth, double momEEDepth );
  
  void setVerbose        ( bool v ) { m_verbose = v; }

  //------------------------------------------------------
  // Main public algorithm functions
  //------------------------------------------------------
  
  void process(const EcalTrigPrimDigiCollection& EcalTrigPrimDigis);
  void process(const HcalTrigPrimDigiCollection& HcalTrigPrimDigis);

  void finish (CaloTowerCollection& FinalCaloTowerCollection);

 private:

  //------------------------------------------------------
  // Verbosity
  //------------------------------------------------------

  bool m_verbose;

  //------------------------------------------------------
  // Energy scales
  //------------------------------------------------------
  
  EcalTPGScale              m_ecalTPGScale;
  const CaloTPGTranscoder * m_caloTPGTranscoder;  

  //------------------------------------------------------
  // Geometry/mapping info
  //------------------------------------------------------

  HcalTrigTowerGeometry               m_hcalTrigTowerGeometry;
  // CaloTrigTowerMap                   *m_caloTrigTowerMap;

  const CaloGeometry                 *m_geometry;
  const CaloTowerConstituentsMap     *m_caloTowerConstituentsMap;
  const EcalTrigTowerConstituentsMap *m_ecalTrigTowerConstituentsMap;
  const CaloSubdetectorGeometry      *m_ecalBarrelGeometry;
  const CaloSubdetectorGeometry      *m_ecalEndcapGeometry;
  const CaloSubdetectorGeometry      *m_caloTowerGeometry;
  const HcalTopology                  m_hcalTopology;

  //------------------------------------------------------
  // Detector depth information 
  //------------------------------------------------------

  double m_momHBDepth;
  double m_momHEDepth;
  double m_momEBDepth;
  double m_momEEDepth;

  //------------------------------------------------------
  // MetaTower setup
  //------------------------------------------------------
  
  template <typename TriggerPrimitiveDigi>
    void assignEnergy(const TriggerPrimitiveDigi* trigPrimDigi);
      
  struct MetaTower {
    MetaTower();
    double E, E_em, E_had, E_outer;
    std::vector< std::pair<DetId, double> > metaConstituents;
    double emSumTimeTimesE, hadSumTimeTimesE, emSumEForTime, hadSumEForTime;
  };
  
  typedef std::map<CaloTowerDetId, MetaTower> MetaTowerMap;
  MetaTowerMap m_metaTowerMap;
  
  MetaTower & find(const CaloTowerDetId &id);

 private:

  //------------------------------------------------------
  // Private algorithm functions
  //------------------------------------------------------
  
  CaloTower convert(const CaloTowerDetId& id, const MetaTower& mt);

  int compactTime(float time);

  //------------------------------------------------------
  // Unpack energy from TPG's
  //------------------------------------------------------

  double getTrigTowerEnergy(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);
  double getTrigTowerEnergy(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);
  
  //------------------------------------------------------
  // Mapping from Trigger Towers to CaloTowers
  //------------------------------------------------------

  std::vector<CaloTowerDetId> getCaloTowers(HcalTrigTowerDetId hcalTrigTowerDetId);
  std::vector<CaloTowerDetId> getCaloTowers(EcalTrigTowerDetId ecalTrigTowerDetId);

  //------------------------------------------------------
  // Shower position determination functions
  //------------------------------------------------------
  
  GlobalPoint emCrystalShwrPos (DetId detId, float fracDepth); 
  GlobalPoint hadSegmentShwrPos(DetId detId, float fracDepth);

  GlobalPoint hadShwrPos(std::vector<std::pair<DetId,double> >& metaContains, float fracDepth, double hadE  );
  GlobalPoint emShwrPos (std::vector<std::pair<DetId,double> >& metaContains, float fracDepth, double totEmE);

    
};

#endif
