#ifndef CALOTOWERSFROMTRIGPRIMSALGO_H
#define CALOTOWERSFROMTRIGPRIMSALGO_H

// ROOT functions
#include <TMath.h>
#include <map>

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
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloEcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

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

  void setGeometry        ( const CaloGeometry *geometry, 
			    const CaloTowerConstituentsMap *caloTowerConstituentsMap,
			    const EcalTrigTowerConstituentsMap *ecalTrigTowerConstituentsMap,
			    const CaloSubdetectorGeometry *ecalBarrelGeometry,
			    const CaloSubdetectorGeometry *ecalEndcapGeometry );
  
  void setDetectorDepths  ( double momHBDepth, double momHEDepth,
			    double momEBDepth, double momEEDepth );
  
  void setL1CaloScales    ( const L1CaloEcalScale *ecalScale, const L1CaloHcalScale *hcalScale);
       
  void setVerbose         ( bool  v ) { m_verbose      = v; }
  void setHadThreshold    ( float t ) { m_hadThreshold = t; }
  void setEmThreshold     ( float t ) { m_emThreshold  = t; }

  //------------------------------------------------------
  // Main public algorithm functions
  //------------------------------------------------------
  
  void setDefaultCaloTowers (const CaloTowerCollection& defaultCaloTowers);
  void process(const EcalTrigPrimDigiCollection& EcalTrigPrimDigis);
  void process(const HcalTrigPrimDigiCollection& HcalTrigPrimDigis);

  void finish (CaloTowerCollection& FinalCaloTowerCollection);

 private:

  //------------------------------------------------------
  // Verbosity
  //------------------------------------------------------

  bool m_verbose;

  //------------------------------------------------------
  // Default CaloTowers
  //------------------------------------------------------

  CaloTowerCollection m_defaultCaloTowers;

  //------------------------------------------------------
  // Energy scales
  //------------------------------------------------------
  
  const L1CaloEcalScale   * m_l1CaloEcalScale;
  const L1CaloHcalScale   * m_l1CaloHcalScale;

  //------------------------------------------------------
  // Geometry/mapping info
  //------------------------------------------------------

  HcalTrigTowerGeometry               m_hcalTrigTowerGeometry;

  const CaloGeometry                 *m_geometry;
  const CaloTowerConstituentsMap     *m_caloTowerConstituentsMap;
  const EcalTrigTowerConstituentsMap *m_ecalTrigTowerConstituentsMap;
  const CaloSubdetectorGeometry      *m_ecalBarrelGeometry;
  const CaloSubdetectorGeometry      *m_ecalEndcapGeometry;
  const CaloSubdetectorGeometry      *m_caloTowerGeometry;
  const HcalTopology                  m_hcalTopology;  

  //------------------------------------------------------
  // Threshold energies
  //------------------------------------------------------
 
  float m_hadThreshold, m_emThreshold;

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

  double getTrigTowerET(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);
  double getTrigTowerET(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);

  double getMeanEta(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);
  double getMeanEta(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);
  
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
