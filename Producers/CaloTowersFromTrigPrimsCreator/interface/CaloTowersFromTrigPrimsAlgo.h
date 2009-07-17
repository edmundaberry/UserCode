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
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

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
    
  void setL1CaloScales    ( const L1CaloEcalScale *ecalScale, const L1CaloHcalScale *hcalScale);
       
  void useHF              ( bool  b ) { m_useHF        = b; }
  void setVerbose         ( bool  v ) { m_verbose      = v; }
  void setHadThreshold    ( float t ) { m_hadThreshold = t; }
  void setEmThreshold     ( float t ) { m_emThreshold  = t; }
  void setDefaultCaloTowers ( CaloTowerCollection* c ) { m_defaultCaloTowers = c; }
  double getTPGEnergy     ()          { return m_tpgEnergy; }
  double getCCTEnergy     ()          { return m_cctEnergy; }

  void resetEnergy() { 
    m_tpgEnergy = 0.0;
    m_cctEnergy = 0.0;
  }

  void checkEnergy(){
    if (m_tpgEnergy - m_cctEnergy > 0.001 ){
      edm::LogError("CaloTowersFromTrigPrimsAlgo") << "Energy distribution error: " 
						   << "TrigPrims  had " << m_tpgEnergy << " GeV in this event, but "
						   << "CaloTowers had " << m_cctEnergy << " GeV in this event";
    }
  }

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

  CaloTowerCollection *m_defaultCaloTowers;

  //------------------------------------------------------
  // Total energy counters
  //------------------------------------------------------

  double m_tpgEnergy;
  double m_cctEnergy;

  //------------------------------------------------------
  // Energy scales
  //------------------------------------------------------
  
  const L1CaloEcalScale * m_l1CaloEcalScale;
  const L1CaloHcalScale * m_l1CaloHcalScale;

  //------------------------------------------------------
  // Geometry/mapping info
  //------------------------------------------------------

  HcalTrigTowerGeometry               m_hcalTrigTowerGeometry;

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
  // Use the HF?
  //------------------------------------------------------
  
  bool m_useHF;

  //------------------------------------------------------
  // MetaTower setup
  //------------------------------------------------------
  
  template <typename TriggerPrimitiveDigi>
    double assignEnergy(const TriggerPrimitiveDigi* trigPrimDigi);
      
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

  double getTrigTowerET    (const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);
  double getTrigTowerET    (const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);

  double getMeanEta        (const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);
  double getMeanEta        (const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);

  EcalTrigTowerDetId getEcalPseudoTowerPartner ( EcalTrigTowerDetId id );
  
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
