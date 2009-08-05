#ifndef CALOTOWERSFROMTRIGPRIMSALGO_H
#define CALOTOWERSFROMTRIGPRIMSALGO_H

// ROOT functions
#include <TMath.h>
#include <map>

// Data collections
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowerNonSortedCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// Mapping/Geometry info
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"

// Energy scaling
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
 
// Framework
#include "FWCore/MessageLogger/interface/MessageLogger.h"


class CaloTowersFromTrigPrimsAlgo {

 public:

  //------------------------------------------------------
  // Constructor & Destructor
  //------------------------------------------------------

  CaloTowersFromTrigPrimsAlgo(double hadThreshold,
			      double emThreshold,
			      bool useHF,
			      bool verbose);
  
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

  //------------------------------------------------------
  // Main public algorithm functions
  //------------------------------------------------------
  
  void process(const EcalTrigPrimDigiCollection& EcalTrigPrimDigis);
  void process(const HcalTrigPrimDigiCollection& HcalTrigPrimDigis);
  
  void finish (CaloTowerCollection& FinalCaloTowerCollection);
  //void finish (CaloTowerNonSortedCollection& FinalCaloTowerCollection);
  
  //------------------------------------------------------
  // Optional energy checking functions
  //------------------------------------------------------

  double getTPGEnergy(){ return m_tpgEnergy; }
  double getCCTEnergy(){ return m_cctEnergy; }

  void resetEnergy();
  void checkEnergy();
  
 private:

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
  // Instruction bools
  //------------------------------------------------------
  
  bool m_useHF, m_verbose;

  //------------------------------------------------------
  // Max ieta, iphi
  //------------------------------------------------------

  const int MAX_ECALTT_IETA, MAX_ECALTT_IPHI;
  const int MAX_HCALTT_IETA, MAX_HCALTT_IPHI;  
  const int MIN_HFTT_IETA;
  const int MAX_CT_IETA, MAX_CT_IPHI;

  //------------------------------------------------------
  // Std::Vectors of values only needed from first event
  //------------------------------------------------------
  
  std::vector < std::vector < int > > m_nDetIdsInCTFromEcal;
  std::vector < std::vector < int > > m_nDetIdsInCTFromHcal;

  std::vector < std::vector < double > > m_ecalTTCoshMeanEta;
  std::vector < std::vector < double > > m_hcalTTCoshMeanEta;

  std::vector < std::vector < std::vector < CaloTowerDetId > > > m_ecalTTMap;
  std::vector < std::vector < std::vector < CaloTowerDetId > > > m_hcalTTMap;

  std::vector < std::vector < GlobalPoint > > m_ctGlobalPoints;
  std::vector < std::vector < double > > m_ctCoshEta;

  //------------------------------------------------------
  // MetaTower setup
  //------------------------------------------------------
  
  struct MetaTower {
    MetaTower();
    double E, E_em, E_had, E_outer;
    std::vector< std::pair<DetId, double> > metaConstituents;
    std::vector<DetId> metaConstituents_withoutEnergyValues;
    double emSumTimeTimesE, hadSumTimeTimesE, emSumEForTime, hadSumEForTime;
  };
  
  typedef std::map<CaloTowerDetId, MetaTower> MetaTowerMap;
  MetaTowerMap m_metaTowerMap;

  //------------------------------------------------------
  // Private algorithm functions
  //------------------------------------------------------

  MetaTower & find(const CaloTowerDetId &id);

  template <typename TriggerPrimitiveDigi>
    void assignEnergy(const TriggerPrimitiveDigi* trigPrimDigi);
  
  CaloTower convert(const CaloTowerDetId& id, const MetaTower& mt);

  int compactTime(float time);

  int getNDetIdsFromThisCalo (const CaloTowerDetId& id, const DetId::Detector& det);

  //------------------------------------------------------
  // Unpack energy from TPG's
  //------------------------------------------------------

  double getTrigTowerEnergy(const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);
  double getTrigTowerEnergy(const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);

  double getTrigTowerET    (const EcalTriggerPrimitiveDigi * ecalTrigPrimDigi);
  double getTrigTowerET    (const HcalTriggerPrimitiveDigi * hcalTrigPrimDigi);

  double getCoshMeanEta    (const EcalTrigTowerDetId & id);
  double getCoshMeanEta    (const HcalTrigTowerDetId & id);
  
  EcalTrigTowerDetId getEcalPseudoTowerPartner ( const EcalTrigTowerDetId& id );
  
  //------------------------------------------------------
  // Mapping from Trigger Towers to CaloTowers
  //------------------------------------------------------

  void getCaloTowers(const HcalTrigTowerDetId& hcalTrigTowerDetId);
  void getCaloTowers(const EcalTrigTowerDetId& ecalTrigTowerDetId);
    
};

#endif
