#ifndef PRODUCERS_CALOTOWERSFROMTRIGPRIMSCREATOR_H
#define PRODUCERS_CALOTOWERSFROMTRIGPRIMSCREATOR_H

// System include files
#include <memory>
#include <TMath.h>

// Framework header files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

// Energy scale info
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloEcalScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

// Mapping/Geometry info
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

// My ROOT tree
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAnalyzerTree.h"
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/FillCaloTowersFromTrigPrimsAnalyzerTree.h"

class CaloTowersFromTrigPrimsAnalyzer : public edm::EDAnalyzer {
public:
  explicit CaloTowersFromTrigPrimsAnalyzer(const edm::ParameterSet&);
  ~CaloTowersFromTrigPrimsAnalyzer();
  
  
private:

  //------------------------------------------------------
  // Predefined analysis functions
  //------------------------------------------------------

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  //------------------------------------------------------
  // Analyzer functions
  //------------------------------------------------------

  void analyzeTPGs();
  void analyzeCaloTowers();
  void analyzeJets();

  void analyzeJets( const reco::GenJetCollection & Jets );
  void analyzeJets( const reco::CaloJetCollection & Jets, bool fromTPGs, bool corrected );
  
  float analyzeCaloTowers( const edm::View<reco::Candidate> & towers, bool isMine);
  
  template <typename TrigPrimDigiCollection>  
  float analyzeTPGs( const TrigPrimDigiCollection& trigPrimDigiCollection );

  //------------------------------------------------------
  // Helper functions for analyzing TPG energies
  //------------------------------------------------------
  
  EcalTrigTowerDetId getEcalPseudoTowerPartner (EcalTrigTowerDetId id);

  float getTrigTowerET (const HcalTriggerPrimitiveDigi& hcalTrigPrimDigi);
  float getTrigTowerET (const EcalTriggerPrimitiveDigi& ecalTrigPrimDigi);

  float getTrigTowerMeanEta (const HcalTriggerPrimitiveDigi& hcalTrigPrimDigi);
  float getTrigTowerMeanEta (const EcalTriggerPrimitiveDigi& ecalTrigPrimDigi);

  //------------------------------------------------------
  // TPG tags
  //------------------------------------------------------
  
  edm::InputTag m_hcalTrigPrimTag;
  edm::InputTag m_ecalTrigPrimTag;

  edm::InputTag m_createdCaloTowerTag;
  edm::InputTag m_defaultCaloTowerTag;

  edm::InputTag m_genJetTag;	
  edm::InputTag m_tpgJetTag;	
  edm::InputTag m_caloJetTag;			                
  edm::InputTag m_tpgCorJetTag;	
  edm::InputTag m_caloCorJetTag;

  //------------------------------------------------------
  // Energy scales
  //------------------------------------------------------

  edm::ESHandle<L1CaloEcalScale>   m_l1CaloEcalScale;
  edm::ESHandle<L1CaloHcalScale>   m_l1CaloHcalScale;

  //------------------------------------------------------
  // Geometry objects
  //------------------------------------------------------

  edm::ESHandle<CaloSubdetectorGeometry>      m_ecalBarrelGeometry;
  edm::ESHandle<CaloSubdetectorGeometry>      m_ecalEndcapGeometry;
  HcalTrigTowerGeometry                       m_hcalTrigTowerGeometry;

  //------------------------------------------------------
  // Constituent mappings
  //------------------------------------------------------
  
  edm::ESHandle<EcalTrigTowerConstituentsMap> m_ecalTrigTowerConstituentsMap;
  
  //---------------------------------------------
  // ROOT tree objects
  //---------------------------------------------
  
  CaloTowersFromTrigPrimsAnalyzerTree m_caloTowersFromTrigPrimsAnalyzerTree;
  FillCaloTowersFromTrigPrimsAnalyzerTree m_fillCaloTowersFromTrigPrimsAnalyzerTree;
  
  //---------------------------------------------
  // Event and EventSetup pointers
  //---------------------------------------------
  
  const edm::Event*      m_event;
  const edm::EventSetup* m_setup;

  bool m_verbose;
  bool m_newEvent;
  int m_tpgNumber;

  float m_tpgEnergy;
  float m_cctEnergy;


};

#endif
