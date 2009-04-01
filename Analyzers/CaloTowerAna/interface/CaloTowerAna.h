// system include files
#include <memory>

// Framework files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Data format files
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

// Hcal object files
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CondFormats/HcalObjects/interface/HcalQIEShape.h"

// ROOT file filler objects
#include "Analyzers/CaloTowerAna/interface/FillCaloTowerTree.h"
#include "Analyzers/CaloTowerAna/interface/CaloTowerTree.h"

// Namespaces
using namespace std; 
using namespace edm; 

class CaloTowerAna : public edm::EDAnalyzer {
   public:
      explicit CaloTowerAna(const edm::ParameterSet&);
      ~CaloTowerAna();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      void setupHandles();
      void analyzeCaloTowers();
      void analyzeJets();

      bool matchHcalDetId(Handle<HBHEDigiCollection> & hbhedigis,
			  Handle<HODigiCollection>   & hodigis,
			  Handle<HFDigiCollection>   & hfdigis,
			  HcalDetId hcalDetId);

      void processHcalDetId(Handle<HBHEDigiCollection> & hbhedigis,
			    Handle<HODigiCollection>   & hodigis,
			    Handle<HFDigiCollection>   & hfdigis,
			    HcalDetId hcalDetId, 
			    int nrjet, int rjet_nct, int rjet_ct_ndigi);

      template <class Digi> void analyzeDigi(Digi& digi, int nrjet, int rjet_nct, int rjet_ct_ndigi);

      //---------------------------------------------
      // Event and EventSetup pointers
      //---------------------------------------------
      
      const edm::Event*      m_event;
      const edm::EventSetup* m_setup;

      //---------------------------------------------
      // edm::InputTag entries
      //---------------------------------------------

      edm::InputTag m_caloTowerTag;
      edm::InputTag m_jetTag;
      edm::InputTag m_hcalDigiTag;

      //---------------------------------------------
      // Event setup handles
      //---------------------------------------------      

      ESHandle<HcalDbService> m_conditions;

      const HcalQIEShape* m_shape;

      //---------------------------------------------
      // ROOT tree objects
      //---------------------------------------------

      FillCaloTowerTree m_fillTree;
      CaloTowerTree     m_caloTowerTree;
};

//---------------------------------------------
// Constructor
//---------------------------------------------

CaloTowerAna::CaloTowerAna(const edm::ParameterSet& iConfig){
  
  //---------------------------------------------
  // Get tags from parameter set
  //---------------------------------------------

  const edm::InputTag d_caloTowerTag("towerMaker");
  m_caloTowerTag = iConfig.getUntrackedParameter<edm::InputTag>("caloTowerTag",d_caloTowerTag);

  const edm::InputTag d_jetTag("iterativeCone5CaloJets");
  m_jetTag = iConfig.getUntrackedParameter<edm::InputTag>("jetTag",d_jetTag);

  const edm::InputTag d_hcalDigiTag("mix");
  m_hcalDigiTag = iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag",d_hcalDigiTag);

  //---------------------------------------------
  // Initialize the ROOT tree
  //---------------------------------------------

  m_fillTree.init("test.root",&m_caloTowerTree);

}

CaloTowerAna::~CaloTowerAna(){}
