 // Libraries
#include <cmath>
#include <memory>
#include <string>
#include <map>

// CMSSW Framework, etc
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// ROOT objects
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TROOT.h"

// My objects
#include "Analyzers/L1SkimAnalyzer/interface/L1SkimTree.h"
#include "Analyzers/L1SkimAnalyzer/interface/FillL1SkimTree.h"
#include "Analyzers/L1SkimAnalyzer/interface/L1RegionLookup.h"

// Namespaces
using namespace std; 
using namespace edm; 

class L1SkimAnalyzer : public edm::EDAnalyzer {
public:
  explicit L1SkimAnalyzer(const edm::ParameterSet&);
  ~L1SkimAnalyzer();
  
private:

  //-----------------------------------------------
  // Primary analysis functions
  //-----------------------------------------------

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  //-----------------------------------------------
  // Subordinate analysis functions
  //-----------------------------------------------
  
  // Test EM candidates
  void getL1CaloEmCandsInfo();

  // Finding position info on GCT regions
  int  getMinHcalIeta(unsigned int etaIndex, bool isForward);
  int  getMaxHcalIeta(unsigned int etaIndex, bool isForward);

  // GCT jets
  void getL1GctJetCands_cenJets();
  void getL1GctJetCands_forJets();
  void getL1GctJetCands_tauJets();

  // HLT jets
  void getHLTRecoCaloJetCands();

  // L1 HT output 
  void getL1GctEtHads();
  
  //-----------------------------------------------
  // Event and event setup pointers
  //-----------------------------------------------
  
  const edm::Event*      m_event;
  const edm::EventSetup* m_setup;

  //-----------------------------------------------
  // ES handles
  //-----------------------------------------------

  // L1 geometry
  ESHandle<L1CaloGeometry> m_l1CaloGeometry;
  
  // Energy scales from the event setup
  ESHandle<L1CaloEtScale> m_jetScale;
  ESHandle<L1CaloEtScale> m_emScale;
  ESHandle<L1GctJetEtCalibrationFunction> m_jetCalibrationFunction;
  
  //-----------------------------------------------
  // My ROOT objects
  //-----------------------------------------------

  // Tree
  L1SkimTree     m_skimTree;
  FillL1SkimTree m_fillTree;

  // Region lookup
  L1RegionLookup m_l1RegionLookup;
  
  //-----------------------------------------------
  // Tags
  //-----------------------------------------------

  edm::InputTag m_l1CaloEmCandsTag;
  
  edm::InputTag m_l1GctJetCands_cenJets_Tag;
  edm::InputTag m_l1GctJetCands_forJets_Tag;
  edm::InputTag m_l1GctJetCands_tauJets_Tag;

  edm::InputTag m_hltRecoCaloJetCandsTag;

  edm::InputTag m_l1GctEtHadsTag;
    
};


