
// Libraries
#include <cmath>
#include <memory>
#include <string>

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

// Sorted collections
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"

// Energy scales
#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"
#include "CondFormats/DataRecord/interface/L1JetEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1EmEtScaleRcd.h"

// My tree objects
#include "Analyzers/L1SkimAnalyzer/interface/L1SkimTree.h"
#include "Analyzers/L1SkimAnalyzer/interface/FillL1SkimTree.h"

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
  
  void getL1CaloEmCandsInfo();
  void getL1GctJetCands_cenJets();
  void getL1GctJetCands_forJets();
  void getL1GctJetCands_tauJets();

  //-----------------------------------------------
  // Event and event setup pointers
  //-----------------------------------------------
  
  const edm::Event*      m_event;
  const edm::EventSetup* m_setup;

  //-----------------------------------------------
  // Energy scales from the event setup
  //-----------------------------------------------

  ESHandle<L1CaloEtScale> m_jetScale;
  ESHandle<L1CaloEtScale> m_emScale;
  
  //-----------------------------------------------
  // Root tree objects
  //-----------------------------------------------

  L1SkimTree     m_skimTree;
  FillL1SkimTree m_fillTree;

  //-----------------------------------------------
  // Tags
  //-----------------------------------------------

  edm::InputTag m_l1CaloEmCandsTag;
  
  edm::InputTag m_l1GctJetCands_cenJets_Tag;
  edm::InputTag m_l1GctJetCands_forJets_Tag;
  edm::InputTag m_l1GctJetCands_tauJets_Tag;
  
  
};

//-----------------------------------------------
// Constructor/destructor
//-----------------------------------------------

L1SkimAnalyzer::L1SkimAnalyzer(const edm::ParameterSet& iConfig)
{

  //-----------------------------------------------
  // Declare tags
  //-----------------------------------------------

  const edm::InputTag dL1CaloEmCandsTag("gctDigis");
  m_l1CaloEmCandsTag = iConfig.getUntrackedParameter<edm::InputTag>("l1CaloEmCandsTag",dL1CaloEmCandsTag);

  const edm::InputTag dL1GctJetCands_cenJets_Tag("gctDigis","cenJets");
  m_l1GctJetCands_cenJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>(" l1GctJetCands_cenJets_Tag",dL1GctJetCands_cenJets_Tag);

  const edm::InputTag dL1GctJetCands_tauJets_Tag("gctDigis","tauJets");
  m_l1GctJetCands_tauJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>(" l1GctJetCands_tauJets_Tag",dL1GctJetCands_tauJets_Tag);
  
  const edm::InputTag dL1GctJetCands_forJets_Tag("gctDigis","forJets");
  m_l1GctJetCands_forJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>(" l1GctJetCands_forJets_Tag",dL1GctJetCands_forJets_Tag);
 
  //-----------------------------------------------
  // Initialize the root tree
  //-----------------------------------------------

  m_fillTree.init("test.root",&m_skimTree);

}


L1SkimAnalyzer::~L1SkimAnalyzer() {}

//-----------------------------------------------
// Main analysis
//-----------------------------------------------

void
L1SkimAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //-----------------------------------------------
  // Set up the event / event setup pointers
  //-----------------------------------------------
  
  m_event = &iEvent;
  m_setup = &iSetup;

  //-----------------------------------------------
  // Get energy scales from database
  //-----------------------------------------------
  
  m_setup -> get<L1JetEtScaleRcd>().get(m_jetScale);
  m_setup -> get<L1EmEtScaleRcd> ().get(m_emScale );

  //-----------------------------------------------
  // Initialize the tree values (-999 = empty)
  //-----------------------------------------------

  m_skimTree.init();

  //-----------------------------------------------
  // Run subordinate analysis functions
  //-----------------------------------------------

  getL1CaloEmCandsInfo();
  
  getL1GctJetCands_cenJets();
  getL1GctJetCands_tauJets();
  getL1GctJetCands_forJets();

  //-----------------------------------------------
  // Fill the tree
  //-----------------------------------------------
  
  m_fillTree.fill();
  
  return;
}

void L1SkimAnalyzer::beginJob(const edm::EventSetup&) {}

void L1SkimAnalyzer::endJob(){

  m_fillTree.finalize();

}

//-----------------------------------------------
// L1 EM Calorimeter candidates
//-----------------------------------------------

void L1SkimAnalyzer::getL1CaloEmCandsInfo() {
  
  edm::Handle<L1CaloEmCollection> l1CaloEmCandsHandle;
  bool l1CaloEmCands_exist = m_event -> getByLabel(m_l1CaloEmCandsTag,l1CaloEmCandsHandle);   
  if (!l1CaloEmCands_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 Calorimeter EM Candidates! (L1CaloEmCands)";
    return;
  }

  int nL1CaloEmCands = 0;
  
  int ieta, iphi;
  int rank, raw;
  int rctCrate, rctRegion, rctCard;
  int cindex;
  int bx;
  int iso = 10;
  
  for (L1CaloEmCollection::const_iterator iL1CaloEmCand = l1CaloEmCandsHandle->begin();
       iL1CaloEmCand != l1CaloEmCandsHandle->end();
       iL1CaloEmCand++){
    
    ieta      = (int) (*iL1CaloEmCand).regionId().ieta();
    iphi      = (int) (*iL1CaloEmCand).regionId().iphi();    	      
    rank      = (int) (*iL1CaloEmCand).rank     ();
    raw       = (int) (*iL1CaloEmCand).raw      ();
    rctCrate  = (int) (*iL1CaloEmCand).rctCrate ();
    rctRegion = (int) (*iL1CaloEmCand).rctRegion();
    rctCard   = (int) (*iL1CaloEmCand).rctCard  ();    
    cindex    = (int) (*iL1CaloEmCand).index    ();
    bx        = (int) (*iL1CaloEmCand).bx       ();

    m_skimTree.caloEm_ieta      [nL1CaloEmCands] = ieta;     
    m_skimTree.caloEm_iphi      [nL1CaloEmCands] = iphi;     
    m_skimTree.caloEm_rank      [nL1CaloEmCands] = rank;     
    m_skimTree.caloEm_raw       [nL1CaloEmCands] = raw;      
    m_skimTree.caloEm_rctCrate  [nL1CaloEmCands] = rctCrate; 
    m_skimTree.caloEm_rctRegion [nL1CaloEmCands] = rctRegion;
    m_skimTree.caloEm_rctCard   [nL1CaloEmCands] = rctCard;  
    m_skimTree.caloEm_cindex    [nL1CaloEmCands] = cindex;   
    m_skimTree.caloEm_bx        [nL1CaloEmCands] = bx;       
    m_skimTree.caloEm_iso       [nL1CaloEmCands] = iso;      
    
    if ( (*iL1CaloEmCand).isolated()) iso = 1;
    else iso = 0;

    nL1CaloEmCands++;
 
  }

  m_skimTree.nCaloEm = nL1CaloEmCands;

}

//-----------------------------------------------
// Get L1 Central Jet Candidates
//-----------------------------------------------

void L1SkimAnalyzer::getL1GctJetCands_cenJets(){

  Handle<L1GctJetCandCollection> l1GctJetCands_cenJets_Handle;
  bool l1GctJetCands_cenJets_exist = m_event -> getByLabel(m_l1GctJetCands_cenJets_Tag,l1GctJetCands_cenJets_Handle);
  if (!l1GctJetCands_cenJets_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 GCT Central Jet Candidates! (L1GctJetCands)";
    return;
  }
  
  unsigned rank;
  unsigned ieta, iphi;
  unsigned capBlock, capIndex;
  int16_t bx;
  int ietaSignBit;
  int ietaSign = 10;
  double et;
  
  int nGctCenJets = 0;

  for (L1GctJetCandCollection::const_iterator iL1GctJetCand_cenJet = l1GctJetCands_cenJets_Handle -> begin();
       iL1GctJetCand_cenJet != l1GctJetCands_cenJets_Handle -> end();
       iL1GctJetCand_cenJet++){

    ietaSignBit = (int) (*iL1GctJetCand_cenJet).etaSign();
    
    if (ietaSignBit == 0) ietaSign = -1;
    if (ietaSignBit == 1) ietaSign =  1;
    
    ieta     = (*iL1GctJetCand_cenJet).etaIndex();        
    iphi     = (*iL1GctJetCand_cenJet).phiIndex();
    rank     = (*iL1GctJetCand_cenJet).rank();
    capBlock = (*iL1GctJetCand_cenJet).capBlock();
    capIndex = (*iL1GctJetCand_cenJet).capIndex();
    bx       = (*iL1GctJetCand_cenJet).bx();
    
    et = m_jetScale -> et(rank);
    
    m_skimTree.gctCenJet_ietaSign[nGctCenJets] = (int) ietaSign;
    m_skimTree.gctCenJet_ieta    [nGctCenJets] = (int) ieta;
    m_skimTree.gctCenJet_iphi    [nGctCenJets] = (int) iphi;
    m_skimTree.gctCenJet_capBlock[nGctCenJets] = (int) capBlock;
    m_skimTree.gctCenJet_capIndex[nGctCenJets] = (int) capIndex;
    m_skimTree.gctCenJet_bx      [nGctCenJets] = (int) bx;
    m_skimTree.gctCenJet_rank    [nGctCenJets] = (int) rank;
    m_skimTree.gctCenJet_et      [nGctCenJets] = (float) et;
    			      
    nGctCenJets++;
  }

  m_skimTree.nGctCenJet = nGctCenJets; 
  
}

//-----------------------------------------------
// Get L1 Forward Jet Candidates
//-----------------------------------------------

void L1SkimAnalyzer::getL1GctJetCands_forJets(){

  Handle<L1GctJetCandCollection> l1GctJetCands_forJets_Handle;
  bool l1GctJetCands_forJets_exist = m_event -> getByLabel(m_l1GctJetCands_forJets_Tag,l1GctJetCands_forJets_Handle);
  if (!l1GctJetCands_forJets_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 GCT For Jet Candidates! (L1GctJetCands)";
    return;
  }
  
  unsigned rank;
  unsigned ieta, iphi;
  unsigned capBlock, capIndex;
  int16_t bx;
  int ietaSignBit;
  int ietaSign = 10;
  double et;
  
  int nGctForJets = 0;

  for (L1GctJetCandCollection::const_iterator iL1GctJetCand_forJet = l1GctJetCands_forJets_Handle -> begin();
       iL1GctJetCand_forJet != l1GctJetCands_forJets_Handle -> end();
       iL1GctJetCand_forJet++){

    ietaSignBit = (int) (*iL1GctJetCand_forJet).etaSign();
    
    if (ietaSignBit == 0) ietaSign = -1;
    if (ietaSignBit == 1) ietaSign =  1;
    
    ieta     = (*iL1GctJetCand_forJet).etaIndex();
    iphi     = (*iL1GctJetCand_forJet).phiIndex();
    rank     = (*iL1GctJetCand_forJet).rank();
    capBlock = (*iL1GctJetCand_forJet).capBlock();
    capIndex = (*iL1GctJetCand_forJet).capIndex();
    bx       = (*iL1GctJetCand_forJet).bx();
    
    et = m_jetScale -> et(rank);
    
    m_skimTree.gctForJet_ietaSign[nGctForJets] = (int) ietaSign;
    m_skimTree.gctForJet_ieta    [nGctForJets] = (int) ieta;
    m_skimTree.gctForJet_iphi    [nGctForJets] = (int) iphi;
    m_skimTree.gctForJet_capBlock[nGctForJets] = (int) capBlock;
    m_skimTree.gctForJet_capIndex[nGctForJets] = (int) capIndex;
    m_skimTree.gctForJet_bx      [nGctForJets] = (int) bx;
    m_skimTree.gctForJet_rank    [nGctForJets] = (int) rank;
    m_skimTree.gctForJet_et      [nGctForJets] = (float) et;
    			      
    nGctForJets++;
  }

  m_skimTree.nGctForJet = nGctForJets; 
  
}

//-----------------------------------------------
// Get L1 Tau Jet Candidates
//-----------------------------------------------

void L1SkimAnalyzer::getL1GctJetCands_tauJets(){

  Handle<L1GctJetCandCollection> l1GctJetCands_tauJets_Handle;
  bool l1GctJetCands_tauJets_exist = m_event -> getByLabel(m_l1GctJetCands_tauJets_Tag,l1GctJetCands_tauJets_Handle);
  if (!l1GctJetCands_tauJets_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 GCT Tau Jet Candidates! (L1GctJetCands)";
    return;
  }
  
  unsigned rank;
  unsigned ieta, iphi;
  unsigned capBlock, capIndex;
  int16_t bx;
  int ietaSignBit;
  int ietaSign = 10;
  double et;
  
  int nGctTauJets = 0;

  for (L1GctJetCandCollection::const_iterator iL1GctJetCand_tauJet = l1GctJetCands_tauJets_Handle -> begin();
       iL1GctJetCand_tauJet != l1GctJetCands_tauJets_Handle -> end();
       iL1GctJetCand_tauJet++){

    ietaSignBit = (int) (*iL1GctJetCand_tauJet).etaSign();
    
    if (ietaSignBit == 0) ietaSign = -1;
    if (ietaSignBit == 1) ietaSign =  1;
    
    ieta     = (*iL1GctJetCand_tauJet).etaIndex();
    iphi     = (*iL1GctJetCand_tauJet).phiIndex();
    rank     = (*iL1GctJetCand_tauJet).rank();
    capBlock = (*iL1GctJetCand_tauJet).capBlock();
    capIndex = (*iL1GctJetCand_tauJet).capIndex();
    bx       = (*iL1GctJetCand_tauJet).bx();
    
    et = m_jetScale -> et(rank);
    
    m_skimTree.gctTauJet_ietaSign[nGctTauJets] = (int) ietaSign;
    m_skimTree.gctTauJet_ieta    [nGctTauJets] = (int) ieta;
    m_skimTree.gctTauJet_iphi    [nGctTauJets] = (int) iphi;
    m_skimTree.gctTauJet_capBlock[nGctTauJets] = (int) capBlock;
    m_skimTree.gctTauJet_capIndex[nGctTauJets] = (int) capIndex;
    m_skimTree.gctTauJet_bx      [nGctTauJets] = (int) bx;
    m_skimTree.gctTauJet_rank    [nGctTauJets] = (int) rank;
    m_skimTree.gctTauJet_et      [nGctTauJets] = (float) et;
    			      
    nGctTauJets++;
  }

  m_skimTree.nGctTauJet = nGctTauJets; 
  
}


DEFINE_FWK_MODULE(L1SkimAnalyzer);
