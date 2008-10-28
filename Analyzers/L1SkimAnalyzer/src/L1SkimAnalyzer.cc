#include <cmath>
#include <memory>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

#include "Analyzers/L1SkimAnalyzer/interface/L1SkimTree.h"
#include "Analyzers/L1SkimAnalyzer/interface/FillL1SkimTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TROOT.h"

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

  //-----------------------------------------------
  // Event and setup pointers
  //-----------------------------------------------
  
  const edm::Event* m_event;
  const edm::EventSetup* m_setup;

  //-----------------------------------------------
  // Root tree objects
  //-----------------------------------------------

  L1SkimTree     m_skimTree;
  FillL1SkimTree m_fillTree;

  //-----------------------------------------------
  // Tags
  //-----------------------------------------------

  edm::InputTag l1CaloEmCandsTag_;
  
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
  l1CaloEmCandsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("l1CaloEmCandsTag",dL1CaloEmCandsTag);

  //-----------------------------------------------
  // Initialize the root tree
  //-----------------------------------------------

  m_fillTree.init("test.root",&m_skimTree);

}


L1SkimAnalyzer::~L1SkimAnalyzer()
{}

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
  // Initialize the tree values (-999 = empty)
  //-----------------------------------------------

  m_skimTree.init();

  //-----------------------------------------------
  // Run subordinate analysis functions
  //-----------------------------------------------
  
  getL1CaloEmCandsInfo();

  //-----------------------------------------------
  // Fill the tree
  //-----------------------------------------------
  
  m_fillTree.fill();
  
  return;
}

void 
L1SkimAnalyzer::beginJob(const edm::EventSetup&)
{}

void 
L1SkimAnalyzer::endJob(){

  m_fillTree.finalize();

}

//-----------------------------------------------
// L1 EM Calorimeter candidates
//-----------------------------------------------

void
L1SkimAnalyzer::getL1CaloEmCandsInfo() {
  
  edm::Handle<L1CaloEmCollection> l1CaloEmCandsHandle;
  bool l1CaloEmCands_exist = m_event -> getByLabel(l1CaloEmCandsTag_,l1CaloEmCandsHandle);   
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

DEFINE_FWK_MODULE(L1SkimAnalyzer);
