// Sorted collections
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

// Energy scales
#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"
#include "CondFormats/DataRecord/interface/L1JetEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1EmEtScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1GctJetEtCalibrationFunction.h"
#include "CondFormats/DataRecord/interface/L1GctJetCalibFunRcd.h"

// L1 Calo Geometry
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"

// CaloTower Boundaries
#include "RecoJets/JetAnalyzers/interface/CaloTowerBoundries.h"

// Header file
#include "Analyzers/L1SkimAnalyzer/interface/L1SkimAnalyzer.h"

// Namespaces
using namespace std; 
using namespace edm; 

//-----------------------------------------------
// Constructor/destructor
//-----------------------------------------------

L1SkimAnalyzer::L1SkimAnalyzer(const edm::ParameterSet& iConfig)
{

  //-----------------------------------------------
  // Declare tags
  //-----------------------------------------------

  const edm::InputTag dL1CaloEmCandsTag("l1GctHwDigis");
  m_l1CaloEmCandsTag = iConfig.getUntrackedParameter<edm::InputTag>("l1CaloEmCandsTag",dL1CaloEmCandsTag);

  const edm::InputTag dL1GctJetCands_cenJets_Tag("l1GctHwDigis","cenJets");
  m_l1GctJetCands_cenJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>(" l1GctJetCands_cenJets_Tag",dL1GctJetCands_cenJets_Tag);

  const edm::InputTag dL1GctJetCands_tauJets_Tag("l1GctHwDigis","tauJets");
  m_l1GctJetCands_tauJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>(" l1GctJetCands_tauJets_Tag",dL1GctJetCands_tauJets_Tag);
  
  const edm::InputTag dL1GctJetCands_forJets_Tag("l1GctHwDigis","forJets");
  m_l1GctJetCands_forJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>(" l1GctJetCands_forJets_Tag",dL1GctJetCands_forJets_Tag);
 
  const edm::InputTag dHLTRecoCaloJetCandsTag("hltIterativeCone5CaloJets");
  m_hltRecoCaloJetCandsTag = iConfig.getUntrackedParameter<edm::InputTag>("hltRecoCaloJetCandsTag",dHLTRecoCaloJetCandsTag);

  const edm::InputTag dL1GctEtHadsTag("l1GctHwDigis");
  m_l1GctEtHadsTag = iConfig.getUntrackedParameter<edm::InputTag>("l1GctEtHadsTag",dL1GctEtHadsTag);

  //-----------------------------------------------
  // Initialize the root tree
  //-----------------------------------------------

  m_fillTree.init("/uscms/home/eberry/3DayLifetime/triggerOutput.root",&m_skimTree);

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
  // Get ESHandles from database records
  //-----------------------------------------------
  
  // Get L1 Calo Geometry
  m_setup -> get<L1CaloGeometryRecord>().get(m_l1CaloGeometry);

  // Get energy scales 
  m_setup -> get<L1JetEtScaleRcd>().get(m_jetScale);
  m_setup -> get<L1EmEtScaleRcd> ().get(m_emScale );
  m_setup -> get<L1GctJetCalibFunRcd>().get(m_jetCalibrationFunction);

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

  getHLTRecoCaloJetCands();

  getL1GctEtHads();

  //-----------------------------------------------
  // Fill the tree
  //-----------------------------------------------
  
  m_fillTree.fill();
  
  return;
}

//-----------------------------------------------
// Tell the analyzer what to do at beginning 
// and end of jobs
//-----------------------------------------------

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
  unsigned etaIndex,phiIndex;
  unsigned capBlock, capIndex;
  int16_t bx;
  int ietaMin,ietaMax;
  int iphiMin,iphiMax;
  int ietaSignBit;
  int ietaSign = 10;
  double et;

  bool isForward = false;
  
  int nGctCenJets = 0;

  for (L1GctJetCandCollection::const_iterator iL1GctJetCand_cenJet = l1GctJetCands_cenJets_Handle -> begin();
       iL1GctJetCand_cenJet != l1GctJetCands_cenJets_Handle -> end();
       iL1GctJetCand_cenJet++){

    ietaSignBit = (int) (*iL1GctJetCand_cenJet).etaSign();
    
    if (ietaSignBit == 0) ietaSign = -1;
    if (ietaSignBit == 1) ietaSign =  1;    
    
    etaIndex = (*iL1GctJetCand_cenJet).etaIndex();        
    phiIndex = (*iL1GctJetCand_cenJet).phiIndex();
    rank     = (*iL1GctJetCand_cenJet).rank();
    capBlock = (*iL1GctJetCand_cenJet).capBlock();
    capIndex = (*iL1GctJetCand_cenJet).capIndex();
    bx       = (*iL1GctJetCand_cenJet).bx();

    ietaMin  = m_l1RegionLookup.getMinHcalIeta(etaIndex,isForward);
    ietaMax  = m_l1RegionLookup.getMaxHcalIeta(etaIndex,isForward);
    iphiMin  = (int) ((phiIndex*4) + 1);
    iphiMax  = (int) (iphiMin + 3);
    
    et = m_jetScale -> et(rank);
    
    m_skimTree.gctCenJet_ietaMin [nGctCenJets] = (int) ietaMin;
    m_skimTree.gctCenJet_ietaMax [nGctCenJets] = (int) ietaMax;
    m_skimTree.gctCenJet_iphiMin [nGctCenJets] = (int) iphiMin;
    m_skimTree.gctCenJet_iphiMax [nGctCenJets] = (int) iphiMax;
    m_skimTree.gctCenJet_ietaSign[nGctCenJets] = (int) ietaSign;
    m_skimTree.gctCenJet_etaIndex[nGctCenJets] = (int) etaIndex;
    m_skimTree.gctCenJet_phiIndex[nGctCenJets] = (int) phiIndex;
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
  unsigned etaIndex,phiIndex;
  unsigned capBlock, capIndex;
  int16_t bx;
  int ietaSignBit;
  int ietaMin,ietaMax;
  int iphiMin,iphiMax;
  int ietaSign = 10;
  double et;

  bool isForward = true;
    
  int nGctForJets = 0;

  for (L1GctJetCandCollection::const_iterator iL1GctJetCand_forJet = l1GctJetCands_forJets_Handle -> begin();
       iL1GctJetCand_forJet != l1GctJetCands_forJets_Handle -> end();
       iL1GctJetCand_forJet++){

    ietaSignBit = (int) (*iL1GctJetCand_forJet).etaSign();
    
    if (ietaSignBit == 0) ietaSign =  1;
    if (ietaSignBit == 1) ietaSign = -1;
    
    etaIndex = (*iL1GctJetCand_forJet).etaIndex();
    phiIndex = (*iL1GctJetCand_forJet).phiIndex();
    rank     = (*iL1GctJetCand_forJet).rank();
    capBlock = (*iL1GctJetCand_forJet).capBlock();
    capIndex = (*iL1GctJetCand_forJet).capIndex();
    bx       = (*iL1GctJetCand_forJet).bx();

    ietaMin  = m_l1RegionLookup.getMinHcalIeta(etaIndex,isForward);
    ietaMax  = m_l1RegionLookup.getMaxHcalIeta(etaIndex,isForward);
    iphiMin  = (int) ((phiIndex*4) + 1);
    iphiMax  = (int) (iphiMin + 3);
    
    et = m_jetScale -> et(rank);
    
    m_skimTree.gctForJet_ietaMin [nGctForJets] = (int) ietaMin;
    m_skimTree.gctForJet_ietaMax [nGctForJets] = (int) ietaMax;
    m_skimTree.gctForJet_iphiMin [nGctForJets] = (int) iphiMin;
    m_skimTree.gctForJet_iphiMax [nGctForJets] = (int) iphiMax;
    m_skimTree.gctForJet_ietaSign[nGctForJets] = (int) ietaSign;
    m_skimTree.gctForJet_etaIndex[nGctForJets] = (int) etaIndex;
    m_skimTree.gctForJet_phiIndex[nGctForJets] = (int) phiIndex;
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
  unsigned etaIndex,phiIndex;
  unsigned capBlock, capIndex;
  int16_t bx;
  int ietaSignBit;
  int ietaMin,ietaMax;
  int iphiMin,iphiMax;
  int ietaSign = 10;
  double et;

  bool isForward = false;
  
  int nGctTauJets = 0;

  for (L1GctJetCandCollection::const_iterator iL1GctJetCand_tauJet = l1GctJetCands_tauJets_Handle -> begin();
       iL1GctJetCand_tauJet != l1GctJetCands_tauJets_Handle -> end();
       iL1GctJetCand_tauJet++){

    ietaSignBit = (int) (*iL1GctJetCand_tauJet).etaSign();
    
    if (ietaSignBit == 0) ietaSign =  1;
    if (ietaSignBit == 1) ietaSign = -1;
    
    etaIndex = (*iL1GctJetCand_tauJet).etaIndex();
    phiIndex = (*iL1GctJetCand_tauJet).phiIndex();
    rank     = (*iL1GctJetCand_tauJet).rank();
    capBlock = (*iL1GctJetCand_tauJet).capBlock();
    capIndex = (*iL1GctJetCand_tauJet).capIndex();
    bx       = (*iL1GctJetCand_tauJet).bx();

    ietaMin  = m_l1RegionLookup.getMinHcalIeta(etaIndex,isForward);
    ietaMax  = m_l1RegionLookup.getMaxHcalIeta(etaIndex,isForward);
    iphiMin  = (int) ((phiIndex*4) + 1);
    iphiMax  = (int) (iphiMin + 3);
    
    et = m_jetScale -> et(rank);
    
    m_skimTree.gctTauJet_ietaMin [nGctTauJets] = (int) ietaMin;
    m_skimTree.gctTauJet_ietaMax [nGctTauJets] = (int) ietaMax;
    m_skimTree.gctTauJet_iphiMin [nGctTauJets] = (int) iphiMin;
    m_skimTree.gctTauJet_iphiMax [nGctTauJets] = (int) iphiMax;
    m_skimTree.gctTauJet_ietaSign[nGctTauJets] = (int) ietaSign;
    m_skimTree.gctTauJet_etaIndex[nGctTauJets] = (int) etaIndex;
    m_skimTree.gctTauJet_phiIndex[nGctTauJets] = (int) phiIndex;
    m_skimTree.gctTauJet_capBlock[nGctTauJets] = (int) capBlock;
    m_skimTree.gctTauJet_capIndex[nGctTauJets] = (int) capIndex;
    m_skimTree.gctTauJet_bx      [nGctTauJets] = (int) bx;
    m_skimTree.gctTauJet_rank    [nGctTauJets] = (int) rank;
    m_skimTree.gctTauJet_et      [nGctTauJets] = (float) et;
    			      
    nGctTauJets++;
  }

  m_skimTree.nGctTauJet = nGctTauJets; 
  
}

//-----------------------------------------------
// Get HLT reco::caloJet info
//-----------------------------------------------

void L1SkimAnalyzer::getHLTRecoCaloJetCands(){

  Handle<reco::CaloJetCollection> hltRecoCaloJetCands_Handle;
  bool hltRecoCaloJetCAnds_exist = m_event -> getByLabel(m_hltRecoCaloJetCandsTag,hltRecoCaloJetCands_Handle);
  if (!hltRecoCaloJetCAnds_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract HLT Reco Calo Jet Candidates! (CaloJets)";
    return;
  }

  int nHLTRecoCaloJets = 0;
  int nJetTowers;
  
  float pt, et;
  float phi, eta;
  int   ieta, iphi;
  
  for (reco::CaloJetCollection::const_iterator iHLTRecoCaloJetCand = hltRecoCaloJetCands_Handle -> begin();
       iHLTRecoCaloJetCand != hltRecoCaloJetCands_Handle -> end();
       iHLTRecoCaloJetCand++){

    pt  = (float) (*iHLTRecoCaloJetCand).pt();
    et  = (float) (*iHLTRecoCaloJetCand).et();
    phi = (float) (*iHLTRecoCaloJetCand).phi();
    eta = (float) (*iHLTRecoCaloJetCand).eta(); 

    m_skimTree.hltJet_pt [nHLTRecoCaloJets] = pt;
    m_skimTree.hltJet_et [nHLTRecoCaloJets] = et;
    m_skimTree.hltJet_phi[nHLTRecoCaloJets] = phi;
    m_skimTree.hltJet_eta[nHLTRecoCaloJets] = eta;
    
    std::vector <CaloTowerPtr> jetTowers = (*iHLTRecoCaloJetCand).getCaloConstituents();

    nJetTowers = 0;

    for (std::vector<CaloTowerPtr>::iterator iJetTower = jetTowers.begin();
	 iJetTower != jetTowers.end();
	 iJetTower++){
  
      ieta = (int) (**iJetTower).id().ieta();
      iphi = (int) (**iJetTower).id().iphi();

      m_skimTree.hltJetTower_ieta[nHLTRecoCaloJets][nJetTowers] = ieta;
      m_skimTree.hltJetTower_iphi[nHLTRecoCaloJets][nJetTowers] = iphi;

      nJetTowers++;
    }
    
    m_skimTree.nHLTJetTowers[nHLTRecoCaloJets] = nJetTowers;
    
    nHLTRecoCaloJets++;

  }
  
  m_skimTree.nHLTJetCands = nHLTRecoCaloJets;
}

void L1SkimAnalyzer::getL1GctEtHads(){

  Handle<L1GctEtHadCollection> l1GctEtHads_Handle;
  bool l1GctEtHads_exist = m_event -> getByLabel(m_l1GctEtHadsTag,l1GctEtHads_Handle);
  if (!l1GctEtHads_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 GCT hadronic E_T sum (l1GctEtHads)";
    return;
  }

  float ht;
  float htScaleLSB = (float)  m_jetCalibrationFunction -> getHtScaleLSB();
  
  int nL1GctEtHad = 0;

  for (L1GctEtHadCollection::const_iterator iL1GctEtHad  = l1GctEtHads_Handle -> begin();
       iL1GctEtHad != l1GctEtHads_Handle -> end();
       iL1GctEtHad++){

    // This is the method for getting HT used within
    // L1ExtraParticlesProd

    ht = ( (*iL1GctEtHad).overFlow() ?
	   ( float ) L1GctEtHad::kEtHadMaxValue :
	   ( float ) (*iL1GctEtHad).et() ) * htScaleLSB + 1.e-6 ;
    
    cout << "***New HT values" << endl;
    cout << "   MaxValue   = " << L1GctEtHad::kEtHadMaxValue << endl;
    cout << "   etValue    = " << (*iL1GctEtHad).et()        << endl;
    cout << "   htScaleLSB = " << htScaleLSB                 << endl;
    
    m_skimTree.gctHT[nL1GctEtHad] = ht;

    nL1GctEtHad++;

  }

  m_skimTree.nL1GctEtHads = nL1GctEtHad;

}

DEFINE_FWK_MODULE(L1SkimAnalyzer);
