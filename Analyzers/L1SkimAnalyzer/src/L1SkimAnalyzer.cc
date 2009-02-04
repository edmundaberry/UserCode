// Sorted collections
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"

// Energy scales
#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"
#include "CondFormats/DataRecord/interface/L1JetEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1EmEtScaleRcd.h"
#include "CondFormats/L1TObjects/interface/L1GctJetEtCalibrationFunction.h"
#include "CondFormats/DataRecord/interface/L1GctJetCalibFunRcd.h"

// Trigger menu and record
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

// L1 Calo Geometry
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"

// Header file
#include "Analyzers/L1SkimAnalyzer/interface/L1SkimAnalyzer.h"

// Namespaces
using namespace std; 
using namespace edm; 
using namespace l1extra;

//-----------------------------------------------
// Constructor/destructor
//-----------------------------------------------

L1SkimAnalyzer::L1SkimAnalyzer(const edm::ParameterSet& iConfig){

  //-----------------------------------------------
  // Declare tags
  //-----------------------------------------------

  const edm::InputTag dL1ExtraJetParticles_cenJets_Tag("l1extraParticles","Central");
  m_l1ExtraJetParticles_cenJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>("l1ExtraJetParticlesCenJetsTag",dL1ExtraJetParticles_cenJets_Tag);
                                   
  const edm::InputTag dL1ExtraJetParticles_tauJets_Tag("l1extraParticles","Tau");
  m_l1ExtraJetParticles_tauJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>("l1ExtraJetParticlesTauJetsTag",dL1ExtraJetParticles_tauJets_Tag);
                                   
  const edm::InputTag dL1ExtraJetParticles_forJets_Tag("l1extraParticles","Forward");
  m_l1ExtraJetParticles_forJets_Tag = iConfig.getUntrackedParameter<edm::InputTag>("l1ExtraJetParticlesForJetsTag",dL1ExtraJetParticles_forJets_Tag);
 
  const edm::InputTag dHLTRecoCaloJetCandsTag("hltIterativeCone5CaloJets");
  m_hltRecoCaloJetCandsTag = iConfig.getUntrackedParameter<edm::InputTag>("hltRecoCaloJetCandsTag",dHLTRecoCaloJetCandsTag);

  const edm::InputTag dL1GctEtHadsTag("l1GctHwDigis");
  m_l1GctEtHadsTag = iConfig.getUntrackedParameter<edm::InputTag>("l1GctEtHadsTag",dL1GctEtHadsTag);

  const edm::InputTag dL1DecisionWordTag("simGtDigis");
  m_l1DecisionWordTag = iConfig.getUntrackedParameter<edm::InputTag>("l1DecisionWordTag",dL1DecisionWordTag);

  const edm::InputTag dGenJetsTag("iterativeCone5GenJets");
  m_genJetsTag = iConfig.getUntrackedParameter<edm::InputTag>("genJetsTag",dGenJetsTag);

  //-----------------------------------------------
  // Do HLT studies?
  //-----------------------------------------------

  m_doHLT = iConfig.getUntrackedParameter<bool>("doHLT",true);

  m_hltJetThreshold = (float) iConfig.getUntrackedParameter<double>("hltJetThreshold",30.0);

  //-----------------------------------------------
  // Where should we save the root tree?
  //-----------------------------------------------
  
  m_outPath   = iConfig.getUntrackedParameter<std::string>  ("outPath","/uscms/home/eberry/data/");
  m_outSuffix = iConfig.getUntrackedParameter<std::string>  ("outSuffix","");
  m_rootFile  = m_outPath + "L1SkimAnalyzerOutput" + m_outSuffix + ".root";

  //-----------------------------------------------
  // Initialize the root tree
  //-----------------------------------------------

  m_fillTree.init(m_rootFile,&m_skimTree);

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
  // Get the run & event numbers
  //-----------------------------------------------
  
  m_event = &iEvent;
  m_setup = &iSetup;

  int run   = iEvent.id().run();
  int event = iEvent.id().event();
  
  m_skimTree.run   = run;
  m_skimTree.event = event;

  //-----------------------------------------------
  // Get ESHandles from database records
  //-----------------------------------------------
  
  // Get L1 Calo Geometry
  m_setup -> get<L1CaloGeometryRecord>().get(m_l1CaloGeometry_temp);

  // Get energy scales 
  m_setup -> get<L1JetEtScaleRcd>().get(m_jetScale);
  m_setup -> get<L1EmEtScaleRcd> ().get(m_emScale );
  m_setup -> get<L1GctJetCalibFunRcd>().get(m_jetCalibrationFunction);

  // Get trigger menu
  m_setup -> get<L1GtTriggerMenuRcd>().get(m_triggerMenu);

  //-----------------------------------------------
  // Initialize the tree values (-999 = empty)
  //-----------------------------------------------

  m_skimTree.init();

  //-----------------------------------------------
  // Run subordinate analysis functions
  //-----------------------------------------------

  getL1ExtraJetParticles_cenJets();
  getL1ExtraJetParticles_tauJets();
  getL1ExtraJetParticles_forJets();

  getL1GctEtHads();
  
  getL1DecisionWord();

  getHLTRecoCaloJetCands();

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

void L1SkimAnalyzer::getL1ExtraJetParticles_cenJets(){

  Handle< L1JetParticleCollection > l1JetParticles_cenJets_Handle;
  bool l1JetParticles_cenJets_exist = m_event -> getByLabel(m_l1ExtraJetParticles_cenJets_Tag,l1JetParticles_cenJets_Handle);
  if (!l1JetParticles_cenJets_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 Extra CenJets (L1JetParticleCollection)";
    return;
  }

  const L1CaloGeometry* l1CaloGeometry_const = &(*m_l1CaloGeometry_temp);

  unsigned rank;
  unsigned etaIndex,phiIndex;
  unsigned capBlock, capIndex;
  int16_t bx;
  int  etaSign;
  float etaMin,etaMax;
  float phiMin,phiMax;
  float px, py, pz, pt;
  float et;
  float eta, phi;

  bool isCentral = true;

  int nL1CenJets = 0;

  if (l1JetParticles_cenJets_exist){
  
    const L1GctJetCand *tempGctJetCand;
    
    for( L1JetParticleCollection::const_iterator iL1JetParticle_cenJet = l1JetParticles_cenJets_Handle -> begin();
	 iL1JetParticle_cenJet != l1JetParticles_cenJets_Handle -> end();
	 iL1JetParticle_cenJet++ ){
      
      px  = (*iL1JetParticle_cenJet).px();
      py  = (*iL1JetParticle_cenJet).py();
      pz  = (*iL1JetParticle_cenJet).pz();
      pt  = (*iL1JetParticle_cenJet).pt();
      et  = (*iL1JetParticle_cenJet).et();
      eta = (*iL1JetParticle_cenJet).eta();
      phi = (*iL1JetParticle_cenJet).phi();
      
      tempGctJetCand = (*iL1JetParticle_cenJet).gctJetCand();
      
      etaSign  = (*tempGctJetCand).etaSign();        
      etaIndex = (*tempGctJetCand).etaIndex();        
      phiIndex = (*tempGctJetCand).phiIndex();
      rank     = (*tempGctJetCand).rank();
      capBlock = (*tempGctJetCand).capBlock();
      capIndex = (*tempGctJetCand).capIndex();
      bx       = (*tempGctJetCand).bx();
      
      etaMin   = l1CaloGeometry_const -> etaBinLowEdge      ((*tempGctJetCand).etaIndex(), isCentral);
      etaMax   = l1CaloGeometry_const -> etaBinHighEdge     ((*tempGctJetCand).etaIndex(), isCentral);
      phiMin   = l1CaloGeometry_const -> emJetPhiBinLowEdge ((*tempGctJetCand).phiIndex());
      phiMax   = l1CaloGeometry_const -> emJetPhiBinHighEdge((*tempGctJetCand).phiIndex());
      
      // Physical info
      m_skimTree.l1CenJet_px      [nL1CenJets] = (float) px;
      m_skimTree.l1CenJet_py      [nL1CenJets] = (float) py;
      m_skimTree.l1CenJet_pz      [nL1CenJets] = (float) pz;
      m_skimTree.l1CenJet_pt      [nL1CenJets] = (float) pt;
      m_skimTree.l1CenJet_et      [nL1CenJets] = (float) et;
      m_skimTree.l1CenJet_eta     [nL1CenJets] = (float) eta;
      m_skimTree.l1CenJet_phi     [nL1CenJets] = (float) phi;
      
      // Digi info
      m_skimTree.l1CenJet_etaMin  [nL1CenJets] = (float) etaMin;
      m_skimTree.l1CenJet_etaMax  [nL1CenJets] = (float) etaMax;
      m_skimTree.l1CenJet_phiMin  [nL1CenJets] = (float) phiMin;
      m_skimTree.l1CenJet_phiMax  [nL1CenJets] = (float) phiMax;
      m_skimTree.l1CenJet_ietaSign[nL1CenJets] = (int) etaSign;
      m_skimTree.l1CenJet_etaIndex[nL1CenJets] = (int) etaIndex;
      m_skimTree.l1CenJet_phiIndex[nL1CenJets] = (int) phiIndex;
      m_skimTree.l1CenJet_capBlock[nL1CenJets] = (int) capBlock;
      m_skimTree.l1CenJet_capIndex[nL1CenJets] = (int) capIndex;
      m_skimTree.l1CenJet_bx      [nL1CenJets] = (int) bx;
      m_skimTree.l1CenJet_rank    [nL1CenJets] = (int) rank;
      
      nL1CenJets++;
      
    }
  }

  m_skimTree.nL1CenJet = nL1CenJets;
  
}


void L1SkimAnalyzer::getL1ExtraJetParticles_tauJets(){

  Handle< L1JetParticleCollection > l1JetParticles_tauJets_Handle;
  bool l1JetParticles_tauJets_exist = m_event -> getByLabel(m_l1ExtraJetParticles_tauJets_Tag,l1JetParticles_tauJets_Handle);
  if (!l1JetParticles_tauJets_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 Extra TauJets (L1JetParticleCollection)";
    return;
  }

  const L1CaloGeometry* l1CaloGeometry_const = &(*m_l1CaloGeometry_temp);

  unsigned rank;
  unsigned etaIndex,phiIndex;
  unsigned capBlock, capIndex;
  int16_t bx;
  int etaSign;
  float etaMin,etaMax;
  float phiMin,phiMax;
  float px, py, pz, pt;
  float et;
  float eta, phi;

  bool isCentral = true;

  int nL1TauJets = 0;
  
  if (l1JetParticles_tauJets_exist){

    const L1GctJetCand *tempGctJetCand;
    
    for( L1JetParticleCollection::const_iterator iL1JetParticle_tauJet = l1JetParticles_tauJets_Handle -> begin();
	 iL1JetParticle_tauJet != l1JetParticles_tauJets_Handle -> end();
	 iL1JetParticle_tauJet++ ){

      px  = (*iL1JetParticle_tauJet).px();
      py  = (*iL1JetParticle_tauJet).py();
      pz  = (*iL1JetParticle_tauJet).pz();
      pt  = (*iL1JetParticle_tauJet).pt();
      et  = (*iL1JetParticle_tauJet).et();
      eta = (*iL1JetParticle_tauJet).eta();
      phi = (*iL1JetParticle_tauJet).phi();
      
      tempGctJetCand = (*iL1JetParticle_tauJet).gctJetCand();
      
      etaSign  = (*tempGctJetCand).etaSign();        
      etaIndex = (*tempGctJetCand).etaIndex();        
      phiIndex = (*tempGctJetCand).phiIndex();
      rank     = (*tempGctJetCand).rank();
      capBlock = (*tempGctJetCand).capBlock();
      capIndex = (*tempGctJetCand).capIndex();
      bx       = (*tempGctJetCand).bx();
      
      etaMin   = l1CaloGeometry_const -> etaBinLowEdge      ((*tempGctJetCand).etaIndex(), isCentral);
      etaMax   = l1CaloGeometry_const -> etaBinHighEdge     ((*tempGctJetCand).etaIndex(), isCentral);
      phiMin   = l1CaloGeometry_const -> emJetPhiBinLowEdge ((*tempGctJetCand).phiIndex());
      phiMax   = l1CaloGeometry_const -> emJetPhiBinHighEdge((*tempGctJetCand).phiIndex());
      
      // Physical info
      m_skimTree.l1TauJet_px      [nL1TauJets] = (float) px;
      m_skimTree.l1TauJet_py      [nL1TauJets] = (float) py;
      m_skimTree.l1TauJet_pz      [nL1TauJets] = (float) pz;
      m_skimTree.l1TauJet_pt      [nL1TauJets] = (float) pt;
      m_skimTree.l1TauJet_et      [nL1TauJets] = (float) et;
      m_skimTree.l1TauJet_eta     [nL1TauJets] = (float) eta;
      m_skimTree.l1TauJet_phi     [nL1TauJets] = (float) phi;
      
      // Digi info
      m_skimTree.l1TauJet_etaMin  [nL1TauJets] = (float) etaMin;
      m_skimTree.l1TauJet_etaMax  [nL1TauJets] = (float) etaMax;
      m_skimTree.l1TauJet_phiMin  [nL1TauJets] = (float) phiMin;
      m_skimTree.l1TauJet_phiMax  [nL1TauJets] = (float) phiMax;
      m_skimTree.l1TauJet_ietaSign[nL1TauJets] = (int) etaSign;
      m_skimTree.l1TauJet_etaIndex[nL1TauJets] = (int) etaIndex;
      m_skimTree.l1TauJet_phiIndex[nL1TauJets] = (int) phiIndex;
      m_skimTree.l1TauJet_capBlock[nL1TauJets] = (int) capBlock;
      m_skimTree.l1TauJet_capIndex[nL1TauJets] = (int) capIndex;
      m_skimTree.l1TauJet_bx      [nL1TauJets] = (int) bx;
      m_skimTree.l1TauJet_rank    [nL1TauJets] = (int) rank;
      
      nL1TauJets++;      
    }
  }
  
  m_skimTree.nL1TauJet = nL1TauJets;
 
}

void L1SkimAnalyzer::getL1ExtraJetParticles_forJets(){

  Handle< L1JetParticleCollection > l1JetParticles_forJets_Handle;
  bool l1JetParticles_forJets_exist = m_event -> getByLabel(m_l1ExtraJetParticles_forJets_Tag,l1JetParticles_forJets_Handle);
  if (!l1JetParticles_forJets_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 Extra ForJets (L1JetParticleCollection)";
    return;
  }

  const L1CaloGeometry* l1CaloGeometry_const = &(*m_l1CaloGeometry_temp);

  unsigned rank;
  unsigned etaIndex,phiIndex;
  unsigned capBlock, capIndex;
  int16_t bx;
  int   etaSign;
  float etaMin,etaMax;
  float phiMin,phiMax;
  float px, py, pz, pt;
  float et;
  float eta, phi;

  bool isCentral = false;

  int nL1ForJets = 0;

  if (l1JetParticles_forJets_exist){
  
    const L1GctJetCand *tempGctJetCand;
    
    for( L1JetParticleCollection::const_iterator iL1JetParticle_forJet = l1JetParticles_forJets_Handle -> begin();
	 iL1JetParticle_forJet != l1JetParticles_forJets_Handle -> end();
	 iL1JetParticle_forJet++ ){
      
      px  = (*iL1JetParticle_forJet).px();
      py  = (*iL1JetParticle_forJet).py();
      pz  = (*iL1JetParticle_forJet).pz();
      pt  = (*iL1JetParticle_forJet).pt();
      et  = (*iL1JetParticle_forJet).et();
      eta = (*iL1JetParticle_forJet).eta();
      phi = (*iL1JetParticle_forJet).phi();
      
      tempGctJetCand = (*iL1JetParticle_forJet).gctJetCand();
      
      etaSign  = (*tempGctJetCand).etaSign();        
      etaIndex = (*tempGctJetCand).etaIndex();        
      phiIndex = (*tempGctJetCand).phiIndex();
      rank     = (*tempGctJetCand).rank();
      capBlock = (*tempGctJetCand).capBlock();
      capIndex = (*tempGctJetCand).capIndex();
      bx       = (*tempGctJetCand).bx();
      
      etaMin   = l1CaloGeometry_const -> etaBinLowEdge      ((*tempGctJetCand).etaIndex(), isCentral);
      etaMax   = l1CaloGeometry_const -> etaBinHighEdge     ((*tempGctJetCand).etaIndex(), isCentral);
      phiMin   = l1CaloGeometry_const -> emJetPhiBinLowEdge ((*tempGctJetCand).phiIndex());
      phiMax   = l1CaloGeometry_const -> emJetPhiBinHighEdge((*tempGctJetCand).phiIndex());

      // Physical info
      m_skimTree.l1ForJet_px      [nL1ForJets] = (float) px;
      m_skimTree.l1ForJet_py      [nL1ForJets] = (float) py;
      m_skimTree.l1ForJet_pz      [nL1ForJets] = (float) pz;
      m_skimTree.l1ForJet_pt      [nL1ForJets] = (float) pt;
      m_skimTree.l1ForJet_et      [nL1ForJets] = (float) et;
      m_skimTree.l1ForJet_eta     [nL1ForJets] = (float) eta;
      m_skimTree.l1ForJet_phi     [nL1ForJets] = (float) phi;
      
      // Digi info
      m_skimTree.l1ForJet_etaMin  [nL1ForJets] = (float) etaMin;
      m_skimTree.l1ForJet_etaMax  [nL1ForJets] = (float) etaMax;
      m_skimTree.l1ForJet_phiMin  [nL1ForJets] = (float) phiMin;
      m_skimTree.l1ForJet_phiMax  [nL1ForJets] = (float) phiMax;
      m_skimTree.l1ForJet_ietaSign[nL1ForJets] = (int) etaSign;
      m_skimTree.l1ForJet_etaIndex[nL1ForJets] = (int) etaIndex;
      m_skimTree.l1ForJet_phiIndex[nL1ForJets] = (int) phiIndex;
      m_skimTree.l1ForJet_capBlock[nL1ForJets] = (int) capBlock;
      m_skimTree.l1ForJet_capIndex[nL1ForJets] = (int) capIndex;
      m_skimTree.l1ForJet_bx      [nL1ForJets] = (int) bx;
      m_skimTree.l1ForJet_rank    [nL1ForJets] = (int) rank;
      
      nL1ForJets++;      
    }
  }
  
  m_skimTree.nL1ForJet = nL1ForJets;
 
}

//-----------------------------------------------
// Get HLT reco::caloJet info
//-----------------------------------------------

void L1SkimAnalyzer::getHLTRecoCaloJetCands(){

  //-----------------------------------------------
  // Are we even bothering?
  //-----------------------------------------------

  if (!m_doHLT) return;

  //-----------------------------------------------
  // Normal code:
  //-----------------------------------------------
  
  Handle<reco::CaloJetCollection> hltRecoCaloJetCands_Handle;
  bool hltRecoCaloJetCands_exist = m_event -> getByLabel(m_hltRecoCaloJetCandsTag,hltRecoCaloJetCands_Handle);
  if (!hltRecoCaloJetCands_exist){
    LogWarning("L1SkimAnalyzer") << "Could not extract HLT Reco Calo Jet Candidates! (CaloJets)";
  }

  Handle<reco::GenJetCollection> genJets_handle;
  bool genJets_exist =  m_event->getByLabel(m_genJetsTag,genJets_handle);
  if (!genJets_exist) {
    LogWarning("Z") << "Could not extract generated jets! ";
  }

  //-----------------------------------------------
  // Loop over the generator jet candidates
  //-----------------------------------------------

  int nGenJets = 0;

  if (genJets_exist){
  
    for (reco::GenJetCollection::const_iterator iGenJet = genJets_handle -> begin();
	 iGenJet != genJets_handle -> end();
	 iGenJet++){
      
      float p   = (float) (*iGenJet).p  ();
      float px  = (float) (*iGenJet).px ();
      float py  = (float) (*iGenJet).py ();    
      float pz  = (float) (*iGenJet).pz ();
      float pt  = (float) (*iGenJet).pt ();
      float et  = (float) (*iGenJet).et ();
      float eta = (float) (*iGenJet).eta();
      float phi = (float) (*iGenJet).phi();
      
      m_skimTree.genJet_p  [nGenJets] = p  ;
      m_skimTree.genJet_px [nGenJets] = px ;
      m_skimTree.genJet_py [nGenJets] = py ;
      m_skimTree.genJet_pz [nGenJets] = pz ;
      m_skimTree.genJet_pt [nGenJets] = pt ;
      m_skimTree.genJet_et [nGenJets] = et ;
      m_skimTree.genJet_eta[nGenJets] = eta;
      m_skimTree.genJet_phi[nGenJets] = phi;

      nGenJets++;
      
    }
  }
  
  m_skimTree.nGenJets = nGenJets;

  //-----------------------------------------------
  // Loop over the reco jet candidates
  //-----------------------------------------------
  
  int nHLTRecoCaloJets = 0;
  int nJetTowers = 0;
  
  if (hltRecoCaloJetCands_exist){

    float pt, et;
    float phi, eta;
    int   ieta, iphi;
    
    float hltHT = 0.0;
        
    for (reco::CaloJetCollection::const_iterator iHLTRecoCaloJetCand = hltRecoCaloJetCands_Handle -> begin();
	 iHLTRecoCaloJetCand != hltRecoCaloJetCands_Handle -> end();
	 iHLTRecoCaloJetCand++){
      
      pt  = (float) (*iHLTRecoCaloJetCand).pt();
      et  = (float) (*iHLTRecoCaloJetCand).et();
      phi = (float) (*iHLTRecoCaloJetCand).phi();
      eta = (float) (*iHLTRecoCaloJetCand).eta(); 
      
      if (pt > m_hltJetThreshold) hltHT += pt;

      m_skimTree.hltJet_pt [nHLTRecoCaloJets] = pt;
      m_skimTree.hltJet_et [nHLTRecoCaloJets] = et;
      m_skimTree.hltJet_phi[nHLTRecoCaloJets] = phi;
      m_skimTree.hltJet_eta[nHLTRecoCaloJets] = eta;
      
      std::vector <CaloTowerPtr> jetTowers = (*iHLTRecoCaloJetCand).getCaloConstituents();
      
      nJetTowers = 0;
      
      //-----------------------------------------------
      // Loop over the jet towers
      //-----------------------------------------------
      
      for (std::vector<CaloTowerPtr>::iterator iJetTower = jetTowers.begin();
	   iJetTower != jetTowers.end();
	   iJetTower++){
	
	ieta = (int) (**iJetTower).id().ieta();
	iphi = (int) (**iJetTower).id().iphi();
	
	m_skimTree.hltJetTower_ieta[nHLTRecoCaloJets][nJetTowers] = ieta;
	m_skimTree.hltJetTower_iphi[nHLTRecoCaloJets][nJetTowers] = iphi;
	
	nJetTowers++;
      }
      
      m_skimTree.hltHT = hltHT;
      m_skimTree.nHLTJetTowers[nHLTRecoCaloJets] = nJetTowers;
      
      nHLTRecoCaloJets++;
      
    }
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

  float ht, htUncorr;
  float htScaleLSB = (float)  m_jetCalibrationFunction -> getHtScaleLSB();
  
  int nL1GctEtHad = 0;

  if (l1GctEtHads_exist){
  
    for (L1GctEtHadCollection::const_iterator iL1GctEtHad  = l1GctEtHads_Handle -> begin();
	 iL1GctEtHad != l1GctEtHads_Handle -> end();
	 iL1GctEtHad++){
      
      // This is the method for getting HT used within
      // L1ExtraParticlesProd
      
      ht       = ( (*iL1GctEtHad).overFlow() ?
		   ( float ) L1GctEtHad::kEtHadMaxValue :
		   ( float ) (*iL1GctEtHad).et() ) * htScaleLSB + 1.e-6 ;
      
      htUncorr = ( (*iL1GctEtHad).overFlow() ?
		   ( float ) L1GctEtHad::kEtHadMaxValue :
		   ( float ) (*iL1GctEtHad).et() ) + 1.e-6 ;
      
      m_skimTree.gctHT[nL1GctEtHad] = ht;
      m_skimTree.gctHT_UnCorr[nL1GctEtHad] = htUncorr;
      
      nL1GctEtHad++;
      
    }
  }

  m_skimTree.nL1GctEtHads = nL1GctEtHad;

}

void L1SkimAnalyzer::getL1DecisionWord(){

  Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord_Handle;
  bool l1GtReadoutRecord_exists = m_event -> getByLabel(m_l1DecisionWordTag,l1GtReadoutRecord_Handle);
  if (!l1GtReadoutRecord_exists){
    LogWarning("L1SkimAnalyzer") << "Could not extract L1 Global Trigger Readout Record (L1 Decision Word)";
    return;      
  }

  const DecisionWord     decisionWord = l1GtReadoutRecord_Handle -> decisionWord(); 
  const L1GtTriggerMenu *triggerMenu  = m_triggerMenu.product();

  bool l1_HTT100 = triggerMenu -> gtAlgorithmResult("L1_HTT100",decisionWord);
  bool l1_HTT200 = triggerMenu -> gtAlgorithmResult("L1_HTT200",decisionWord);
  bool l1_HTT300 = triggerMenu -> gtAlgorithmResult("L1_HTT300",decisionWord);
  bool l1_HTT400 = triggerMenu -> gtAlgorithmResult("L1_HTT400",decisionWord);
  bool l1_HTT500 = triggerMenu -> gtAlgorithmResult("L1_HTT500",decisionWord);

  bool l1_SingleJet15 = triggerMenu -> gtAlgorithmResult("L1_SingleJet15",decisionWord);
  bool l1_SingleJet20 = triggerMenu -> gtAlgorithmResult("L1_SingleJet20",decisionWord);
  bool l1_SingleJet30 = triggerMenu -> gtAlgorithmResult("L1_SingleJet30",decisionWord);
  bool l1_SingleJet50 = triggerMenu -> gtAlgorithmResult("L1_SingleJet50",decisionWord);
  
  m_skimTree.l1_SingleJet15 = 0; if (l1_SingleJet15) m_skimTree.l1_SingleJet15 = 1;
  m_skimTree.l1_SingleJet20 = 0; if (l1_SingleJet20) m_skimTree.l1_SingleJet20 = 1;
  m_skimTree.l1_SingleJet30 = 0; if (l1_SingleJet30) m_skimTree.l1_SingleJet30 = 1;
  m_skimTree.l1_SingleJet50 = 0; if (l1_SingleJet50) m_skimTree.l1_SingleJet50 = 1;

  m_skimTree.l1_HTT100      = 0; if (l1_HTT100     ) m_skimTree.l1_HTT100      = 1;
  m_skimTree.l1_HTT200      = 0; if (l1_HTT200     ) m_skimTree.l1_HTT200      = 1;
  m_skimTree.l1_HTT300      = 0; if (l1_HTT300     ) m_skimTree.l1_HTT300      = 1;
  m_skimTree.l1_HTT400      = 0; if (l1_HTT400     ) m_skimTree.l1_HTT400      = 1;
  m_skimTree.l1_HTT500      = 0; if (l1_HTT500     ) m_skimTree.l1_HTT500      = 1;

} 

DEFINE_FWK_MODULE(L1SkimAnalyzer);
