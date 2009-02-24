// -*- C++ -*-
//
// Package:    DijetAnalyzer
// Class:      DijetAnalyzer
// 
/**\class DijetAnalyzer DijetAnalyzer.cc Analyzers/DijetAnalyzer/src/DijetAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Edmund Berry"
//         Created:  Wed Jun 11 14:17:06 CDT 2008
// $Id: DijetAnalyzer.cc,v 1.2 2008/09/12 16:32:42 eberry Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Analyzers/DijetAnalyzer/interface/FillMyTree.h"
#include "Analyzers/DijetAnalyzer/interface/MyTree.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtSums.h"

#include "TH1.h"
#include "TFile.h"

class DijetAnalyzer : public edm::EDAnalyzer {
public:
  explicit DijetAnalyzer(const edm::ParameterSet&);
  ~DijetAnalyzer();
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  
  FillMyTree m_fillTree;
  MyTree     m_myTree;

  edm::InputTag calJetsTag_;
  edm::InputTag caloTowerTag_;

  edm::InputTag hbheDigiTag_;
  edm::InputTag hoDigiTag_;
  edm::InputTag hfDigiTag_;

  edm::InputTag l1ForJetTag_;
  edm::InputTag l1CenJetTag_;
  edm::InputTag l1TauJetTag_;

  edm::InputTag hbheRecHitTag_;
  edm::InputTag hoRecHitTag_;
  edm::InputTag hfRecHitTag_;

  std::string rootFile_;
  
  TFile *outFile;
  
};


DijetAnalyzer::DijetAnalyzer(const edm::ParameterSet& iConfig):
  rootFile_(iConfig.getUntrackedParameter<std::string>("rootFile","DijetAnalyzer_output.root"))
{

  const edm::InputTag dCalJetTag("iterativeCone5CaloJets");
  calJetsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("calJetsTag",dCalJetTag);

  const edm::InputTag dCaloTowerTag("towerMaker");
  caloTowerTag_ = iConfig.getUntrackedParameter<edm::InputTag>("caloTowerTag",dCaloTowerTag);

  const edm::InputTag dHBHEDigiTag("hcalDigis");
  hbheDigiTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag",dHBHEDigiTag);

  const edm::InputTag dHODigiTag  ("hcalDigis");
  hoDigiTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag",dHODigiTag);

  const edm::InputTag dHFDigiTag  ("hcalDigis");
  hfDigiTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag",dHFDigiTag);

  const edm::InputTag dGctCenJetsTag ("l1GctHwDigis","cenJets");
  l1CenJetTag_ = iConfig.getUntrackedParameter<edm::InputTag>("cenJetTag",dGctCenJetsTag);
  
  const edm::InputTag dGctForJetsTag ("l1GctHwDigis","forJets");
  l1ForJetTag_ = iConfig.getUntrackedParameter<edm::InputTag>("forJetTag",dGctForJetsTag);
  
  const edm::InputTag dGctTauJetsTag ("l1GctHwDigis","tauJets");
  l1TauJetTag_ = iConfig.getUntrackedParameter<edm::InputTag>("tauJetTag",dGctTauJetsTag);
  
  const edm::InputTag dHBHERecHitTag ("hbhereco");
  hbheRecHitTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hbheRecHitTag",dHBHERecHitTag);
  
  const edm::InputTag dHORecHitTag ("horeco");
  hoRecHitTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hoRecHitTag",dHORecHitTag);

  const edm::InputTag dHFRecHitTag ("hfreco");
  hfRecHitTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hfRecHitTag",dHFRecHitTag);

  m_fillTree.init(rootFile_,&m_myTree);

}

DijetAnalyzer::~DijetAnalyzer()
{}

void
DijetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  //-----------------------------------------------------
  // Declare namespaces and initialize tree
  //-----------------------------------------------------

  using namespace edm;
  using namespace reco;
  using namespace std;

  m_myTree.init();

  //-----------------------------------------------------
  // Get conditions for ADC2fC coding
  //-----------------------------------------------------

  ESHandle<HcalDbService> conditions;
  iSetup.get<HcalDbRecord>().get(conditions);
  
  const HcalQIEShape* shape = conditions->getHcalShape();
  CaloSamples tool;
 
  //-----------------------------------------------------
  // Get all handles for data types. Catch exceptions.
  //-----------------------------------------------------
  
  bool gotHandle;

  Handle<CaloJetCollection> jets;
  gotHandle =  iEvent.getByLabel(calJetsTag_,jets);
  if (!gotHandle){
    LogWarning("DijetAnalyzer") << "Could not extract reco jets with the tag" << calJetsTag_;
  }
  
  Handle<CaloTowerCollection> caloTowers;
  gotHandle =  iEvent.getByLabel(caloTowerTag_,caloTowers);
  if (!gotHandle){
    LogWarning("DijetAnalyzer") << "Could not extract calo towers";
  }
  
  Handle<HBHERecHitCollection> HBHERecHits;
  gotHandle =  iEvent.getByLabel(hbheRecHitTag_ ,HBHERecHits);
  if (!gotHandle){
    LogWarning("DijetAnalyzer") << "Could not extract HBHE Rec Hits";
  }

  Handle<HORecHitCollection> HORecHits;
  gotHandle =  iEvent.getByLabel(hoRecHitTag_ ,HORecHits);
  if (!gotHandle){
    LogWarning("DijetAnalyzer") << "Could not extract HO Rec Hits";
  }
  
  Handle<HFRecHitCollection> HFRecHits;
  gotHandle =  iEvent.getByLabel(hfRecHitTag_ ,HFRecHits);
  if (!gotHandle){
    LogWarning("DijetAnalyzer") << "Could not extract HF Rec Hits";
  }
  
  Handle<HBHEDigiCollection> HBHEDigis;
  gotHandle =  iEvent.getByLabel(hbheDigiTag_,HBHEDigis);
  if (!gotHandle){
    LogWarning("DijetAnalyzer") << "Could not extract HB/HE Digis";
  }
  
  Handle<HODigiCollection> HODigis;
  gotHandle =  iEvent.getByLabel(hoDigiTag_,HODigis);
  if (!gotHandle){
    LogWarning("DijetAnalyzer") << "Could not extract HO Digis";
  }
  
  Handle<HFDigiCollection> HFDigis;
  gotHandle =  iEvent.getByLabel(hfDigiTag_,HFDigis);
  if (!gotHandle){
    LogWarning("DijetAnalyzer") << "Could not extract HF Digis";
  }

  //-----------------------------------------------------
  // Declare important values
  //-----------------------------------------------------

  int run        = iEvent.id().run();
  int event      = iEvent.id().event();
  int nrjet      = 0;
  int ngtower    = 0;
  int ngtower_con = 0;

  float hadEnergy,hadEt,emEnergy,emEt;
  float p, px, py, pz, pt;
  float e, et;
  float eta, phi;
  float fC_fillVal;

  int ieta, ietaAbs, iphi, zside, depth;
  int ADC_fillVal;  
  int ntslice;
  int njtower;
  int nrhits;
  int subdet;

  int genTower_con_subdet, genTower_con_ieta, 
    genTower_con_iphi, genTower_con_zside, genTower_con_depth;

  //-----------------------------------------------------
  // This is the kind of iteration statement you should
  // use to loop over the L1GctJetCandCollection
  //-----------------------------------------------------

  //L1GctJetCandCollection::const_iterator iL1ForJet = l1ForJets -> begin();
  //iL1ForJet != l1ForJets -> end();
  //iL1ForJet ++ ){}
  
  //-----------------------------------------------------
  // Loop over all towers
  //-----------------------------------------------------

  if (caloTowers.isValid()){
    
    for (CaloTowerCollection::const_iterator iGenTower = caloTowers->begin();
	 iGenTower != caloTowers->end();
	 iGenTower ++){
      
      ieta      = (int)   (*iGenTower).id().ieta();
      ietaAbs   = (int)   (*iGenTower).id().ietaAbs();
      iphi      = (int)   (*iGenTower).id().iphi();
      zside     = (int)   (*iGenTower).id().zside();
      
      eta       = (float) (*iGenTower).eta();
      phi       = (float) (*iGenTower).phi();
      
      hadEnergy = (float) (*iGenTower).hadEnergy();
      hadEt     = (float) (*iGenTower).hadEt();
      
      emEnergy  = (float) (*iGenTower).emEnergy();
      emEt      = (float) (*iGenTower).emEnergy();
      
      m_myTree.genTower_eta      [ngtower] = eta;
      m_myTree.genTower_phi      [ngtower] = phi;
      
      m_myTree.genTower_ieta     [ngtower] = ieta;
      m_myTree.genTower_ietaAbs  [ngtower] = ietaAbs;
      m_myTree.genTower_iphi     [ngtower] = iphi;
      
      m_myTree.genTower_zside    [ngtower] = zside    ;
      m_myTree.genTower_hadEnergy[ngtower] = hadEnergy;
      m_myTree.genTower_hadEt    [ngtower] = hadEt    ;
      m_myTree.genTower_emEnergy [ngtower] = emEnergy ;
      m_myTree.genTower_emEt     [ngtower] = emEt     ;
      
      //-----------------------------------------------------
      // I want depth information, so I need to loop over
      // the tower constituents.
      //-----------------------------------------------------

      size_t nTowerCons = (*iGenTower).constituentsSize();

      for (size_t iTowerCon = 0; iTowerCon < nTowerCons; iTowerCon++){
	
	DetId towerConId = (*iGenTower).constituent(iTowerCon);

	//-----------------------------------------------------
	// Is this in the HCAL?
	//-----------------------------------------------------

	DetId::Detector DetNum = towerConId.det();
       	
	if (DetNum == DetId::Hcal){

	  HcalDetId towerConHcalId = HcalDetId(towerConId);  

	  genTower_con_subdet = (int) towerConHcalId.subdet();
	  genTower_con_ieta   = (int) towerConHcalId.ieta();
	  genTower_con_iphi   = (int) towerConHcalId.iphi();
	  genTower_con_zside  = (int) towerConHcalId.zside();
	  genTower_con_depth  = (int) towerConHcalId.depth();

	  m_myTree.genTower_con_subdet[ngtower][iTowerCon] = genTower_con_subdet;
	  m_myTree.genTower_con_ieta  [ngtower][iTowerCon] = genTower_con_ieta;
	  m_myTree.genTower_con_iphi  [ngtower][iTowerCon] = genTower_con_iphi;
	  m_myTree.genTower_con_zside [ngtower][iTowerCon] = genTower_con_zside;
	  m_myTree.genTower_con_depth [ngtower][iTowerCon] = genTower_con_depth;

	  ngtower_con++;
	  
	}	  
      }
	

      m_myTree.genTower_ncon[ngtower] = ngtower_con;

      ngtower++;
      
    }
    
    m_myTree.ngtower = ngtower;
  }
  
  else
    printf("NO Calo Towers!\n");
    
  //-----------------------------------------------------
  // Loop over jets and get towers, then tower digis from
  // them.
  //-----------------------------------------------------

  if (jets.isValid()){
    
    m_myTree.run   = run;
    m_myTree.event = event;
    
    //-----------------------------------------------------
    // Loop over the reco jets
    //-----------------------------------------------------
    
    for (CaloJetCollection::const_iterator iJet = jets->begin();
	 iJet != jets->end();
	 iJet ++){

      p   = (float) (*iJet).p();
      px  = (float) (*iJet).px();
      py  = (float) (*iJet).py();
      pz  = (float) (*iJet).pz();
      pt  = (float) (*iJet).pt();
      
      e   = (float) (*iJet).energy();
      et  = (float) (*iJet).et();
      
      phi = (float) (*iJet).phi();
      eta = (float) (*iJet).eta();
      
      m_myTree.jet_p [nrjet] = p;
      m_myTree.jet_pt[nrjet] = pt;
      m_myTree.jet_px[nrjet] = px;
      m_myTree.jet_py[nrjet] = py;
      m_myTree.jet_pz[nrjet] = pz;
      
      m_myTree.jet_e [nrjet] = e;
      m_myTree.jet_et[nrjet] = et;
      
      m_myTree.jet_eta[nrjet] = eta;
      m_myTree.jet_phi[nrjet] = phi;
      
      //-----------------------------------------------------
      // Loop over the towers for each reco jet
      // Minor changes have been made to run this in CMSSW_2_2_X
      // (see commented lines)
      //-----------------------------------------------------
      
      //std::vector <CaloTowerRef> jetTowers = (*iJet).getConstituents();
      std::vector <CaloTowerPtr> jetTowers = (*iJet).getCaloConstituents();
      
      njtower = 0;
      
      
      //for(std::vector<CaloTowerRef>::iterator iJetTower = jetTowers.begin();
      for(std::vector<CaloTowerPtr>::iterator iJetTower = jetTowers.begin();
	  iJetTower != jetTowers.end();
	  iJetTower++){	   

	ieta      = (int)   (**iJetTower).id().ieta();
	ietaAbs   = (int)   (**iJetTower).id().ietaAbs();
	iphi      = (int)   (**iJetTower).id().iphi();
	zside     = (int)   (**iJetTower).id().zside();
	
	eta       = (float) (**iJetTower).eta();
	phi       = (float) (**iJetTower).phi();
	
	hadEnergy = (float) (**iJetTower).hadEnergy();
	hadEt     = (float) (**iJetTower).hadEt();
	emEnergy  = (float) (**iJetTower).emEnergy();
	emEt      = (float) (**iJetTower).emEt();
	
	m_myTree.jetTower_ieta     [nrjet][njtower] = ieta;
	m_myTree.jetTower_ietaAbs  [nrjet][njtower] = ietaAbs;
	m_myTree.jetTower_iphi     [nrjet][njtower] = iphi;
	m_myTree.jetTower_zside    [nrjet][njtower] = zside;
	
	m_myTree.jetTower_eta      [nrjet][njtower] = eta;
	m_myTree.jetTower_phi      [nrjet][njtower] = phi;
	m_myTree.jetTower_emEnergy [nrjet][njtower] = emEnergy;
	m_myTree.jetTower_emEt     [nrjet][njtower] = emEt;
	m_myTree.jetTower_hadEnergy[nrjet][njtower] = hadEnergy;
	m_myTree.jetTower_hadEt    [nrjet][njtower] = hadEt;
	
	//-----------------------------------------------------
	// Loop over the hits for each tower
	//-----------------------------------------------------
	
	int nRecHits = (int) (**iJetTower).constituentsSize();
	
	nrhits = 0;

	for (int iRecHit = 0; iRecHit < nRecHits; iRecHit++){
	  
	  DetId RecHitDetID = (**iJetTower).constituent(nrhits);
	  DetId::Detector DetNum = RecHitDetID.det();

	  //-----------------------------------------------------
	  // Are we in the HCAL? (not trivial)
	  //-----------------------------------------------------

	  if (DetNum == DetId::Hcal){
	    
	    HcalDetId HcalID = RecHitDetID;
	    
	    const HcalQIECoder* channelCoder = conditions->getHcalCoder(HcalID);
	    HcalCoderDb coder (*channelCoder, *shape);

	    HcalSubdetector HcalSubdet = HcalID.subdet();

	    subdet = (int) HcalSubdet;
	    depth  = (int) HcalID.depth();

	    m_myTree.jetTower_rhSubdet[nrjet][njtower][nrhits] = (int) subdet;
	    m_myTree.jetTower_depth   [nrjet][njtower][nrhits] = (int) depth;
	    
	    //-----------------------------------------------------
	    // Decode the Digi, depending on the subdetector
	    //-----------------------------------------------------

	    ntslice = 0;

	    if (HcalSubdet == HcalBarrel || HcalSubdet == HcalEndcap){

	      HBHEDigiCollection::const_iterator   iHBHEDigi   = HBHEDigis   -> find(HcalID);
	      HBHERecHitCollection::const_iterator iHBHERecHit = HBHERecHits -> find(HcalID);

	      coder.adc2fC(*iHBHEDigi,tool);

	      for (int iTS = 0; iTS < tool.size(); iTS++){
		fC_fillVal  = (float) tool[iTS];
		ADC_fillVal = (int) (*iHBHEDigi)[iTS].adc();

		m_myTree.jetTower_digi_fC  [nrjet][njtower][nrhits][iTS] = fC_fillVal;
		m_myTree.jetTower_digi_ADC [nrjet][njtower][nrhits][iTS] = ADC_fillVal;
		
		ntslice++;
	      }

	      m_myTree.jetTower_rhEnergy[nrjet][njtower][nrhits] = (*iHBHERecHit).energy();

	    }
	    
	    else if (HcalSubdet == HcalOuter){

	      HODigiCollection::const_iterator   iHODigi   = HODigis   -> find(HcalID);
	      HORecHitCollection::const_iterator iHORecHit = HORecHits -> find(HcalID);

	      coder.adc2fC(*iHODigi,tool);

	      for (int iTS = 0; iTS < tool.size(); iTS++){
		fC_fillVal  = (float) tool[iTS];
		ADC_fillVal = (int) (*iHODigi)[iTS].adc();

		m_myTree.jetTower_digi_fC  [nrjet][njtower][nrhits][iTS] = fC_fillVal;
		m_myTree.jetTower_digi_ADC [nrjet][njtower][nrhits][iTS] = ADC_fillVal;

		ntslice++;

	      }

	      m_myTree.jetTower_rhEnergy[nrjet][njtower][nrhits] = (*iHORecHit).energy();

	    }
	    
	    else if (HcalSubdet == HcalForward){

	      HFDigiCollection::const_iterator   iHFDigi   = HFDigis   -> find(HcalID);
	      HFRecHitCollection::const_iterator iHFRecHit = HFRecHits -> find(HcalID);

	      coder.adc2fC(*iHFDigi,tool);

	      for (int iTS = 0; iTS < tool.size(); iTS++){
		fC_fillVal  = (float) tool[iTS];
		ADC_fillVal = (int) (*iHFDigi)[iTS].adc();

		m_myTree.jetTower_digi_fC  [nrjet][njtower][nrhits][iTS] = fC_fillVal;
		m_myTree.jetTower_digi_ADC [nrjet][njtower][nrhits][iTS] = ADC_fillVal;

		ntslice++;
	      }
	      
	      m_myTree.jetTower_rhEnergy[nrjet][njtower][nrhits] = (*iHFRecHit).energy();

	    }

	    else {
	      printf ("Weird subdet!\n");
	      exit(0);
	    } 
	    
	    m_myTree.jetTower_ntslice[nrjet][njtower][nrhits] = (int) ntslice;	    
	    
	    nrhits++;
	  }  
	  
	}	
	m_myTree.jetTower_nrhits[nrjet][njtower] = (int) nrhits;
	njtower++;	
      }      
      m_myTree.jet_njtowers[nrjet] = (int) njtower;
      nrjet++;      
    }
    m_myTree.nrjet = (int) nrjet;
  }
  
  else printf("NO jets!\n");
  
  //-----------------------------------------------------
  // Fill the tree and finish
  //-----------------------------------------------------

  m_fillTree.fill();
  
}

void DijetAnalyzer::beginJob(const edm::EventSetup&){
  
}

void DijetAnalyzer::endJob() {
  
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  m_fillTree.finalize();
}

//define this as a plug-in
DEFINE_FWK_MODULE(DijetAnalyzer);
