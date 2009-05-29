// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Data format files
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

// Root tree info
#include "Analyzers/CaloTowerComp/interface/FillCaloTowerCompTree.h"
#include "Analyzers/CaloTowerComp/interface/CaloTowerCompTree.h"

class CaloTowerComp : public edm::EDAnalyzer {
public:
  explicit CaloTowerComp(const edm::ParameterSet&);
  ~CaloTowerComp();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  void analyzeCaloTowers( edm::InputTag tag, CaloTowerCompTree& tree );

  //---------------------------------------------
  // Event and EventSetup pointers
  //---------------------------------------------
  
  const edm::Event*      m_event;
  const edm::EventSetup* m_setup;

  //---------------------------------------------
  // edm::InputTag entries
  //---------------------------------------------
  
  edm::InputTag m_defaultCaloTowerTag;
  edm::InputTag m_createdCaloTowerTag;

  //---------------------------------------------
  // ROOT tree objects
  //---------------------------------------------

  FillCaloTowerCompTree m_fillCreatedCaloTowerCompTree;
  FillCaloTowerCompTree m_fillDefaultCaloTowerCompTree;
  
  CaloTowerCompTree     m_createdCaloTowerCompTree;
  CaloTowerCompTree     m_defaultCaloTowerCompTree;
  

};

CaloTowerComp::CaloTowerComp(const edm::ParameterSet& iConfig){

  //---------------------------------------------
  // Get user input
  //---------------------------------------------

  const edm::InputTag d_defaultCaloTowerTag("towerMaker");		      
  const edm::InputTag d_createdCaloTowerTag("caloTowersFromTrigPrimsCreator");

  m_defaultCaloTowerTag = iConfig.getUntrackedParameter<edm::InputTag>("defaultCaloTowerTag",d_defaultCaloTowerTag);
  m_createdCaloTowerTag = iConfig.getUntrackedParameter<edm::InputTag>("createdCaloTowerTag",d_createdCaloTowerTag);

  const std::string d_defaultFileName("data/CaloTowerInfo_Default.root");
  const std::string d_createdFileName("data/CaloTowerInfo_Created.root");

  std::string defaultCaloTowerFileName = iConfig.getUntrackedParameter<std::string>("defaultCaloTowerFileName", d_defaultFileName);
  std::string createdCaloTowerFileName = iConfig.getUntrackedParameter<std::string>("createdCaloTowerFileName", d_createdFileName);

  //---------------------------------------------
  // Initialize the ROOT trees
  //---------------------------------------------
  
  m_fillCreatedCaloTowerCompTree.init(createdCaloTowerFileName , &m_createdCaloTowerCompTree);
  m_fillDefaultCaloTowerCompTree.init(defaultCaloTowerFileName , &m_defaultCaloTowerCompTree);
  
}


CaloTowerComp::~CaloTowerComp(){}

void CaloTowerComp::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  using namespace edm;
  
  m_createdCaloTowerCompTree.init();
  m_defaultCaloTowerCompTree.init();

  //-----------------------------------------------------
  // Pass on necessary pointers
  //-----------------------------------------------------

  m_event = &iEvent;
  m_setup = &iSetup;

  //-----------------------------------------------------
  // Get run and event information
  //-----------------------------------------------------
  
  int run   = (int) m_event -> id().run();
  int event = (int) m_event -> id().event();

  m_createdCaloTowerCompTree.run   = run;
  m_defaultCaloTowerCompTree.run   = run;

  m_createdCaloTowerCompTree.event = event;
  m_defaultCaloTowerCompTree.event = event;

  //-----------------------------------------------------
  // Run helper analysis functions
  //-----------------------------------------------------

  analyzeCaloTowers( m_defaultCaloTowerTag, m_defaultCaloTowerCompTree );
  analyzeCaloTowers( m_createdCaloTowerTag, m_createdCaloTowerCompTree );

  //-----------------------------------------------------
  // Fill ROOT tree
  //-----------------------------------------------------
  
  m_fillCreatedCaloTowerCompTree.fill();
  m_fillDefaultCaloTowerCompTree.fill();
  
}

void CaloTowerComp::analyzeCaloTowers( edm::InputTag tag, CaloTowerCompTree& tree ){

  using namespace reco;
  using namespace std;
  using namespace edm;
  
  //---------------------------------------------
  // Get the 'view' of CaloTowers from Candidate
  //---------------------------------------------

  edm::Handle<edm::View<Candidate> > towers;
  bool gotHandle =  m_event -> getByLabel(tag ,towers);
  if (!gotHandle) {
    LogWarning("CaloTowerComp") << "Could not extract calo towers with tag: " << m_createdCaloTowerTag << endl;
    return;
  }

  //---------------------------------------------
  // Loop over CaloTowers
  //---------------------------------------------

  int   nct = 0;

  int   ct_ieta, ct_iphi, ct_nhcon, ct_necon;
  float ct_eEm, ct_eHad, ct_eOut;
  
  for (unsigned int iTower = 0; iTower < towers -> size(); iTower++){
    
    CaloTowerRef calotower = ((*towers).refAt(iTower)).castTo<CaloTowerRef>();
    if (!calotower.isNull()){
      
      ct_ieta = calotower->id().ieta();
      ct_iphi = calotower->id().iphi();

      ct_eEm  = calotower->emEnergy();
      ct_eHad = calotower->hadEnergy();
      ct_eOut = calotower->outerEnergy();

      vector<DetId> constituents = calotower->constituents();
      vector<DetId>::iterator constituent = constituents.begin();

      ct_nhcon = 0;
      ct_necon = 0;

      for (; constituent != constituents.end(); constituent++){

	if ((*constituent).det() == DetId::Hcal) ct_nhcon++;
	if ((*constituent).det() == DetId::Ecal) ct_necon++;

      }

      tree.ct_ieta [nct] = ct_ieta;
      tree.ct_iphi [nct] = ct_iphi;

      tree.ct_nhcon[nct] = ct_nhcon;
      tree.ct_necon[nct] = ct_necon;

      tree.ct_eEm  [nct] = ct_eEm ;
      tree.ct_eHad [nct] = ct_eHad;
      tree.ct_eOut [nct] = ct_eOut;
      
      nct++;
      
    } // end if(calotower)    

    tree.nct = nct;

  } // end loop over candidates
}

void CaloTowerComp::beginJob(const edm::EventSetup&){}

void CaloTowerComp::endJob(){

  m_fillCreatedCaloTowerCompTree.finalize();
  m_fillDefaultCaloTowerCompTree.finalize();

}


DEFINE_FWK_MODULE(CaloTowerComp);
