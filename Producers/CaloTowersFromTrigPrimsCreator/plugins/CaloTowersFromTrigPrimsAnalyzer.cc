#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAnalyzer.h"

//---------------------------------------------
// Constructor
//---------------------------------------------

CaloTowersFromTrigPrimsAnalyzer::CaloTowersFromTrigPrimsAnalyzer(const edm::ParameterSet& iConfig){

  //---------------------------------------------
  // Get user input
  //---------------------------------------------

  edm::InputTag d_defaultCaloTowerTag ("towerMaker");
  edm::InputTag d_createdCaloTowerTag ("caloTowersFromTrigPrimsCreator");
  edm::InputTag d_hcalTrigPrimTag("simHcalTriggerPrimitiveDigis");
  edm::InputTag d_ecalTrigPrimTag("simEcalTriggerPrimitiveDigis");

  m_hcalTrigPrimTag     = iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigPrimTag",d_hcalTrigPrimTag);
  m_ecalTrigPrimTag     = iConfig.getUntrackedParameter<edm::InputTag>("ecalTrigPrimTag",d_ecalTrigPrimTag);
  m_createdCaloTowerTag = iConfig.getUntrackedParameter<edm::InputTag>("createdCaloTowerTag",d_createdCaloTowerTag);
  m_defaultCaloTowerTag = iConfig.getUntrackedParameter<edm::InputTag>("defaultCaloTowerTag",d_defaultCaloTowerTag);

  std::string d_outputFileName("CaloTowersFromTrigPrimsAnalyzerOutput.root");
  std::string outputFileName = iConfig.getUntrackedParameter<std::string>("outputFileName", d_outputFileName);

  //---------------------------------------------
  // Initialize ROOT trees
  //---------------------------------------------

  m_fillCaloTowersFromTrigPrimsAnalyzerTree.init(outputFileName,&m_caloTowersFromTrigPrimsAnalyzerTree);

}

//---------------------------------------------
// Destructor
//---------------------------------------------

CaloTowersFromTrigPrimsAnalyzer::~CaloTowersFromTrigPrimsAnalyzer() {}

//---------------------------------------------
// Methods for getting trig tower E_T
//---------------------------------------------

float CaloTowersFromTrigPrimsAnalyzer::getTrigTowerET ( const HcalTriggerPrimitiveDigi& hcalTrigPrimDigi ){

  unsigned short hcalRctInput = hcalTrigPrimDigi.SOI_compressedEt() * 2 + hcalTrigPrimDigi.SOI_fineGrain();
  hcalRctInput /= 2;
  
  unsigned short ietaAbs = hcalTrigPrimDigi.id().ietaAbs();
  short sign             = hcalTrigPrimDigi.id().zside();

  double et = (double) m_l1CaloHcalScale -> et (hcalRctInput, ietaAbs, sign);

  return et;
}


float CaloTowersFromTrigPrimsAnalyzer::getTrigTowerET ( const EcalTriggerPrimitiveDigi& ecalTrigPrimDigi ){

  unsigned short ecalRctInput = ecalTrigPrimDigi.compressedEt() * 2 + ecalTrigPrimDigi.fineGrain();
  ecalRctInput /= 2;
  
  unsigned short ietaAbs = ecalTrigPrimDigi.id().ietaAbs();
  short sign             = ecalTrigPrimDigi.id().zside();
  
  double et = (double) m_l1CaloEcalScale -> et (ecalRctInput, ietaAbs, sign);
  
  return et;

}

//---------------------------------------------
// Methods for getting trig tower position
//---------------------------------------------

float CaloTowersFromTrigPrimsAnalyzer::getTrigTowerMeanEta(const HcalTriggerPrimitiveDigi& hcalTrigPrimDigi){

  double etaMin, etaMax, etaMean;

  m_hcalTrigTowerGeometry.towerEtaBounds( hcalTrigPrimDigi.id().ieta(), etaMin, etaMax);

  etaMean = 0.5 * ( etaMin + etaMax );
  
  return etaMean;
  
}


float CaloTowersFromTrigPrimsAnalyzer::getTrigTowerMeanEta(const EcalTriggerPrimitiveDigi& ecalTrigPrimDigi){
  
  double meanTheta = 0.0;
  
  std::vector<DetId> EcalDetIds = m_ecalTrigTowerConstituentsMap -> constituentsOf(ecalTrigPrimDigi.id());
  std::vector<DetId>::iterator ecalDetId = EcalDetIds.begin();
  
  for (; ecalDetId != EcalDetIds.end(); ecalDetId++){

    if ((*ecalDetId).subdetId() == EcalBarrel) 
      meanTheta += (double) m_ecalBarrelGeometry -> getGeometry( *ecalDetId ) -> getPosition().theta();
    if ((*ecalDetId).subdetId() == EcalEndcap) 
      meanTheta += (double) m_ecalEndcapGeometry -> getGeometry( *ecalDetId ) -> getPosition().theta();

  }

  if (EcalDetIds.size() != 0 ) meanTheta /= (double) EcalDetIds.size();
  if (EcalDetIds.size() == 0 ) return -999.0;
  
  double meanEta = (-1.0) * TMath::Log( TMath::Tan ( (meanTheta / 2.0) ) );

  return meanEta;

}

//---------------------------------------------
// Main analysis function
//---------------------------------------------

void CaloTowersFromTrigPrimsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  
  //-----------------------------------------------------
  // Get all of your ESHandles
  //-----------------------------------------------------
  
  iSetup.get<IdealGeometryRecord>     ().get(m_ecalTrigTowerConstituentsMap    );
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel", m_ecalBarrelGeometry);
  iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap", m_ecalEndcapGeometry);
  iSetup.get<L1CaloEcalScaleRcd>      ().get(m_l1CaloEcalScale                 );
  iSetup.get<L1CaloHcalScaleRcd>      ().get(m_l1CaloHcalScale                 );
  
  //-----------------------------------------------------
  // Pass on necessary pointers
  //-----------------------------------------------------

  m_event = &iEvent;
  m_setup = &iSetup;
  
  //-----------------------------------------------------
  // Initialize the ROOT trees
  //-----------------------------------------------------
  
  m_caloTowersFromTrigPrimsAnalyzerTree.init();
  
  //-----------------------------------------------------
  // Get run and event information
  //-----------------------------------------------------
  
  int run   = (int) m_event -> id().run();
  int event = (int) m_event -> id().event();
  
  m_caloTowersFromTrigPrimsAnalyzerTree.run   = run;
  m_caloTowersFromTrigPrimsAnalyzerTree.event = event;
  
  //-----------------------------------------------------
  // Perform the analysis
  //-----------------------------------------------------

  analyzeTPGs();
  analyzeCaloTowers();

  //-----------------------------------------------------
  // Finalize the ROOT trees
  //-----------------------------------------------------

  m_fillCaloTowersFromTrigPrimsAnalyzerTree.fill();

}

//-----------------------------------------------------
// Get TPG handles and analyze
//-----------------------------------------------------

void CaloTowersFromTrigPrimsAnalyzer::analyzeTPGs(){
  
  //-----------------------------------------------------
  // Get trigger primitives from the event
  //-----------------------------------------------------
  
  edm::Handle<HcalTrigPrimDigiCollection> HCALTrigPrimDigis;
  bool hcalTrigPrimDigiTagExists = m_event -> getByLabel(m_hcalTrigPrimTag,HCALTrigPrimDigis);
  if (!hcalTrigPrimDigiTagExists){
    edm::LogWarning("TPGAnalyzer") << "Could not extract HCAL trigger primitives with " << m_hcalTrigPrimTag;
    return;
  }

  edm::Handle<EcalTrigPrimDigiCollection> ECALTrigPrimDigis;
  bool ecalTrigPrimDigiTagExists = m_event -> getByLabel(m_ecalTrigPrimTag,ECALTrigPrimDigis);
  if (!ecalTrigPrimDigiTagExists){
    edm::LogWarning("TPGAnalyzer") << "Could not extract ECAL trigger primitives with " << m_ecalTrigPrimTag;
    return;
  }

  //-----------------------------------------------------
  // Analyze the individual handles
  //-----------------------------------------------------

  analyzeTPGs(*HCALTrigPrimDigis);
  analyzeTPGs(*ECALTrigPrimDigis);

}

template <typename TrigPrimDigiCollection>
void CaloTowersFromTrigPrimsAnalyzer::analyzeTPGs(const TrigPrimDigiCollection& trigPrimDigiCollection){

  typedef typename TrigPrimDigiCollection::const_iterator Iterator;
  
  Iterator trigPrimDigi = trigPrimDigiCollection.begin();

  int nTPG = m_caloTowersFromTrigPrimsAnalyzerTree.ntpg;
  if ( nTPG == -999 ) nTPG = 0;

  for(; trigPrimDigi != trigPrimDigiCollection.end(); ++trigPrimDigi)    {

    int ietaAbs = (*trigPrimDigi).id().ietaAbs();

    int ieta = (*trigPrimDigi).id().ieta();
    int iphi = (*trigPrimDigi).id().iphi();
    
    int isEcal = 0;
    int isHcal = 0;
    int isHF   = 0;
    
    if ( (*trigPrimDigi).id().det() == DetId::Ecal ) isEcal = 1;
    if ( (*trigPrimDigi).id().det() == DetId::Hcal ) isHcal = 1;    
    if ( isHcal == 1 && ietaAbs >= 29              ) isHF   = 1;

    float et      = getTrigTowerET     (*trigPrimDigi);
    float etaMean = getTrigTowerMeanEta(*trigPrimDigi);
    float energy  = et * cosh(etaMean);
    
    m_caloTowersFromTrigPrimsAnalyzerTree.tpg_ieta   [nTPG] = ieta;
    m_caloTowersFromTrigPrimsAnalyzerTree.tpg_iphi   [nTPG] = iphi;

    m_caloTowersFromTrigPrimsAnalyzerTree.tpg_isHcal [nTPG] = isHcal;
    m_caloTowersFromTrigPrimsAnalyzerTree.tpg_isEcal [nTPG] = isEcal;
    m_caloTowersFromTrigPrimsAnalyzerTree.tpg_isHF   [nTPG] = isHF  ;

    m_caloTowersFromTrigPrimsAnalyzerTree.tpg_et     [nTPG] = et;
    m_caloTowersFromTrigPrimsAnalyzerTree.tpg_energy [nTPG] = energy;
    m_caloTowersFromTrigPrimsAnalyzerTree.tpg_meanEta[nTPG] = etaMean;
    
    nTPG++;
  }

  m_caloTowersFromTrigPrimsAnalyzerTree.ntpg = nTPG;

}

//-----------------------------------------------------
// Get calotower handles and analyze
//-----------------------------------------------------

void CaloTowersFromTrigPrimsAnalyzer::analyzeCaloTowers(){

  //---------------------------------------------
  // Get the 'view' of CaloTowers from Candidate
  //---------------------------------------------

  edm::Handle<edm::View<reco::Candidate> > created_towers;
  bool gotHandle =  m_event -> getByLabel( m_createdCaloTowerTag , created_towers);
  if (!gotHandle) {
    edm::LogWarning("CaloTowersFromTrigPrimsAnalyzer") << "Could not extract created calo towers with tag: " << m_createdCaloTowerTag;
    return;
  }


  edm::Handle<edm::View<reco::Candidate> > default_towers;
  gotHandle =  m_event -> getByLabel( m_defaultCaloTowerTag , default_towers);
  if (!gotHandle) {
    edm::LogWarning("CaloTowersFromTrigPrimsAnalyzer") << "Could not extract default calo towers with tag: " << m_defaultCaloTowerTag;
    return;
  }

  bool isMine = true;

  analyzeCaloTowers (*created_towers, isMine);
  analyzeCaloTowers (*default_towers,!isMine);

}

void CaloTowersFromTrigPrimsAnalyzer::analyzeCaloTowers( const edm::View<reco::Candidate> & towers, bool isMine){

  using namespace reco;
  using namespace std;
  using namespace edm;
  
  //---------------------------------------------
  // Loop over CaloTowers
  //---------------------------------------------

  int   nct = m_caloTowersFromTrigPrimsAnalyzerTree.nct;
  if ( nct == -999 ) nct = 0;

  int   ct_ieta, ct_iphi, ct_nhcon, ct_necon;
  float ct_eEm, ct_eHad, ct_eOut;
  float ct_etEm, ct_etHad, ct_etOut;
  
  for (unsigned int iTower = 0; iTower < towers.size(); iTower++){
    
    CaloTowerRef calotower = (towers.refAt(iTower)).castTo<CaloTowerRef>();
    if (!calotower.isNull()){
      
      ct_ieta  = calotower->id().ieta();
      ct_iphi  = calotower->id().iphi();
      
      ct_eEm   = (float) calotower->emEnergy();
      ct_eHad  = (float) calotower->hadEnergy();
      ct_eOut  = (float) calotower->outerEnergy();

      ct_etEm  = (float) calotower->emEt();
      ct_etHad = (float) calotower->hadEt();
      ct_etOut = (float) calotower->outerEt();

      vector<DetId> constituents = calotower->constituents();
      vector<DetId>::iterator constituent = constituents.begin();

      ct_nhcon = 0;
      ct_necon = 0;

      for (; constituent != constituents.end(); constituent++){

	if ((*constituent).det() == DetId::Hcal) ct_nhcon++;
	if ((*constituent).det() == DetId::Ecal) ct_necon++;

      }

      m_caloTowersFromTrigPrimsAnalyzerTree.ct_isMine [nct] = 1;
      if (!isMine) m_caloTowersFromTrigPrimsAnalyzerTree.ct_isMine [nct] = 0;
 						      
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_ieta   [nct] = ct_ieta;
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_iphi   [nct] = ct_iphi;
		    				      
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_nhcon  [nct] = ct_nhcon;
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_necon  [nct] = ct_necon;
		    				      
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_eEm    [nct] = ct_eEm ;
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_eHad   [nct] = ct_eHad;
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_eOut   [nct] = ct_eOut;
      						      
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_etEm   [nct] = ct_etEm ;
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_etHad  [nct] = ct_etHad;
      m_caloTowersFromTrigPrimsAnalyzerTree.ct_etOut  [nct] = ct_etOut;

      nct++;
      
    } // end if(!calotower.isNull())    

  } // end loop over candidates


  m_caloTowersFromTrigPrimsAnalyzerTree.nct = nct;
}

//-----------------------------------------------------
// Required analysis methods
//-----------------------------------------------------

void CaloTowersFromTrigPrimsAnalyzer::beginJob(const edm::EventSetup&){}

void CaloTowersFromTrigPrimsAnalyzer::endJob() {
  
  m_fillCaloTowersFromTrigPrimsAnalyzerTree.finalize();
  
}


DEFINE_FWK_MODULE(CaloTowersFromTrigPrimsAnalyzer);
