#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsCreator.h"

//------------------------------------------------------
// Constructor, 
//   get tags for data collections
//   and instruction bools here.
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::CaloTowersFromTrigPrimsCreator(const ParameterSet& iConfig) :
m_caloTowersFromTrigPrimsAlgo(){

  InputTag d_hcalTrigPrimTag("simHcalTriggerPrimitiveDigis");
  m_hcalTrigPrimTag = iConfig.getUntrackedParameter<InputTag>("hcalTrigPrimTag",d_hcalTrigPrimTag);

  InputTag d_ecalTrigPrimTag("simEcalTriggerPrimitiveDigis");
  m_ecalTrigPrimTag = iConfig.getUntrackedParameter<InputTag>("ecalTrigPrimTag",d_ecalTrigPrimTag);
  
  InputTag d_hoTrigPrimTag("");
  m_hoTrigPrimTag = iConfig.getUntrackedParameter<InputTag>("hoTrigPrimTag",d_hoTrigPrimTag);

  bool d_verbose = true;
  m_verbose = iConfig.getUntrackedParameter("verbose",d_verbose);

  produces<CaloTowerCollection>();

}

//------------------------------------------------------
// Destructor
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::~CaloTowersFromTrigPrimsCreator(){}

//------------------------------------------------------
// Main production function
//------------------------------------------------------

void CaloTowersFromTrigPrimsCreator::produce(Event& iEvent, const EventSetup& iSetup) {

  //-----------------------------------------------------
  // Get all of your ESHandles
  //-----------------------------------------------------

  iSetup.get<CaloGeometryRecord>      ().get(m_geometry                        );
  iSetup.get<IdealGeometryRecord>     ().get(m_caloTowerConstituentsMap        );
  iSetup.get<IdealGeometryRecord>     ().get(m_ecalTrigTowerConstituentsMap    );
  iSetup.get<CaloTPGRecord>           ().get(m_caloTPGTranscoder               );
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel", m_ecalBarrelGeometry);
  iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap", m_ecalEndcapGeometry);

  m_ecalTPGScale.setEventSetup(iSetup);

  //-----------------------------------------------------
  // Get trigger primitives from the event
  //-----------------------------------------------------
  
  Handle<HcalTrigPrimDigiCollection> HCALTrigPrimDigis;
  bool hcalTrigPrimDigiTagExists = iEvent.getByLabel(m_hcalTrigPrimTag,HCALTrigPrimDigis);
  if (!hcalTrigPrimDigiTagExists){
    LogWarning("CaloTowersFromTrigPrimsCreator") << "Could not extract HCAL trigger primitives with " << m_hcalTrigPrimTag;
    return;
  }

  Handle<EcalTrigPrimDigiCollection> ECALTrigPrimDigis;
  bool ecalTrigPrimDigiTagExists = iEvent.getByLabel(m_ecalTrigPrimTag,ECALTrigPrimDigis);
  if (!ecalTrigPrimDigiTagExists){
    LogWarning("CaloTowersFromTrigPrimsCreator") << "Could not extract ECAL trigger primitives with " << m_ecalTrigPrimTag;
    return;
  }
  
  //------------------------------------------------------
  // Create an empty collection
  //------------------------------------------------------
  
  auto_ptr<CaloTowerCollection> caloTowerCollection(new CaloTowerCollection());

  //-----------------------------------------------------
  // Initialize the algorithm
  //-----------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo.setGeometry(m_geometry.product(),
					    m_caloTowerConstituentsMap.product(), 
					    m_ecalTrigTowerConstituentsMap.product(),
					    m_ecalBarrelGeometry.product(),
					    m_ecalEndcapGeometry.product() 	       );
  
  m_caloTowersFromTrigPrimsAlgo.setHcalTPGCoder (m_caloTPGTranscoder.product());
  m_caloTowersFromTrigPrimsAlgo.setEcalTPGScale (m_ecalTPGScale);
  m_caloTowersFromTrigPrimsAlgo.setMapping();


  //------------------------------------------------------
  // Apply the algorithm
  //------------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo.process(*HCALTrigPrimDigis);
  m_caloTowersFromTrigPrimsAlgo.process(*ECALTrigPrimDigis);

  //------------------------------------------------------
  // Finish the algorithm
  //------------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo.finish(*caloTowerCollection);

  //------------------------------------------------------
  // Add the final CaloTowerCollection to the event
  //------------------------------------------------------

  iEvent.put(caloTowerCollection);

}

void CaloTowersFromTrigPrimsCreator::beginJob(const EventSetup&){}

void CaloTowersFromTrigPrimsCreator::endJob() {}

DEFINE_FWK_MODULE(CaloTowersFromTrigPrimsCreator);
