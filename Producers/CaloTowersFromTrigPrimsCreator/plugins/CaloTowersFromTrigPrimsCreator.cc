#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsCreator.h"

//------------------------------------------------------
// Constructor, 
//   get tags for data collections
//   and instruction bools here.
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::CaloTowersFromTrigPrimsCreator(const edm::ParameterSet& iConfig) :
  m_caloTowersFromTrigPrimsAlgo(){
  
  //-----------------------------------------------------
  // Get input from the user
  //-----------------------------------------------------
  
  double d_momHBDepth = 0.2;
  double d_momHEDepth = 0.4;
  double d_momEBDepth = 0.3;
  double d_momEEDepth = 0.0;
  m_momHBDepth = iConfig.getUntrackedParameter<double> ("momHBDepth",d_momHBDepth);  
  m_momHEDepth = iConfig.getUntrackedParameter<double> ("momHEDepth",d_momHEDepth);
  m_momEBDepth = iConfig.getUntrackedParameter<double> ("momEBDepth",d_momEBDepth);  
  m_momEEDepth = iConfig.getUntrackedParameter<double> ("momEEDepth",d_momEEDepth);  

  edm::InputTag d_hcalTrigPrimTag("simHcalTriggerPrimitiveDigis");
  edm::InputTag d_ecalTrigPrimTag("simEcalTriggerPrimitiveDigis");
  m_hcalTrigPrimTag = iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigPrimTag",d_hcalTrigPrimTag);
  m_ecalTrigPrimTag = iConfig.getUntrackedParameter<edm::InputTag>("ecalTrigPrimTag",d_ecalTrigPrimTag);
  
  bool d_verbose = true;
  m_verbose = iConfig.getUntrackedParameter<bool>("verbose",d_verbose);

  //-----------------------------------------------------
  // Tell producer what to produce
  //-----------------------------------------------------

  produces<CaloTowerCollection>();

}

//------------------------------------------------------
// Destructor
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::~CaloTowersFromTrigPrimsCreator(){}

//------------------------------------------------------
// Main production function
//------------------------------------------------------

void CaloTowersFromTrigPrimsCreator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

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
  
  edm::Handle<HcalTrigPrimDigiCollection> HCALTrigPrimDigis;
  bool hcalTrigPrimDigiTagExists = iEvent.getByLabel(m_hcalTrigPrimTag,HCALTrigPrimDigis);
  if (!hcalTrigPrimDigiTagExists){
    edm::LogWarning("CaloTowersFromTrigPrimsCreator") << "Could not extract HCAL trigger primitives with " << m_hcalTrigPrimTag;
    return;
  }

  edm::Handle<EcalTrigPrimDigiCollection> ECALTrigPrimDigis;
  bool ecalTrigPrimDigiTagExists = iEvent.getByLabel(m_ecalTrigPrimTag,ECALTrigPrimDigis);
  if (!ecalTrigPrimDigiTagExists){
    edm::LogWarning("CaloTowersFromTrigPrimsCreator") << "Could not extract ECAL trigger primitives with " << m_ecalTrigPrimTag;
    return;
  }
  
  //------------------------------------------------------
  // Create an empty collection
  //------------------------------------------------------
  
  std::auto_ptr<CaloTowerCollection> caloTowerCollection(new CaloTowerCollection());

  //-----------------------------------------------------
  // Initialize the algorithm
  //-----------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo.setVerbose       (m_verbose);
  m_caloTowersFromTrigPrimsAlgo.setGeometry      (m_geometry.product(),
						  m_caloTowerConstituentsMap.product(), 
						  m_ecalTrigTowerConstituentsMap.product(),
						  m_ecalBarrelGeometry.product(),
						  m_ecalEndcapGeometry.product() 	    );
  
  m_caloTowersFromTrigPrimsAlgo.setDetectorDepths(m_momHBDepth, m_momHEDepth,
						  m_momEBDepth, m_momEEDepth );
  
  m_caloTowersFromTrigPrimsAlgo.setHcalTPGCoder  (m_caloTPGTranscoder.product());
  m_caloTowersFromTrigPrimsAlgo.setEcalTPGScale  (m_ecalTPGScale);

  //m_caloTowersFromTrigPrimsAlgo.setMapping();
  
  //------------------------------------------------------
  // Apply the algorithm
  //------------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo.process(*ECALTrigPrimDigis);
  m_caloTowersFromTrigPrimsAlgo.process(*HCALTrigPrimDigis);
  
  //------------------------------------------------------
  // Finish the algorithm
  //------------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo.finish(*caloTowerCollection);

  //------------------------------------------------------
  // Add the final CaloTowerCollection to the event
  //------------------------------------------------------

  iEvent.put(caloTowerCollection);

}

void CaloTowersFromTrigPrimsCreator::beginJob(const edm::EventSetup&){}

void CaloTowersFromTrigPrimsCreator::endJob() {}

DEFINE_FWK_MODULE(CaloTowersFromTrigPrimsCreator);
