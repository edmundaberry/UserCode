#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsCreator.h"
#include <time.h>

//------------------------------------------------------
// Constructor, 
//   get tags for data collections
//   and instruction bools here.
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::CaloTowersFromTrigPrimsCreator(const edm::ParameterSet& iConfig) :
  m_caloTowersFromTrigPrimsAlgo( iConfig.getUntrackedParameter<double>("hadThreshold"),
				 iConfig.getUntrackedParameter<double>("emThreshold"),
				 iConfig.getUntrackedParameter<bool>("useHF"),
				 iConfig.getUntrackedParameter<bool>("verbose"))

{
  
  //-----------------------------------------------------
  // Get handles using input from the user
  //-----------------------------------------------------
  
  edm::InputTag d_hcalTrigPrimTag("simHcalTriggerPrimitiveDigis");
  edm::InputTag d_ecalTrigPrimTag("simEcalTriggerPrimitiveDigis");
  m_hcalTrigPrimTag = iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigPrimTag",d_hcalTrigPrimTag);
  m_ecalTrigPrimTag = iConfig.getUntrackedParameter<edm::InputTag>("ecalTrigPrimTag",d_ecalTrigPrimTag);
  
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

  using namespace std;

  //-----------------------------------------------------
  // Get all of your ESHandles
  //-----------------------------------------------------

  iSetup.get<CaloGeometryRecord>      ().get(m_geometry                        );
  iSetup.get<IdealGeometryRecord>     ().get(m_caloTowerConstituentsMap        );
  iSetup.get<IdealGeometryRecord>     ().get(m_ecalTrigTowerConstituentsMap    );
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel", m_ecalBarrelGeometry);
  iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap", m_ecalEndcapGeometry);
  iSetup.get<L1CaloEcalScaleRcd>      ().get(m_l1CaloEcalScale                 );
  iSetup.get<L1CaloHcalScaleRcd>      ().get(m_l1CaloHcalScale                 );
  
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

  m_caloTowersFromTrigPrimsAlgo.setGeometry    (m_geometry.product(),
					        m_caloTowerConstituentsMap.product(), 
					        m_ecalTrigTowerConstituentsMap.product(),
					        m_ecalBarrelGeometry.product(),
					        m_ecalEndcapGeometry.product());
  
  m_caloTowersFromTrigPrimsAlgo.setL1CaloScales(m_l1CaloEcalScale.product(),
					        m_l1CaloHcalScale.product() );   
  
  m_caloTowersFromTrigPrimsAlgo.resetEnergy();

  //------------------------------------------------------
  // Apply the algorithm
  //------------------------------------------------------

  m_caloTowersFromTrigPrimsAlgo.process(*ECALTrigPrimDigis);
  m_caloTowersFromTrigPrimsAlgo.process(*HCALTrigPrimDigis);

  //------------------------------------------------------
  // Fill the empty collection
  //------------------------------------------------------

  m_caloTowersFromTrigPrimsAlgo.finish(*caloTowerCollection);
  m_caloTowersFromTrigPrimsAlgo.checkEnergy();

  //------------------------------------------------------
  // Add the final CaloTowerCollection to the event
  //------------------------------------------------------

  iEvent.put(caloTowerCollection);

}

void CaloTowersFromTrigPrimsCreator::beginJob(const edm::EventSetup&){}

void CaloTowersFromTrigPrimsCreator::endJob() {}

DEFINE_FWK_MODULE(CaloTowersFromTrigPrimsCreator);
