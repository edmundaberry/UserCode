#include <memory>
#include <string>

#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsCreator.h"
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowerNonSortedCollection.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "FWCore/Utilities/interface/CPUTimer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1CaloEcalScaleRcd.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoLocalCalo/CaloTowersCreator/interface/EScales.h"

//------------------------------------------------------
// Constructor, 
//   get tags for data collections
//   and instruction bools here.
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::CaloTowersFromTrigPrimsCreator(const edm::ParameterSet& iConfig) :
  m_hcalTrigPrimTag ( iConfig.getParameter<edm::InputTag>("hcalTrigPrimTag")),
  m_ecalTrigPrimTag ( iConfig.getParameter<edm::InputTag>("ecalTrigPrimTag")),
  m_caloTowersFromTrigPrimsAlgo( new CaloTowersFromTrigPrimsAlgo (iConfig.getParameter<double>("hadThreshold"),
								  iConfig.getParameter<double>("emThreshold"),
								  iConfig.getParameter<bool>("useHF"),
								  iConfig.getParameter<bool>("verbose")))
  
  //{ produces<CaloTowerNonSortedCollection>("CaloTowerNonSortedCollection"); }
{ produces<CaloTowerCollection>(""); }

//------------------------------------------------------
// Destructor
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::~CaloTowersFromTrigPrimsCreator(){
  if ( m_caloTowersFromTrigPrimsAlgo != 0 ) delete m_caloTowersFromTrigPrimsAlgo;
}

//------------------------------------------------------
// Main production function
//------------------------------------------------------

void CaloTowersFromTrigPrimsCreator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //-----------------------------------------------------
  // Get all of your ESHandles
  //-----------------------------------------------------

  edm::ESHandle<CaloGeometry>                 geometry;
  edm::ESHandle<CaloSubdetectorGeometry>      ecalBarrelGeometry;
  edm::ESHandle<CaloSubdetectorGeometry>      ecalEndcapGeometry;
  edm::ESHandle<EcalTrigTowerConstituentsMap> ecalTrigTowerConstituentsMap;
  edm::ESHandle<CaloTowerConstituentsMap>     caloTowerConstituentsMap;
  edm::ESHandle<L1CaloEcalScale>              l1CaloEcalScale;
  edm::ESHandle<L1CaloHcalScale>              l1CaloHcalScale;
  
  iSetup.get<CaloGeometryRecord>      ().get(geometry                        );
  iSetup.get<IdealGeometryRecord>     ().get(caloTowerConstituentsMap        );
  iSetup.get<IdealGeometryRecord>     ().get(ecalTrigTowerConstituentsMap    );
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel", ecalBarrelGeometry);
  iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap", ecalEndcapGeometry);
  iSetup.get<L1CaloEcalScaleRcd>      ().get(l1CaloEcalScale                 );
  iSetup.get<L1CaloHcalScaleRcd>      ().get(l1CaloHcalScale                 );
  
  //-----------------------------------------------------
  // Get trigger primitives from the event
  //-----------------------------------------------------

  edm::Handle<HcalTrigPrimDigiCollection> HCALTrigPrimDigis;
  edm::Handle<EcalTrigPrimDigiCollection> ECALTrigPrimDigis;

  iEvent.getByLabel(m_hcalTrigPrimTag, HCALTrigPrimDigis );
  iEvent.getByLabel(m_ecalTrigPrimTag, ECALTrigPrimDigis );
  
  //------------------------------------------------------
  // Create an empty collection
  //------------------------------------------------------
  
  std::auto_ptr<CaloTowerCollection> caloTowerCollection(new CaloTowerCollection());

  //-----------------------------------------------------
  // Initialize the algorithm
  //-----------------------------------------------------

  m_caloTowersFromTrigPrimsAlgo -> setGeometry    (geometry.product(),
						   caloTowerConstituentsMap.product(), 
						   ecalTrigTowerConstituentsMap.product(),
						   ecalBarrelGeometry.product(),
						   ecalEndcapGeometry.product());
  
  m_caloTowersFromTrigPrimsAlgo -> setL1CaloScales(l1CaloEcalScale.product(),
						   l1CaloHcalScale.product() );   
  
  
  m_caloTowersFromTrigPrimsAlgo -> resetEnergy();
  
  //------------------------------------------------------
  // Apply the algorithm
  //------------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo -> process(*ECALTrigPrimDigis);
  m_caloTowersFromTrigPrimsAlgo -> process(*HCALTrigPrimDigis);

  //------------------------------------------------------
  // Fill the empty collection
  //------------------------------------------------------

  m_caloTowersFromTrigPrimsAlgo -> finish(*caloTowerCollection);
  m_caloTowersFromTrigPrimsAlgo -> checkEnergy();

  //------------------------------------------------------
  // Add the final CaloTowerCollection to the event
  //------------------------------------------------------

  //iEvent.put(caloTowerCollection,"CaloTowerNonSortedCollection");
  iEvent.put(caloTowerCollection);

}

void CaloTowersFromTrigPrimsCreator::beginJob(const edm::EventSetup&){}

void CaloTowersFromTrigPrimsCreator::endJob() {}

DEFINE_FWK_MODULE(CaloTowersFromTrigPrimsCreator);
