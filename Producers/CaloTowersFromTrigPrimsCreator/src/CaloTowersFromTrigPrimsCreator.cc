// System files
#include <memory>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <iomanip>

// Framework files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Data collections
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// Geometry
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// Mapping object
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTrigTowerMapBuilder.h"

using namespace edm;
using namespace std;

class CaloTowersFromTrigPrimsCreator : public EDProducer {
public:
  explicit CaloTowersFromTrigPrimsCreator(const ParameterSet&);
  ~CaloTowersFromTrigPrimsCreator();
  
private:
  virtual void beginJob(const EventSetup&) ;
  virtual void produce(Event&, const EventSetup&);
  virtual void endJob() ;
    
  //------------------------------------------------------
  // Instruction bools
  //------------------------------------------------------
    
  bool m_verbose;

  //------------------------------------------------------
  // InputTags 
  //------------------------------------------------------
  
  InputTag m_hcalTrigPrimTag;
  InputTag m_ecalTrigPrimTag;

  //------------------------------------------------------
  // ESHandles
  //------------------------------------------------------

  ESHandle<CaloTowerConstituentsMap> m_caloTowerConstituentsMap;
  ESHandle<CaloGeometry> m_geometry;

  //------------------------------------------------------
  // My helper objects
  //------------------------------------------------------

  CaloIdCompiler          *m_caloIdCompiler;
  CaloTrigTowerMapBuilder *m_caloTrigTowerMapBuilder;

  //------------------------------------------------------
  // Geometry
  //------------------------------------------------------
  
  HcalTrigTowerGeometry    m_trigTowerGeometry;

  //------------------------------------------------------
  // typedefs from CaloTrigTowerMapBuilder header file
  // Do I really need to be re-doing these?
  //------------------------------------------------------
  
  typedef CaloTrigTowerMapBuilder::CaloToTrigTowerMap CaloToTrigTowerMap;
  typedef CaloTrigTowerMapBuilder::TrigToCaloTowerMap TrigToCaloTowerMap;

};

//------------------------------------------------------
// Constructor, get tags for data collections
//   and instruction bools here.
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::CaloTowersFromTrigPrimsCreator(const ParameterSet& iConfig){

  InputTag d_hcalTrigPrimTag("hcalDigis");
  m_hcalTrigPrimTag = iConfig.getUntrackedParameter<InputTag>("hcalTrigPrimTag",d_hcalTrigPrimTag);

  InputTag d_ecalTrigPrimTag("ecalDigis","EcalTriggerPrimitives");
  m_ecalTrigPrimTag = iConfig.getUntrackedParameter<InputTag>("ecalTrigPrimTag",d_ecalTrigPrimTag);
  
  bool d_verbose = true;
  m_verbose = iConfig.getUntrackedParameter("verbose",d_verbose);

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

  iSetup.get<CaloGeometryRecord>().get(m_geometry);
  iSetup.get<IdealGeometryRecord>().get(m_caloTowerConstituentsMap);
  
  //-----------------------------------------------------
  // Set up the mapping objects and give them geometries
  //-----------------------------------------------------

  m_caloIdCompiler          = new CaloIdCompiler();
  m_caloTrigTowerMapBuilder = new CaloTrigTowerMapBuilder();
  m_caloIdCompiler          -> setGeometry(m_geometry.product(),m_caloTowerConstituentsMap.product());
  m_caloTrigTowerMapBuilder -> setGeometry(m_geometry.product(),m_caloTowerConstituentsMap.product());
  
  m_caloTrigTowerMapBuilder -> buildMap();
				 
  //-----------------------------------------------------
  // Retrieve maps
  //-----------------------------------------------------

  CaloToTrigTowerMap caloToTrigTowerMap = m_caloTrigTowerMapBuilder -> getCaloToTrigTowerMap();
  TrigToCaloTowerMap trigToCaloTowerMap = m_caloTrigTowerMapBuilder -> getTrigToCaloTowerMap();
  
  //-----------------------------------------------------
  // Get vectors of all valid TrigTowerDetId's
  //----------------------------------------------------- 
  
  vector<HcalTrigTowerDetId> trigTowerDetIds = m_caloIdCompiler -> getAllHcalTrigTowerDetIds();
  vector<CaloTowerDetId>     caloTowerDetIds = m_caloIdCompiler -> getAllCaloTowerDetIds    ();

  cout << "trigTowerDetIds have size: " << trigTowerDetIds.size() << endl;
  cout << "caloTowerDetIds have size: " << caloTowerDetIds.size() << endl;

  vector<HcalTrigTowerDetId>::iterator trigTowerDetId_iter = trigTowerDetIds.begin();
  vector<CaloTowerDetId>::iterator     caloTowerDetId_iter = caloTowerDetIds.begin();
  
  //-----------------------------------------------------
  // Get real handles of trigger primitives 
  // (This will be more useful later)
  //-----------------------------------------------------
  
  /*
  Handle<HcalTrigPrimDigiCollection> HCALTrigPrimDigis;
  bool hcalTrigPrimDigiTagExists = iEvent.getByLabel(m_hcalTrigPrimTag,HCALTrigPrimDigis);
  if (!hcalTrigPrimDigiTagExists){
    LogWarning("CaloTowersFromTrigPrimsCreator") << "Could not extract HCAL trigger primitives with label: " << m_hcalTrigPrimTag;
    return;
  }

  Handle<EcalTrigPrimDigiCollection> ECALTrigPrimDigis;
  bool ecalTrigPrimDigiTagExists = iEvent.getByLabel(m_ecalTrigPrimTag,ECALTrigPrimDigis);
  if (!ecalTrigPrimDigiTagExists){
    LogWarning("CaloTowersFromTrigPrimsCreator") << "Could not extract ECAL trigger primitives with label: " << m_ecalTrigPrimTag;
    return;
  }
  */

  //------------------------------------------------------
  // Loop over the HCAL trig prims first
  //------------------------------------------------------
  
  /*
  HcalTrigPrimDigiCollection::const_iterator hcalTrigPrimDigi = HCALTrigPrimDigis -> begin();
  
  for (; hcalTrigPrimDigi != HCALTrigPrimDigis -> end(); hcalTrigPrimDigi++){

    vector<HcalDetId> HcalDetIds = m_trigTowerGeometry.detIds((*hcalTrigPrimDigi).id());
  */
  

  //------------------------------------------------------
  // 1.  Loop over all CaloTowerDetId's
  // 2.  See how many different trigger towers each was mapped to
  // 3.  You hope that each CaloTower was only mapped to one trigger tower
  //------------------------------------------------------

  vector<CaloTowerDetId> multiplyMappedCaloTowerDetIds;
  vector<CaloTowerDetId>::iterator thisMultiplyMappedCaloTowerDetId;

  for(; caloTowerDetId_iter != caloTowerDetIds.end(); ++caloTowerDetId_iter){

    int entries = caloToTrigTowerMap.count(*caloTowerDetId_iter);

    if (entries > 1) {   
      multiplyMappedCaloTowerDetIds.push_back(*caloTowerDetId_iter);      
    }    

  }
  
  cout << multiplyMappedCaloTowerDetIds.size() << " CaloTower(s) with constituents that were assigned to different trigger towers" << endl;

  for (thisMultiplyMappedCaloTowerDetId  = multiplyMappedCaloTowerDetIds.begin();
       thisMultiplyMappedCaloTowerDetId != multiplyMappedCaloTowerDetIds.end();
       thisMultiplyMappedCaloTowerDetId++){    
    cout << "  " << *thisMultiplyMappedCaloTowerDetId;
    
    vector<DetId> constituents = m_caloTowerConstituentsMap.product() -> constituentsOf(*thisMultiplyMappedCaloTowerDetId);
    vector<DetId>::iterator thisConstituent = constituents.begin();

    cout << " (" << constituents.size() << ") ";

    for (; thisConstituent != constituents.end(); thisConstituent++){
      if ( (*thisConstituent).det() == DetId::Hcal ){	
	HcalDetId hcalDetId = HcalDetId(*thisConstituent); 

	int nTriggerTowers = m_trigTowerGeometry.towerIds(hcalDetId).size();

	if (nTriggerTowers != 2) cout << hcalDetId << " ";
      }
    }
    cout << endl;
    
    
  }  
}

void CaloTowersFromTrigPrimsCreator::beginJob(const EventSetup&){}

void CaloTowersFromTrigPrimsCreator::endJob() {}

DEFINE_FWK_MODULE(CaloTowersFromTrigPrimsCreator);
