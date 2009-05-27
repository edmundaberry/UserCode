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
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// TPG transcoder info
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibCalorimetry/EcalTPGTools/interface/EcalTPGScale.h"

// Mapping object
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTrigTowerMap.h"

// Algorithm
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAlgo.h"

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
  InputTag m_hoTrigPrimTag;
  InputTag m_ecalTrigPrimTag;

  //------------------------------------------------------
  // ESHandles
  //------------------------------------------------------

  EcalTPGScale m_ecalTPGScale ;

  ESHandle<CaloTPGTranscoder>            m_caloTPGTranscoder;
  ESHandle<EcalTrigTowerConstituentsMap> m_ecalTrigTowerConstituentsMap;
  ESHandle<CaloTowerConstituentsMap>     m_caloTowerConstituentsMap;
  ESHandle<CaloGeometry>                 m_geometry;
  ESHandle<CaloSubdetectorGeometry>      m_ecalBarrelGeometry;
  ESHandle<CaloSubdetectorGeometry>      m_ecalEndcapGeometry;


  //------------------------------------------------------
  // Mapping from trig towers to calo towers
  //------------------------------------------------------

  CaloTrigTowerMap *m_caloTrigTowerMap;

  //------------------------------------------------------
  // Geometry
  //------------------------------------------------------

  HcalTrigTowerGeometry    m_trigTowerGeometry;

  //------------------------------------------------------
  // Algorithm
  //------------------------------------------------------

  CaloTowersFromTrigPrimsAlgo m_caloTowersFromTrigPrimsAlgo;

};

//------------------------------------------------------
// Constructor, get tags for data collections
//   and instruction bools here.
//------------------------------------------------------

CaloTowersFromTrigPrimsCreator::CaloTowersFromTrigPrimsCreator(const ParameterSet& iConfig) :
m_caloTowersFromTrigPrimsAlgo()
{

  InputTag d_hcalTrigPrimTag("simHcalTriggerPrimitiveDigis");
  m_hcalTrigPrimTag = iConfig.getUntrackedParameter<InputTag>("hcalTrigPrimTag",d_hcalTrigPrimTag);

  InputTag d_ecalTrigPrimTag("simEcalTriggerPrimitiveDigis");
  m_ecalTrigPrimTag = iConfig.getUntrackedParameter<InputTag>("ecalTrigPrimTag",d_ecalTrigPrimTag);
  
  InputTag d_hoTrigPrimTag("");
  m_hoTrigPrimTag = iConfig.getUntrackedParameter<InputTag>("hoTrigPrimTag",d_hoTrigPrimTag);

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

  iSetup.get<CaloGeometryRecord>      ().get(m_geometry                        );
  iSetup.get<IdealGeometryRecord>     ().get(m_caloTowerConstituentsMap        );
  iSetup.get<IdealGeometryRecord>     ().get(m_ecalTrigTowerConstituentsMap    );
  iSetup.get<CaloTPGRecord>           ().get(m_caloTPGTranscoder               );
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel", m_ecalBarrelGeometry);
  iSetup.get<EcalEndcapGeometryRecord>().get("EcalEndcap", m_ecalEndcapGeometry);

  m_ecalTPGScale.setEventSetup(iSetup);

  //-----------------------------------------------------
  // Set up the mapping objects and give them geometries
  //-----------------------------------------------------

  m_caloTrigTowerMap = new CaloTrigTowerMap();

  m_caloTrigTowerMap -> setGeometry(m_geometry.product(),
				    m_caloTowerConstituentsMap.product(), 
				    m_ecalTrigTowerConstituentsMap.product());
  
  m_caloTowersFromTrigPrimsAlgo.setGeometry(m_geometry.product(),
					    m_caloTowerConstituentsMap.product(), 
					    m_ecalTrigTowerConstituentsMap.product(),
					    m_ecalBarrelGeometry.product(),
					    m_ecalEndcapGeometry.product() 	       );
  
  m_caloTowersFromTrigPrimsAlgo.setCoder       (m_caloTPGTranscoder.product());
  m_caloTowersFromTrigPrimsAlgo.setEcalTPGScale(m_ecalTPGScale);
  m_caloTowersFromTrigPrimsAlgo.setMap         (m_caloTrigTowerMap);

  //-----------------------------------------------------
  // Get trigger primitives 
  //-----------------------------------------------------
  
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

  Handle<HOTrigPrimDigiCollection> HOTrigPrimDigis;
  bool hoTrigPrimDigiTagExists = iEvent.getByType(HOTrigPrimDigis);
  if (!hoTrigPrimDigiTagExists){
    LogWarning("CaloTowersFromTrigPrimsCreator") << "Could not extract HO trigger primitives with label: " << m_hoTrigPrimTag;
    return;
  }

  //------------------------------------------------------
  // Create an empty collection
  //------------------------------------------------------
  
  auto_ptr<CaloTowerCollection> caloTowerCollection(new CaloTowerCollection());

  //------------------------------------------------------
  // Apply the algorithm
  //------------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo.process(*HCALTrigPrimDigis);
  m_caloTowersFromTrigPrimsAlgo.process(*HOTrigPrimDigis);
  m_caloTowersFromTrigPrimsAlgo.process(*ECALTrigPrimDigis);

  //------------------------------------------------------
  // Finish the algorithm
  //------------------------------------------------------
  
  m_caloTowersFromTrigPrimsAlgo.finish(*caloTowerCollection);

  iEvent.put(caloTowerCollection,"TEST123");

}

void CaloTowersFromTrigPrimsCreator::beginJob(const EventSetup&){}

void CaloTowersFromTrigPrimsCreator::endJob() {}

DEFINE_FWK_MODULE(CaloTowersFromTrigPrimsCreator);
