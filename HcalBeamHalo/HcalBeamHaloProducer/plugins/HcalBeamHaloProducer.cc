#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "HcalBeamHalo/HcalBeamHaloProducer/interface/HcalBeamHaloCollection.h"
#include "HcalBeamHalo/HcalBeamHaloProducer/interface/HcalBeamHaloAlgo.h"

class HcalBeamHaloProducer : public edm::EDProducer {
public:
  explicit HcalBeamHaloProducer(const edm::ParameterSet&);
  ~HcalBeamHaloProducer();
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag m_hbheRecHitTag;
  edm::InputTag m_standAloneTrackTag;
  
  std::vector<int> m_rechit_blacklist;

  double m_rechit_minEnergy  ;
  double m_window_minEnergy  ;   
  int    m_window_minCounts  ;
  int    m_window_iphiWidth  ;
  int    m_max_n_windows     ; 
  bool   m_verbose           ;
  
  HcalBeamHaloAlgo * m_algo;

};

HcalBeamHaloProducer::HcalBeamHaloProducer(const edm::ParameterSet& iConfig):
  m_hbheRecHitTag      ( iConfig.getParameter<edm::InputTag>("HBHERecHits") ),
  m_standAloneTrackTag ( iConfig.getParameter<edm::InputTag>("StandAloneTracks") ),
  m_rechit_blacklist   ( iConfig.getParameter<std::vector<int> >("BlackListCells")),
  m_rechit_minEnergy   ( iConfig.getParameter<double>("MinRecHitEnergy") ),            
  m_window_minEnergy   ( iConfig.getParameter<double>("MinWindowEnergy") ),
  m_window_minCounts   ( iConfig.getParameter<int>("MinWindowCounts") ),
  m_window_iphiWidth   ( iConfig.getParameter<int>("Width") ),
  m_max_n_windows      ( iConfig.getParameter<int>("MaxNWindows")),
  m_verbose            ( iConfig.getParameter<bool>("Verbose")),
  m_algo ( new HcalBeamHaloAlgo ( m_rechit_blacklist, m_rechit_minEnergy, m_window_minEnergy, m_window_minCounts, m_window_iphiWidth, m_max_n_windows, m_verbose ))
{
  produces<HcalBeamHaloCollection>("");
}

HcalBeamHaloProducer::~HcalBeamHaloProducer(){
  if (m_algo) delete m_algo;
}

void HcalBeamHaloProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::auto_ptr<HcalBeamHaloCollection> hcalBeamHaloCollection ( new HcalBeamHaloCollection() );

  edm::Handle<HBHERecHitCollection> hbheRecHitCollection;
  edm::Handle<reco::TrackCollection> standAloneTrackCollection;

  edm::ESHandle<CaloGeometry> calo_geometry_handle;
  iSetup.get<CaloGeometryRecord>().get(calo_geometry_handle);
  
  iEvent.getByLabel(m_standAloneTrackTag, standAloneTrackCollection );
  iEvent.getByLabel(m_hbheRecHitTag, hbheRecHitCollection);
  
  m_algo -> setGeom ( calo_geometry_handle.product() );
  m_algo -> process ( *hbheRecHitCollection, * standAloneTrackCollection );
  m_algo -> getHalo ( hcalBeamHaloCollection);

  iEvent.put(hcalBeamHaloCollection);
  
}

void HcalBeamHaloProducer::beginJob(){}
void  HcalBeamHaloProducer::endJob() {}

DEFINE_FWK_MODULE(HcalBeamHaloProducer);
  
