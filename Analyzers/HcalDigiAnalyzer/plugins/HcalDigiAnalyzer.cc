// system include files
#include <memory>
#include <iostream>

// My includes
#include "Analyzers/HcalDigiAnalyzer/interface/HcalDigiAnalyzer.h"

// HCAL includes
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"

// Framework includes
#include "FWCore/Framework/interface/ESHandle.h"

HcalDigiAnalyzer::HcalDigiAnalyzer(const edm::ParameterSet& iConfig):
  m_outPath   (iConfig.getUntrackedParameter<std::string>  ("outPath","./") ),
  m_outSuffix (iConfig.getUntrackedParameter<std::string>  ("outSuffix",""))
{
  m_rootFile = m_outPath + std::string("HcalFlatNtuple") + m_outSuffix + std::string(".root");
  m_fillDigi.init(m_rootFile, &m_digiTree);
}

HcalDigiAnalyzer::~HcalDigiAnalyzer(){}

template <class DigiCollection, class RecoCollection > 
void HcalDigiAnalyzer::analyzeDigiCollection ( const HcalDbService & conditions,  const DigiCollection & digis , const RecoCollection & recos ){

  //-----------------------------------------------------
  // Get iterators
  //-----------------------------------------------------
   
  typename RecoCollection::const_iterator reco;
  typename RecoCollection::const_iterator reco_end = recos.end(); 
  typename DigiCollection::const_iterator digi     = digis.begin();
  typename DigiCollection::const_iterator digi_end = digis.end(); 

  //-----------------------------------------------------
  // Get setup information
  //-----------------------------------------------------
  
  const HcalQIEShape* shape = conditions.getHcalShape();

  //-----------------------------------------------------
  // 
  //-----------------------------------------------------

  int ndigis = (int) digis.size();
  int idigi = 0;

  m_digiTree.digi_timeslice_dv        -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_er        -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_raw       -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_adc       -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_nomFC     -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_fiber     -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_fiberChan -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_capid     -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_allFC     -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_pedFC     -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_gain      -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_rcgain    -> resize ( ndigis ) ;
  m_digiTree.digi_timeslice_energy    -> resize ( ndigis ) ;

  m_digiTree.digi_subdet              -> resize ( ndigis ) ;
  m_digiTree.digi_ieta                -> resize ( ndigis ) ;
  m_digiTree.digi_iphi                -> resize ( ndigis ) ;
  m_digiTree.digi_depth               -> resize ( ndigis ) ;
  m_digiTree.digi_presamples          -> resize ( ndigis ) ;
  m_digiTree.digi_nTS                 -> resize ( ndigis ) ;
  m_digiTree.digi_fiberIdleOffset     -> resize ( ndigis ) ;
  
  //-----------------------------------------------------
  // Loop through digis
  //-----------------------------------------------------
  
  for (; digi != digi_end ; ++digi ) {

    //-----------------------------------------------------
    // Get digi-specific cc objects
    //-----------------------------------------------------
    
    const HcalDetId       * hcalDetId    = & digi -> id();
    const HcalQIECoder    * channelCoder =   conditions.getHcalCoder( *hcalDetId );
    const HcalCalibrations* calibrations = & conditions.getHcalCalibrations ( *hcalDetId ) ;
    
    HcalCoderDb coder (*channelCoder, *shape); 
    CaloSamples tool;
    coder.adc2fC ( * digi, tool );

    //-----------------------------------------------------
    // Get digi-specific values
    //-----------------------------------------------------
        
    int digi_subdet          = hcalDetId -> subdet();
    int digi_ieta            = hcalDetId -> ieta();
    int digi_iphi            = hcalDetId -> iphi();
    int digi_depth           = hcalDetId -> depth();
    int digi_presamples      = digi      -> presamples() ;
    int digi_nTS             = digi      -> size();
    int digi_fiberIdleOffset = digi      -> fiberIdleOffset();

    (*m_digiTree.digi_subdet         )[idigi] = digi_subdet          ;
    (*m_digiTree.digi_ieta           )[idigi] = digi_ieta            ;
    (*m_digiTree.digi_iphi           )[idigi] = digi_iphi            ;
    (*m_digiTree.digi_depth          )[idigi] = digi_depth           ;
    (*m_digiTree.digi_presamples     )[idigi] = digi_presamples      ;
    (*m_digiTree.digi_nTS            )[idigi] = digi_nTS             ;
    (*m_digiTree.digi_fiberIdleOffset)[idigi] = digi_fiberIdleOffset ;
    
    //-----------------------------------------------------
    // Loop through digi time slices
    //-----------------------------------------------------

    (*m_digiTree.digi_timeslice_dv        )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_er        )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_raw       )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_adc       )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_nomFC     )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_fiber     )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_fiberChan )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_capid     )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_allFC     )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_pedFC     )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_gain      )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_rcgain    )[idigi].resize ( digi_nTS );
    (*m_digiTree.digi_timeslice_energy    )[idigi].resize ( digi_nTS );

    for ( int iTS = 0; iTS < digi_nTS ; ++iTS ) {
      
      //-----------------------------------------------------
      // Get slice-specific cc objects
      //-----------------------------------------------------

      const HcalQIESample * qieSample = & digi -> sample (iTS);

      //-----------------------------------------------------
      // Get slice-specific values
      //-----------------------------------------------------

      int    digi_timeslice_dv        = (int) qieSample -> dv();
      int    digi_timeslice_er        = (int) qieSample -> er();
      int    digi_timeslice_raw       = qieSample -> raw();
      int    digi_timeslice_adc       = qieSample -> adc();
      int    digi_timeslice_nomFC     = qieSample -> nominal_fC();
      int    digi_timeslice_fiber     = qieSample -> fiber();
      int    digi_timeslice_fiberChan = qieSample -> fiberChan();
      int    digi_timeslice_capid     = qieSample -> capid();
      double digi_timeslice_allFC     = tool[iTS];
      double digi_timeslice_pedFC     = calibrations -> pedestal     ( digi_timeslice_capid );
      double digi_timeslice_gain      = calibrations -> rawgain      ( digi_timeslice_capid );
      double digi_timeslice_rcgain    = calibrations -> respcorrgain ( digi_timeslice_capid ) ;
      double digi_timeslice_energy    = ( digi_timeslice_allFC - digi_timeslice_pedFC ) * digi_timeslice_rcgain;

      (*m_digiTree.digi_timeslice_dv        )[idigi][iTS] = digi_timeslice_dv       ;
      (*m_digiTree.digi_timeslice_er        )[idigi][iTS] = digi_timeslice_er       ;
      (*m_digiTree.digi_timeslice_raw       )[idigi][iTS] = digi_timeslice_raw      ;
      (*m_digiTree.digi_timeslice_adc       )[idigi][iTS] = digi_timeslice_adc      ;
      (*m_digiTree.digi_timeslice_nomFC     )[idigi][iTS] = digi_timeslice_nomFC    ;
      (*m_digiTree.digi_timeslice_fiber     )[idigi][iTS] = digi_timeslice_fiber    ;
      (*m_digiTree.digi_timeslice_fiberChan )[idigi][iTS] = digi_timeslice_fiberChan;
      (*m_digiTree.digi_timeslice_capid     )[idigi][iTS] = digi_timeslice_capid    ;
      (*m_digiTree.digi_timeslice_allFC     )[idigi][iTS] = digi_timeslice_allFC    ;
      (*m_digiTree.digi_timeslice_pedFC     )[idigi][iTS] = digi_timeslice_pedFC    ;
      (*m_digiTree.digi_timeslice_gain      )[idigi][iTS] = digi_timeslice_gain     ;
      (*m_digiTree.digi_timeslice_rcgain    )[idigi][iTS] = digi_timeslice_rcgain   ;
      (*m_digiTree.digi_timeslice_energy    )[idigi][iTS] = digi_timeslice_energy   ;
      
    }

    //-----------------------------------------------------
    // For each digi, try to find a rechit
    //-----------------------------------------------------

    reco = recos.find ( * hcalDetId ) ;

    double reco_energy = -999.;
    double reco_time   = -999.;

    if ( reco != reco_end ) {
      reco_energy = reco -> energy();
      reco_time   = reco -> time();
    }
    
    (*m_digiTree.reco_energy).push_back ( reco_energy ) ;
    (*m_digiTree.reco_energy).push_back ( reco_time   ) ;

    ++idigi;
    
  } // end of loop over digis
}

void HcalDigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


  //-----------------------------------------------------
  // Fill the run and event number information
  //-----------------------------------------------------
  
  int run   = iEvent.id().run();
  int event = iEvent.id().event();  

  //-----------------------------------------------------
  // Initialize the Tree values
  //-----------------------------------------------------

  m_digiTree.clear();
  
  (*m_digiTree.run  ).push_back ( run );
  (*m_digiTree.event).push_back ( event ) ;
  
  //-----------------------------------------------------
  // edm::ESHandles
  //-----------------------------------------------------

  edm::ESHandle<HcalDbService> conditions;

  //-----------------------------------------------------
  // edm::Handles
  //-----------------------------------------------------

  edm::Handle<HBHEDigiCollection  >  hbheDigis;
  edm::Handle<HBHERecHitCollection>  hbheRecos;

  edm::Handle<HFDigiCollection  >  hfDigis;
  edm::Handle<HFRecHitCollection>  hfRecos;

  edm::Handle<HODigiCollection  >  hoDigis;
  edm::Handle<HORecHitCollection>  hoRecos;

  //-----------------------------------------------------
  // Get EventSetup objects
  //-----------------------------------------------------

  iSetup.get<HcalDbRecord>().get(conditions);

  //-----------------------------------------------------
  // Get EDM event objects
  //-----------------------------------------------------

  bool gotHBHEDigis = iEvent.getByLabel("hcalDigis", hbheDigis);
  if (!gotHBHEDigis ) {
    std::cout << "Could not find HBHE Digis" << std::endl;
    exit(0);
  }

  bool gotHBHERecos = iEvent.getByLabel("hbheprereco", hbheRecos);
  if (!gotHBHERecos ) {
    std::cout << "Could not find HBHE Recos" << std::endl;
    exit(0);
  }

  bool gotHFDigis = iEvent.getByLabel("hcalDigis", hfDigis);
  if (!gotHFDigis ) {
    std::cout << "Could not find HF Digis" << std::endl;
    exit(0);
  }

  bool gotHFRecos = iEvent.getByLabel("hfreco", hfRecos);
  if (!gotHFRecos ) {
    std::cout << "Could not find HF Recos" << std::endl;
    exit(0);
  }

  bool gotHODigis = iEvent.getByLabel("hcalDigis", hoDigis);
  if (!gotHODigis ) {
    std::cout << "Could not find HO Digis" << std::endl;
    exit(0);
  }

  bool gotHORecos = iEvent.getByLabel("horeco", hoRecos);
  if (!gotHORecos ) {
    std::cout << "Could not find HO Recos" << std::endl;
    exit(0);
  }

  //-----------------------------------------------------
  // Analyze digis
  //-----------------------------------------------------

  analyzeDigiCollection ( *conditions  ,  *hbheDigis , *hbheRecos ) ;
  analyzeDigiCollection ( *conditions  ,  *hfDigis   , *hfRecos   ) ;
  analyzeDigiCollection ( *conditions  ,  *hoDigis   , *hoRecos   ) ;

  //-----------------------------------------------------
  // Fill tree
  //-----------------------------------------------------

  m_fillDigi.fill();

}

void HcalDigiAnalyzer::beginJob(){}

void HcalDigiAnalyzer::endJob() {
  m_fillDigi.finalize();
}

DEFINE_FWK_MODULE(HcalDigiAnalyzer);
