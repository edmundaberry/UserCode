#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiComparer.h"
#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiCollection.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

HcalSLHCTrigPrimDigiComparer::HcalSLHCTrigPrimDigiComparer(const edm::ParameterSet& iConfig):
  m_upgradeDigisTag(iConfig.getParameter<edm::InputTag>("upgradeTrigPrimTag")),
  m_defaultDigisTag(iConfig.getParameter<edm::InputTag>("defaultTrigPrimTag"))
{}

HcalSLHCTrigPrimDigiComparer::~HcalSLHCTrigPrimDigiComparer(){}


void HcalSLHCTrigPrimDigiComparer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //------------------------------------------------------
  // First get handles
  //------------------------------------------------------

  edm::Handle<HcalSLHCTrigPrimDigiCollection> UpgradeDigis;
  edm::Handle<HcalTrigPrimDigiCollection>     DefaultDigis;

  bool gotUpgrade = iEvent.getByLabel ( m_upgradeDigisTag , UpgradeDigis );
  bool gotDefault = iEvent.getByLabel ( m_defaultDigisTag , DefaultDigis );
  
  if (!gotUpgrade) { 
    edm::LogWarning("HcalSLHCTrigPrimDigiComparer") << "Cannot get " << m_upgradeDigisTag; 
    return; 
  }
  
  if (!gotDefault) { 
    edm::LogWarning("HcalSLHCTrigPrimDigiComparer") << "Cannot get " << m_defaultDigisTag; 
    return; 
  }

  //------------------------------------------------------
  // Set iterators
  //------------------------------------------------------
  
  HcalTrigPrimDigiCollection::const_iterator defaultDigi     = DefaultDigis -> begin();
  HcalTrigPrimDigiCollection::const_iterator defaultDigi_end = DefaultDigis -> end();

  HcalSLHCTrigPrimDigiCollection::const_iterator upgradeDigi;
  HcalSLHCTrigPrimDigiCollection::const_iterator upgradeDigi_end = UpgradeDigis -> end();

  //------------------------------------------------------
  // Now perform the first of three checks
  //
  // Check #1: Are there the same number of digis in each
  // collection?  There should be.
  //------------------------------------------------------

  if (DefaultDigis -> size() != UpgradeDigis -> size()) {
    std::cout << "***ERROR: Size mismatch!" << std::endl;
    std::cout << "   Default digis have size: " << DefaultDigis -> size() << std::endl;
    std::cout << "   Upgrade digis have size: " << UpgradeDigis -> size() << std::endl;
  }

  //------------------------------------------------------
  // Check #2: Is the compressedEt bit word the same 
  // for upgrade digis as for default digis? It should be.
  //
  // Also, is there an upgrade digi for each default digi?
  //------------------------------------------------------

  int upgrade_SOI_compressedIsoEt;
  int upgrade_SOI_compressedEt   ;
  int default_SOI_compressedEt   ;

  for (; defaultDigi != defaultDigi_end; ++defaultDigi){
    HcalTrigTowerDetId default_id = defaultDigi -> id ();

    upgradeDigi = UpgradeDigis -> find ( default_id );

    if ( upgradeDigi == upgradeDigi_end ){
      std::cout << "***ERROR: No Upgrade digi to match the default digi at: " << default_id << std::endl;
      continue;
    }

    upgrade_SOI_compressedIsoEt = upgradeDigi -> SOI_compressedIsoEt();
    upgrade_SOI_compressedEt    = upgradeDigi -> SOI_compressedEt();
    default_SOI_compressedEt    = defaultDigi -> SOI_compressedEt();

    if ( upgrade_SOI_compressedEt != default_SOI_compressedEt ){
      std::cout << "***ERROR: SOI compressedEt disagreement for digi at: " << default_id << std::endl;
      std::cout << "   Default SOI compressedEt = " << default_SOI_compressedEt << std::endl;
      std::cout << "   Upgrade SOI compressedEt = " << upgrade_SOI_compressedEt << std::endl;
    }
  }

  //------------------------------------------------------
  // Check #3: When you set the minimum isolation depth to
  // zero and the maximum depth to 20, do you see the
  // same values for isolation Et as for total Et?
  // You should.
  // 
  // NB: Ignore the HF, since we're not doing isolation there.
  //
  // Also, is there a default digi for each upgrade digi?
  //------------------------------------------------------

  upgradeDigi = UpgradeDigis -> begin();

  for (; upgradeDigi != upgradeDigi_end; ++upgradeDigi){

    HcalTrigTowerDetId upgrade_id = upgradeDigi -> id ();

    bool isHF = ( upgrade_id.ietaAbs () >= 29 );

    defaultDigi = DefaultDigis -> find (upgrade_id);

    if ( defaultDigi == defaultDigi_end ){
      std::cout << "***ERROR: No Default digi to match the upgrade digi at: " << upgrade_id << std::endl;
      continue;
    }

    upgrade_SOI_compressedIsoEt = upgradeDigi -> SOI_compressedIsoEt();
    upgrade_SOI_compressedEt    = upgradeDigi -> SOI_compressedEt();
    default_SOI_compressedEt    = defaultDigi -> SOI_compressedEt();

    if ( upgrade_SOI_compressedEt != default_SOI_compressedEt ){
      std::cout << "***ERROR: SOI compressedEt disagreement for digi at: " << upgrade_id << std::endl;
      std::cout << "   Default SOI compressedEt = " << default_SOI_compressedEt << std::endl;
      std::cout << "   Upgrade SOI compressedEt = " << upgrade_SOI_compressedEt << std::endl;
    }

    if ( upgrade_SOI_compressedEt != upgrade_SOI_compressedIsoEt && (!isHF) ){
      std::cout << "***ERROR: Upgrade SOI et/isoEt disagreement for digi at: " << upgrade_id << std::endl;
      std::cout << "   Upgrade SOI compressedEt    = " << default_SOI_compressedEt    << std::endl;
      std::cout << "   Upgrade SOI compressedIsoEt = " << upgrade_SOI_compressedIsoEt << std::endl;
    }

  }
       

}


void HcalSLHCTrigPrimDigiComparer::beginJob(){}

void HcalSLHCTrigPrimDigiComparer::endJob() {}

DEFINE_FWK_MODULE(HcalSLHCTrigPrimDigiComparer);
