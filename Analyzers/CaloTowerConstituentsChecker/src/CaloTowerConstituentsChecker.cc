// system include files
#include <memory>
#include <algorithm>
#include <vector>
#include <string>
#include <TH2F.h>
#include <TFile.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// My include files
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

using namespace std;
using namespace edm;

class CaloTowerConstituentsChecker : public edm::EDAnalyzer {
public:
  explicit CaloTowerConstituentsChecker(const edm::ParameterSet&);
  ~CaloTowerConstituentsChecker();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void checkForCaloTowerMapDisagreement(std::vector<DetId> cells);
  void checkForEcalDoubleCounting (std::vector<DetId> cells);
  void checkForHcalDoubleCounting (std::vector<DetId> cells);

  void fillEcalCountingMap(std::vector<CaloTowerDetId> caloTowerIds);
  void fillHcalCountingMap(std::vector<CaloTowerDetId> caloTowerIds);
  
  void getCaloTowerDetIds_fromSubdet(const std::vector<DetId> & cellIds, 
				     std::vector<CaloTowerDetId> & towerIds);
  
  std::vector<CaloTowerDetId> getCaloTowerDetIds(std::vector<DetId> ebCells, std::vector<DetId> eeCells,
						 std::vector<DetId> hbCells, std::vector<DetId> heCells,
						 std::vector<DetId> hoCells, std::vector<DetId> hfCells );
  
  void setup();
  
  typedef multimap<uint32_t, int, less<uint32_t> > CaloTowerDetIdCountMap;
  typedef multimap<uint32_t, int, less<uint32_t> >::value_type CaloTowerDetIdCount;
  CaloTowerDetIdCountMap caloTowerDetIdCountMap;
  
  typedef multimap<DetId, CaloTowerDetId, less<DetId> > EcalDetIdCountMap;
  typedef multimap<DetId, CaloTowerDetId, less<DetId> >::value_type EcalDetIdCount;
  EcalDetIdCountMap ecalDetIdCountMap;

  typedef multimap<DetId, CaloTowerDetId, less<DetId> > HcalDetIdCountMap;
  typedef multimap<DetId, CaloTowerDetId, less<DetId> >::value_type HcalDetIdCount;
  HcalDetIdCountMap hcalDetIdCountMap;

  const CaloTowerConstituentsMap * theCaloTowerConstituentsMap;

  int valid_unassigned_detids;
  int valid_multicount_detids;
  int valid_doublecount_detids;
  int valid_doublecount_ok_detids;
  int valid_singlecount_detids;
  int valid_singlecount_ok_detids;
  int valid_total_detids;

  bool verbose;

  bool check_eb_doublecount;
  bool check_ee_doublecount;
  bool check_hb_doublecount;
  bool check_he_doublecount;
  bool check_ho_doublecount;
  bool check_hf_doublecount;
  
  string filename;

  TFile *plotFile;
  TH2F  *h_doublecount_occ[5];

};

CaloTowerConstituentsChecker::CaloTowerConstituentsChecker(const edm::ParameterSet& iConfig){

  verbose              = iConfig.getUntrackedParameter<bool>  ("verbose"             ,true);

  check_eb_doublecount = iConfig.getUntrackedParameter<bool>  ("check_eb_doublecount",true);
  check_ee_doublecount = iConfig.getUntrackedParameter<bool>  ("check_ee_doublecount",true);

  check_hb_doublecount = iConfig.getUntrackedParameter<bool>  ("check_hb_doublecount",true);
  check_he_doublecount = iConfig.getUntrackedParameter<bool>  ("check_he_doublecount",true);
  check_ho_doublecount = iConfig.getUntrackedParameter<bool>  ("check_ho_doublecount",true);
  check_hf_doublecount = iConfig.getUntrackedParameter<bool>  ("check_hf_doublecount",true);

  filename             = iConfig.getUntrackedParameter<string>("filename", "double_count.root");

  plotFile = new TFile(filename.c_str(), "RECREATE");

  char depthName[500];

  for (int depth = 1; depth <= 4; depth++){
    sprintf(depthName,"h_doublecount_occ_depth%d",depth);
    h_doublecount_occ[depth] = new TH2F(depthName,"",83,-41.5,41.5,72,0.5,72.5);
  }


}

CaloTowerConstituentsChecker::~CaloTowerConstituentsChecker() {}

void CaloTowerConstituentsChecker::setup(){

  valid_unassigned_detids     = 0;
  valid_multicount_detids     = 0;
  valid_doublecount_detids    = 0;
  valid_doublecount_ok_detids = 0;
  valid_singlecount_detids    = 0;
  valid_singlecount_ok_detids = 0;
  valid_total_detids          = 0;

}

void CaloTowerConstituentsChecker::getCaloTowerDetIds_fromSubdet(const std::vector<DetId> & cellIds, 
								 std::vector<CaloTowerDetId> & towerIds){
  
  int entries;
  
  std::vector<DetId>::const_iterator cellId_iter = cellIds.begin();
  std::vector<DetId>::const_iterator cellId_end  = cellIds.end();

  for (; cellId_iter != cellId_end; ++cellId_iter){

    CaloTowerDetId temp = theCaloTowerConstituentsMap -> towerOf(*cellId_iter);       
    
    entries = caloTowerDetIdCountMap.count(temp.rawId());
    
    if (entries == 0) towerIds.push_back(temp);
    
    caloTowerDetIdCountMap.insert(CaloTowerDetIdCount(temp.rawId(),entries+1));

  }
}

std::vector<CaloTowerDetId> CaloTowerConstituentsChecker::getCaloTowerDetIds(std::vector<DetId> ebCells,
									     std::vector<DetId> eeCells,
									     std::vector<DetId> hbCells,
									     std::vector<DetId> heCells,
									     std::vector<DetId> hoCells,
									     std::vector<DetId> hfCells ){
  
  std::vector<CaloTowerDetId> retval;
  
  if (check_eb_doublecount) getCaloTowerDetIds_fromSubdet(ebCells, retval);
  if (check_ee_doublecount) getCaloTowerDetIds_fromSubdet(eeCells, retval);
  if (check_hb_doublecount) getCaloTowerDetIds_fromSubdet(hbCells, retval);
  if (check_he_doublecount) getCaloTowerDetIds_fromSubdet(heCells, retval);
  if (check_ho_doublecount) getCaloTowerDetIds_fromSubdet(hoCells, retval);
  if (check_hf_doublecount) getCaloTowerDetIds_fromSubdet(hfCells, retval);
  
  std::sort   (retval.begin(), retval.end());  
  std::unique (retval.begin(), retval.end());

  return retval;

}

void CaloTowerConstituentsChecker::checkForHcalDoubleCounting (std::vector<DetId> detIds){
  
  int entries;
  std::vector<DetId>::iterator detId_iter = detIds.begin();

  for (; detId_iter != detIds.end(); ++detId_iter){

    valid_total_detids++;

    DetId detId = DetId(*detId_iter);
    
    CaloTowerDetId caloTowerDetId = theCaloTowerConstituentsMap -> towerOf(detId);

    entries = hcalDetIdCountMap.count(detId);
    
    if      (entries == 0)  valid_unassigned_detids++;
    else if (entries == 1)  valid_singlecount_detids++;      
    else if (entries  > 1){
      
      HcalDetId hcalDetId = HcalDetId(detId);
      if (verbose){
	cout << "This HcalDetId showed up " << entries << " times" <<  endl;
	cout << "  " << hcalDetId << endl;	    
      }
      
      cout << "  It was mapped to CaloTowerDetId's at:" << endl;
      HcalDetIdCountMap::const_iterator hcalDetIdCountMap_iter =  hcalDetIdCountMap.begin();
      for (; hcalDetIdCountMap_iter != hcalDetIdCountMap.end(); ++hcalDetIdCountMap_iter){
	DetId tempDetId = hcalDetIdCountMap_iter -> first;
	CaloTowerDetId tempCaloTowerDetId = hcalDetIdCountMap_iter -> second;
	if (tempDetId == detId )
	  cout << "    ieta = " << caloTowerDetId.ieta() << ", iphi = " << caloTowerDetId.iphi() << endl;
      }  
      
      valid_multicount_detids++;
      
    }
    else cout << "Entries = " << entries << endl;
  }  
  
  
}

void CaloTowerConstituentsChecker::checkForEcalDoubleCounting (std::vector<DetId> detIds){

  int entries;
  std::vector<DetId>::iterator detId_iter = detIds.begin();

  for (; detId_iter != detIds.end(); ++detId_iter){

    valid_total_detids++;

    DetId detId = DetId(*detId_iter);
    
    CaloTowerDetId caloTowerDetId = theCaloTowerConstituentsMap -> towerOf(detId);

    entries = ecalDetIdCountMap.count(detId);
    
    if      (entries == 0)  valid_unassigned_detids++;
    else if (entries == 1)  valid_singlecount_detids++;      
    else if (entries  > 1){
      
      if (detId.subdetId() == EcalBarrel){
	EBDetId ebDetId = EBDetId(detId);
	if (verbose){
	  cout << "This EB Det ID showed up " << entries << " times" <<  endl;
	  cout << "  " << ebDetId << endl;	    
	}
      }

      if (detId.subdetId() == EcalEndcap){
	EEDetId eeDetId = EEDetId(detId);
	if (verbose){
	  cout << "This EE Det ID showed up " << entries << " times" <<  endl;
	  cout << "  " << eeDetId << endl;	    
	}
      }
      
      cout << "  It was mapped to CaloTowerDetId's at:" << endl;
      EcalDetIdCountMap::const_iterator ecalDetIdCountMap_iter =  ecalDetIdCountMap.begin();
      for (; ecalDetIdCountMap_iter != ecalDetIdCountMap.end(); ++ecalDetIdCountMap_iter){
	DetId tempDetId = ecalDetIdCountMap_iter -> first;
	CaloTowerDetId tempCaloTowerDetId = ecalDetIdCountMap_iter -> second;
	if (tempDetId == detId )
	  cout << "    ieta = " << caloTowerDetId.ieta() << ", iphi = " << caloTowerDetId.iphi() << endl;
      }  

      valid_multicount_detids++;
      
    }
    else cout << "Entries = " << entries << endl;
  }  
  
  
}

void CaloTowerConstituentsChecker::fillEcalCountingMap(std::vector<CaloTowerDetId> caloTowerDetIds){

  std::vector<CaloTowerDetId>::iterator caloTowerDetId_iter = caloTowerDetIds.begin();
  
  for(; caloTowerDetId_iter != caloTowerDetIds.end(); ++caloTowerDetId_iter){

    std::vector<DetId> DetIds_mapped = theCaloTowerConstituentsMap -> constituentsOf(*caloTowerDetId_iter);
    
    std::vector<DetId>::iterator detId_mapped_iter = DetIds_mapped.begin();
    
    for (; detId_mapped_iter != DetIds_mapped.end(); ++detId_mapped_iter){
      
      ecalDetIdCountMap.insert(EcalDetIdCount(*detId_mapped_iter, *caloTowerDetId_iter));

    }
    
  }

}

void CaloTowerConstituentsChecker::fillHcalCountingMap(std::vector<CaloTowerDetId> caloTowerDetIds){

  std::vector<CaloTowerDetId>::iterator caloTowerDetId_iter = caloTowerDetIds.begin();
  
  for(; caloTowerDetId_iter != caloTowerDetIds.end(); ++caloTowerDetId_iter){

    std::vector<DetId> DetIds_mapped = theCaloTowerConstituentsMap -> constituentsOf(*caloTowerDetId_iter);
    
    std::vector<DetId>::iterator detId_mapped_iter = DetIds_mapped.begin();
    
    for (; detId_mapped_iter != DetIds_mapped.end(); ++detId_mapped_iter){
      
      hcalDetIdCountMap.insert(HcalDetIdCount(*detId_mapped_iter, *caloTowerDetId_iter));

    }
    
  }

}



void CaloTowerConstituentsChecker::checkForCaloTowerMapDisagreement(std::vector<DetId> detIds_original){
  
  bool mapping_successful;
  bool all_mappings_successful = true;

  std::vector<DetId>::iterator detId_original_iter = detIds_original.begin();
  
  for (; detId_original_iter != detIds_original.end(); ++detId_original_iter){

    DetId detId_original = *detId_original_iter;

    CaloTowerDetId caloTowerDetId = theCaloTowerConstituentsMap -> towerOf (*detId_original_iter);
    
    std::vector<DetId> DetIds_mapped = theCaloTowerConstituentsMap -> constituentsOf(caloTowerDetId);
    std::vector<DetId>::iterator detId_mapped_iter = DetIds_mapped.begin();
    
    mapping_successful = false;
    
    for (;  detId_mapped_iter != DetIds_mapped.end(); ++detId_mapped_iter)
      mapping_successful |= (*detId_mapped_iter == detId_original);
    
    all_mappings_successful &= mapping_successful;

    if (!mapping_successful ){

      cout << "-----------------------------------" << endl;

      if ( detId_original.det() == DetId::Ecal){
	if ( detId_original.subdetId() == EcalBarrel) {
	  EBDetId ebDetId = EBDetId(detId_original);
	  cout << "Original ECAL Det ID: " << ebDetId << endl;
	}
	
	if ( detId_original.subdetId() == EcalEndcap) {
	  EEDetId eeDetId = EEDetId(detId_original);
	  cout << "Original ECAL Det ID: " << eeDetId << endl;
	}

	cout << "Mapped to CaloTower at : " << caloTowerDetId << endl;
	cout << "But we were unable to map from the CaloTower back to the EcalDetId" << endl;
	cout << "Mapped ECAL digis included: " << endl;
	
	detId_mapped_iter = DetIds_mapped.begin();
	for (;  detId_mapped_iter != DetIds_mapped.end(); ++detId_mapped_iter){

	  if ( (*detId_mapped_iter).det() != DetId::Ecal) continue;
	  if ( (*detId_mapped_iter).subdetId() == EcalBarrel ){
	    EBDetId mappedEBDetId = EBDetId(*detId_mapped_iter);
	    cout << " " << mappedEBDetId << endl;
	  }
	  
	  if ( (*detId_mapped_iter).subdetId() == EcalEndcap ){
	    EEDetId mappedEEDetId = EEDetId(*detId_mapped_iter);
	    cout << " " << mappedEEDetId << endl;
	  }
	}
      }
      
      if ( detId_original.det() == DetId::Hcal){
	HcalDetId hcalDetId = HcalDetId(detId_original);
	cout << "Original HCAL Det ID: " << hcalDetId << endl;

	cout << "Mapped to CaloTower at : " << caloTowerDetId << endl;
	cout << "But we were unable to map from the CaloTower back to the HcalDetId" << endl;
	cout << "Mapped HCAL digis included: " << endl;
	
	detId_mapped_iter = DetIds_mapped.begin();
	for (;  detId_mapped_iter != DetIds_mapped.end(); ++detId_mapped_iter){
	  if ( (*detId_mapped_iter).det() != DetId::Hcal) continue;
	  HcalDetId mappedHcalDetId = HcalDetId(*detId_mapped_iter);
	  cout << " " << mappedHcalDetId << endl;
	}
      }    
    }
  }

  if (all_mappings_successful) std::cout << "Mapping and reverse mapping agree!" << std::endl;
}

void CaloTowerConstituentsChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;

  setup();

  //-----------------------------------------------------
  // Get vectors of all valid EB/EEDetId's
  //-----------------------------------------------------
  
  edm::ESHandle<CaloGeometry> geometry;
  edm::ESHandle<CaloTowerConstituentsMap> caloTowerConstituentsMap;

  iSetup.get<CaloGeometryRecord >().get(geometry);     
  iSetup.get<IdealGeometryRecord>().get(caloTowerConstituentsMap);
  theCaloTowerConstituentsMap = caloTowerConstituentsMap.product();
  
  std::vector<DetId> ebCells = geometry->getValidDetIds(DetId::Ecal, EcalBarrel );
  std::vector<DetId> eeCells = geometry->getValidDetIds(DetId::Ecal, EcalEndcap );

  std::vector<DetId> hbCells = geometry->getValidDetIds(DetId::Hcal, HcalBarrel );
  std::vector<DetId> heCells = geometry->getValidDetIds(DetId::Hcal, HcalEndcap );
  std::vector<DetId> hoCells = geometry->getValidDetIds(DetId::Hcal, HcalOuter  );
  std::vector<DetId> hfCells = geometry->getValidDetIds(DetId::Hcal, HcalForward);
  
  std::cout << "I have " << ebCells.size() << " ECAL barrel  cells" << std::endl;
  std::cout << "I have " << eeCells.size() << " ECAL endcap  cells" << std::endl;
  std::cout << "I have " << hbCells.size() << " HCAL barrel  cells" << std::endl;
  std::cout << "I have " << heCells.size() << " HCAL endcap  cells" << std::endl;
  std::cout << "I have " << hoCells.size() << " HCAL outer   cells" << std::endl;
  std::cout << "I have " << hfCells.size() << " HCAL forward cells" << std::endl;
  
  //-----------------------------------------------------
  // You can't do a similar method for the CaloTowers
  //-----------------------------------------------------
  
  std::vector<CaloTowerDetId> caloTowerDetIds = getCaloTowerDetIds (ebCells,eeCells,hbCells,heCells,hoCells,hfCells);
  
  //-----------------------------------------------------
  // Make sure that for every EB/EEDetId that gets mapped
  // to a CaloTower, your method can map back from
  // that CaloTower to the original EB/EEDetId
  //-----------------------------------------------------
  
  checkForCaloTowerMapDisagreement(ebCells);
  checkForCaloTowerMapDisagreement(eeCells);
  
  checkForCaloTowerMapDisagreement(hbCells);
  checkForCaloTowerMapDisagreement(heCells);
  checkForCaloTowerMapDisagreement(hoCells);
  checkForCaloTowerMapDisagreement(hfCells);

  //-----------------------------------------------------
  // Check for double-counting.  Make sure that no 
  // single EB/EEDetId is a constituent of more than
  // one CaloTower.
  //-----------------------------------------------------
  
  fillEcalCountingMap (caloTowerDetIds);
  fillHcalCountingMap (caloTowerDetIds);

  if (check_eb_doublecount) checkForEcalDoubleCounting (ebCells);
  if (check_ee_doublecount) checkForEcalDoubleCounting (eeCells);

  if (check_eb_doublecount) checkForHcalDoubleCounting (hbCells);
  if (check_ee_doublecount) checkForHcalDoubleCounting (heCells);
  if (check_eb_doublecount) checkForHcalDoubleCounting (hoCells);
  if (check_ee_doublecount) checkForHcalDoubleCounting (hfCells);

  //-----------------------------------------------------
  // Print out our findings
  //-----------------------------------------------------

  cout << endl;
  cout << "There were a total of " << caloTowerDetIds.size() << " unique CaloTowers " << endl;   
  cout << "Out of " << valid_total_detids       << " total valid DetId's" << endl;
  cout << "  "      << valid_unassigned_detids  << " are not assigned to an CaloTower" << endl;
  cout << "  "      << valid_singlecount_detids << " are assigned to one CaloTower by the detIds function" << endl;
  cout << "  "      << valid_multicount_detids  << " are assigned to more than two CaloTowers" << endl;

}

void CaloTowerConstituentsChecker::beginJob(const edm::EventSetup&) {}

void CaloTowerConstituentsChecker::endJob() {

  plotFile -> cd();
  
  for (int depth = 1; depth <= 4; depth++)
    h_doublecount_occ[depth] -> Write(h_doublecount_occ[depth] -> GetName());

  plotFile -> Close();

}

DEFINE_FWK_MODULE(CaloTowerConstituentsChecker);
