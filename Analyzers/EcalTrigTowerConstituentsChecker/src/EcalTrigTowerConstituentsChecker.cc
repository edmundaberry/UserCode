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
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

using namespace std;
using namespace edm;

class EcalTrigTowerConstituentsChecker : public edm::EDAnalyzer {
public:
  explicit EcalTrigTowerConstituentsChecker(const edm::ParameterSet&);
  ~EcalTrigTowerConstituentsChecker();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void checkForEcalMapDisagreement(std::vector<DetId> cells);
  void checkForEcalDoubleCounting (std::vector<DetId> cells);

  void fillEcalCountingMap(std::vector<EcalTrigTowerDetId> trigTowerIds);

  void getEcalTrigTowerDetIds_fromSubdet(const std::vector<DetId> & cellIds, 
					 std::vector<EcalTrigTowerDetId> & towerIds);
  std::vector<EcalTrigTowerDetId> getEcalTrigTowerDetIds(std::vector<DetId> ebCells, std::vector<DetId> eeCells);
  
  void setup();

  typedef multimap<uint32_t, int, less<uint32_t> > EcalTrigTowerDetIdCountMap;
  typedef multimap<uint32_t, int, less<uint32_t> >::value_type EcalTrigTowerDetIdCount;
  EcalTrigTowerDetIdCountMap ecalTrigTowerDetIdCountMap;
  
  typedef multimap<DetId, EcalTrigTowerDetId, less<DetId> > EcalDetIdCountMap;
  typedef multimap<DetId, EcalTrigTowerDetId, less<DetId> >::value_type EcalDetIdCount;
  EcalDetIdCountMap ecalDetIdCountMap;

  const EcalTrigTowerConstituentsMap * theEcalTrigTowerConstituentsMap;
  edm::ESHandle<EcalTrigTowerConstituentsMap> ecalTrigTowerConstituentsMap_Handle;

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

  string filename;

  TFile *plotFile;
  TH2F  *h_doublecount_occ[5];

};

EcalTrigTowerConstituentsChecker::EcalTrigTowerConstituentsChecker(const edm::ParameterSet& iConfig){

  verbose              = iConfig.getUntrackedParameter<bool>  ("verbose"             ,true);

  check_eb_doublecount = iConfig.getUntrackedParameter<bool>  ("check_eb_doublecount",true);
  check_ee_doublecount = iConfig.getUntrackedParameter<bool>  ("check_ee_doublecount",true);

  filename             = iConfig.getUntrackedParameter<string>("filename", "double_count.root");

  plotFile = new TFile(filename.c_str(), "RECREATE");

  char depthName[500];

  for (int depth = 1; depth <= 4; depth++){
    sprintf(depthName,"h_doublecount_occ_depth%d",depth);
    h_doublecount_occ[depth] = new TH2F(depthName,"",83,-41.5,41.5,72,0.5,72.5);
  }


}

EcalTrigTowerConstituentsChecker::~EcalTrigTowerConstituentsChecker() {}

void EcalTrigTowerConstituentsChecker::setup(){

  valid_unassigned_detids     = 0;
  valid_multicount_detids     = 0;
  valid_doublecount_detids    = 0;
  valid_doublecount_ok_detids = 0;
  valid_singlecount_detids    = 0;
  valid_singlecount_ok_detids = 0;
  valid_total_detids          = 0;

}

void EcalTrigTowerConstituentsChecker::getEcalTrigTowerDetIds_fromSubdet(const std::vector<DetId> & cellIds, 
									 std::vector<EcalTrigTowerDetId> & towerIds){

  int entries;
  
  std::vector<DetId>::const_iterator cellId_iter = cellIds.begin();
  std::vector<DetId>::const_iterator cellId_end  = cellIds.end();

  for (; cellId_iter != cellId_end; ++cellId_iter){

    EcalTrigTowerDetId temp = theEcalTrigTowerConstituentsMap -> towerOf(*cellId_iter);       

    if (temp.ietaAbs() == 28 || temp.ietaAbs() == 27){
      
      EcalSubdetector subdet = temp.subDet();
      int ieta               = temp.ieta();
      int ietaAbs            = temp.ietaAbs();
      int iphi               = temp.iphi();
      int zside              = ieta / ietaAbs;
            
      bool iphi_is_even = ( iphi == (iphi/2)*2);

      if ( iphi_is_even) --iphi;
      if (!iphi_is_even) ++iphi;
  
      EcalTrigTowerDetId pseudo = EcalTrigTowerDetId(zside,subdet,ietaAbs,iphi);
      
      if ( ! (theEcalTrigTowerConstituentsMap -> constituentsOf(pseudo).empty())){
	std::cout << "NOT A PSEUDO TOWER!!!" << std::endl;
      }

      else {
      
	int pseudo_entries = ecalTrigTowerDetIdCountMap.count(pseudo.rawId());
	
	if (pseudo_entries == 0) {
	  towerIds.push_back(pseudo);
	}
	
	ecalTrigTowerDetIdCountMap.insert(EcalTrigTowerDetIdCount(pseudo.rawId(),pseudo_entries+1));
      }
      
    }

    entries = ecalTrigTowerDetIdCountMap.count(temp.rawId());
    
    if (entries == 0) towerIds.push_back(temp);
    
    ecalTrigTowerDetIdCountMap.insert(EcalTrigTowerDetIdCount(temp.rawId(),entries+1));

  }
}

std::vector<EcalTrigTowerDetId> EcalTrigTowerConstituentsChecker::getEcalTrigTowerDetIds(std::vector<DetId> ebCells,
											 std::vector<DetId> eeCells){

  std::vector<EcalTrigTowerDetId> retval;


  if (check_eb_doublecount) getEcalTrigTowerDetIds_fromSubdet(ebCells, retval);
  if (check_ee_doublecount) getEcalTrigTowerDetIds_fromSubdet(eeCells, retval);

  std::sort   (retval.begin(), retval.end());  
  std::unique (retval.begin(), retval.end());

  return retval;

}

void EcalTrigTowerConstituentsChecker::checkForEcalDoubleCounting (std::vector<DetId> detIds){

  int entries;
  std::vector<DetId>::iterator detId_iter = detIds.begin();

  for (; detId_iter != detIds.end(); ++detId_iter){

    valid_total_detids++;

    DetId detId = DetId(*detId_iter);
    
    EcalTrigTowerDetId trigTowerDetId = theEcalTrigTowerConstituentsMap -> towerOf(detId);

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
      
      cout << "  It was mapped to EcalTrigTowerDetId's at:" << endl;
      EcalDetIdCountMap::const_iterator ecalDetIdCountMap_iter =  ecalDetIdCountMap.begin();
      for (; ecalDetIdCountMap_iter != ecalDetIdCountMap.end(); ++ecalDetIdCountMap_iter){
	DetId tempDetId = ecalDetIdCountMap_iter -> first;
	EcalTrigTowerDetId tempEcalTrigTowerDetId = ecalDetIdCountMap_iter -> second;
	if (tempDetId == detId )
	  cout << "    ieta = " << trigTowerDetId.ieta() << ", iphi = " << trigTowerDetId.iphi() << endl;
      }
      
      valid_multicount_detids++;
      
    }
    else cout << "Entries = " << entries << endl;
  }  
  
  
}

void EcalTrigTowerConstituentsChecker::fillEcalCountingMap(std::vector<EcalTrigTowerDetId> trigTowerDetIds){

  //int entries;

  std::vector<EcalTrigTowerDetId>::iterator trigTowerDetId_iter = trigTowerDetIds.begin();

  for(; trigTowerDetId_iter != trigTowerDetIds.end(); ++trigTowerDetId_iter){

    std::vector<DetId> DetIds_mapped = theEcalTrigTowerConstituentsMap -> constituentsOf(*trigTowerDetId_iter);
    
    std::vector<DetId>::iterator detId_mapped_iter = DetIds_mapped.begin();
    
    for (; detId_mapped_iter != DetIds_mapped.end(); ++detId_mapped_iter){
      
      ecalDetIdCountMap.insert(EcalDetIdCount(*detId_mapped_iter, *trigTowerDetId_iter));

    }
    
  }

}

void EcalTrigTowerConstituentsChecker::checkForEcalMapDisagreement(std::vector<DetId> detIds_original){
  
  bool mapping_successful;

  std::vector<DetId>::iterator detId_original_iter = detIds_original.begin();
  
  for (; detId_original_iter != detIds_original.end(); ++detId_original_iter){

    DetId detId_original = *detId_original_iter;

    EcalTrigTowerDetId trigTowerDetId = theEcalTrigTowerConstituentsMap -> towerOf (*detId_original_iter);
    
    std::vector<DetId> DetIds_mapped = theEcalTrigTowerConstituentsMap -> constituentsOf(trigTowerDetId);
    std::vector<DetId>::iterator detId_mapped_iter = DetIds_mapped.begin();
    
    mapping_successful = false;
    
    for (;  detId_mapped_iter != DetIds_mapped.end(); ++detId_mapped_iter)
      mapping_successful |= (*detId_mapped_iter == detId_original);
    
    if (!mapping_successful ){

      cout << "-----------------------------------" << endl;

      if ( detId_original.subdetId() == EcalBarrel) {
	EBDetId ebDetId = EBDetId(detId_original);
	cout << "Original ECAL Det ID: " << ebDetId << endl;
      }

      if ( detId_original.subdetId() == EcalEndcap) {
	EEDetId eeDetId = EEDetId(detId_original);
	cout << "Original ECAL Det ID: " << eeDetId << endl;
      }

      cout << "Mapped to Trigger Tower at : " << trigTowerDetId << endl;
      cout << "But we were unable to map from the Trigger Tower back to the EcalDetId" << endl;
      cout << "Mapped ECAL digis included: " << endl;
      
      detId_mapped_iter = DetIds_mapped.begin();
      for (;  detId_mapped_iter != DetIds_mapped.end(); ++detId_mapped_iter){
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
  }
}

void EcalTrigTowerConstituentsChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;

  setup();

  //-----------------------------------------------------
  // Get vectors of all valid EB/EEDetId's
  //-----------------------------------------------------
  
  edm::ESHandle<CaloGeometry> geometry;
  iSetup.get<CaloGeometryRecord >().get(geometry);     
  iSetup.get<IdealGeometryRecord>().get(ecalTrigTowerConstituentsMap_Handle);
  
  theEcalTrigTowerConstituentsMap = ecalTrigTowerConstituentsMap_Handle.product();
  
  std::vector<DetId> ebCells = geometry->getValidDetIds(DetId::Ecal, EcalBarrel );
  std::vector<DetId> eeCells = geometry->getValidDetIds(DetId::Ecal, EcalEndcap );

  std::cout << "I have " << ebCells.size() << " ECAL barrel cells" << std::endl;
  std::cout << "I have " << eeCells.size() << " ECAL endcap cells" << std::endl;

  //-----------------------------------------------------
  // You can't do a similar method for the 
  // EcalTrigTowerDetId's.
  //-----------------------------------------------------
  
  std::vector<EcalTrigTowerDetId> ecalTrigTowerDetIds = getEcalTrigTowerDetIds (ebCells, eeCells);

  //-----------------------------------------------------
  // Make sure that for every EB/EEDetId that gets mapped
  // to a trigger tower, your method can map back from
  // that trigger tower to the original EB/EEDetId
  //-----------------------------------------------------
  
  checkForEcalMapDisagreement(ebCells);
  checkForEcalMapDisagreement(eeCells);

  //-----------------------------------------------------
  // Check for double-counting.  Make sure that no 
  // single EB/EEDetId is a constituent of more than
  // one TriggerTower.
  //-----------------------------------------------------
  
  fillEcalCountingMap (ecalTrigTowerDetIds);

  if (check_eb_doublecount) checkForEcalDoubleCounting (ebCells);
  if (check_ee_doublecount) checkForEcalDoubleCounting (eeCells);

  //-----------------------------------------------------
  // Print out our findings
  //-----------------------------------------------------

  cout << endl;
  cout << "There were a total of " << ecalTrigTowerDetIds.size() << " unique trigger towers " << endl;   
  cout << "Out of " << valid_total_detids       << " total valid EB/EEDetId's" << endl;
  cout << "  "      << valid_unassigned_detids  << " are not assigned to an EcalTrigTower" << endl;
  cout << "  "      << valid_singlecount_detids << " are assigned to one EcalTrigTower by the detIds function" << endl;
  cout << "  "      << valid_multicount_detids  << " are assigned to more than two EcalTrigTowers" << endl;

}

void EcalTrigTowerConstituentsChecker::beginJob(const edm::EventSetup&) {}

void EcalTrigTowerConstituentsChecker::endJob() {

  plotFile -> cd();
  
  for (int depth = 1; depth <= 4; depth++)
    h_doublecount_occ[depth] -> Write(h_doublecount_occ[depth] -> GetName());

  plotFile -> Close();

}

DEFINE_FWK_MODULE(EcalTrigTowerConstituentsChecker);
