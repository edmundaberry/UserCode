// -*- C++ -*-
//
// Package:    HcalTrigTowerGeometryChecker
// Class:      HcalTrigTowerGeometryChecker
// 
/**\class HcalTrigTowerGeometryChecker HcalTrigTowerGeometryChecker.cc Analyzers/HcalTrigTowerGeometryChecker/src/HcalTrigTowerGeometryChecker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Edmund Berry"
//         Created:  Thu Mar 19 16:42:49 CDT 2009
// $Id$
//


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
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"

using namespace std;
using namespace edm;

class HcalTrigTowerGeometryChecker : public edm::EDAnalyzer {
public:
  explicit HcalTrigTowerGeometryChecker(const edm::ParameterSet&);
  ~HcalTrigTowerGeometryChecker();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void checkForMapDisagreement(std::vector<DetId> cells);
  void checkForDoubleCounting (std::vector<DetId> cells);
  void fillCountingMap(std::vector<HcalTrigTowerDetId> trigTowerIds);
  std::vector<HcalTrigTowerDetId> getTrigTowerDetIds(std::vector<DetId> hbCells, std::vector<DetId> heCells,
						     std::vector<DetId> hoCells, std::vector<DetId> hfCells);
  void getTrigTowerDetIds_fromSubdet(std::vector<DetId> cellIds, 
				     std::vector<HcalTrigTowerDetId> & towerIds);
  void setup();
  HcalTrigTowerGeometry theTrigTowerGeometry;


  typedef multimap<HcalTrigTowerDetId, int, less<HcalTrigTowerDetId> > HcalTrigTowerDetIdCountMap;
  typedef multimap<HcalTrigTowerDetId, int, less<HcalTrigTowerDetId> >::value_type HcalTrigTowerDetIdCount;
  HcalTrigTowerDetIdCountMap hcalTrigTowerDetIdCountMap;
  

  typedef multimap<HcalDetId, int, less<HcalDetId> > HcalDetIdCountMap;
  typedef multimap<HcalDetId, int, less<HcalDetId> >::value_type HcalDetIdCount;
  HcalDetIdCountMap hcalDetIdCountMap;

  int valid_unassigned_detids;
  int valid_multicount_detids;
  int valid_ok_detids;
  int valid_total_detids;

  bool verbose;
  bool check_hb_doublecount;
  bool check_he_doublecount;
  bool check_ho_doublecount;
  bool check_hf_doublecount;
  string filename;

  TFile *plotFile;
  TH2F  *h_doublecount_occ[5];

};

HcalTrigTowerGeometryChecker::HcalTrigTowerGeometryChecker(const edm::ParameterSet& iConfig){

  verbose              = iConfig.getUntrackedParameter<bool>  ("verbose"             ,true);
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

HcalTrigTowerGeometryChecker::~HcalTrigTowerGeometryChecker() {}

void HcalTrigTowerGeometryChecker::setup(){

  valid_unassigned_detids = 0;
  valid_multicount_detids = 0;
  valid_ok_detids         = 0;
  valid_total_detids      = 0;

}

void HcalTrigTowerGeometryChecker::getTrigTowerDetIds_fromSubdet(std::vector<DetId> cellIds, 
								 std::vector<HcalTrigTowerDetId> & towerIds){

  int entries;
  
  std::vector<DetId>::iterator cellId_iter = cellIds.begin();
  
  for (; cellId_iter != cellIds.end(); ++cellId_iter){
    HcalDetId hcalDetId = HcalDetId(*cellId_iter);
    std::vector<HcalTrigTowerDetId> temp = theTrigTowerGeometry.towerIds(*cellId_iter);
    std::vector<HcalTrigTowerDetId>::iterator temp_iter = temp.begin();
    for (; temp_iter != temp.end(); ++temp_iter){
      
      entries = hcalTrigTowerDetIdCountMap.count(*temp_iter);
      
      if (entries == 0) towerIds.push_back(*temp_iter);
      hcalTrigTowerDetIdCountMap.insert(HcalTrigTowerDetIdCount(*temp_iter,entries+1));
	
    }
  }
  

}

std::vector<HcalTrigTowerDetId> HcalTrigTowerGeometryChecker::getTrigTowerDetIds(std::vector<DetId> hbCells,
										 std::vector<DetId> heCells,
										 std::vector<DetId> hoCells,
										 std::vector<DetId> hfCells){

  std::vector<HcalTrigTowerDetId> retval;

  if (check_hb_doublecount) getTrigTowerDetIds_fromSubdet(hbCells, retval);
  if (check_he_doublecount) getTrigTowerDetIds_fromSubdet(heCells, retval);
  if (check_ho_doublecount) getTrigTowerDetIds_fromSubdet(hoCells, retval);
  if (check_hf_doublecount) getTrigTowerDetIds_fromSubdet(hfCells, retval);
  
  std::sort   (retval.begin(), retval.end());  
  std::unique (retval.begin(), retval.end());

  return retval;

}

void HcalTrigTowerGeometryChecker::checkForDoubleCounting (std::vector<DetId> detIds){

  int cell_ieta, cell_iphi, cell_depth, cell_subdet;
  int entries;
  std::vector<DetId>::iterator detId_iter = detIds.begin();

  for (; detId_iter != detIds.end(); ++detId_iter){

    valid_total_detids++;

    HcalDetId hcalDetId = HcalDetId(*detId_iter);
    
    entries = hcalDetIdCountMap.count(hcalDetId);
    
    if      (entries == 0)  valid_unassigned_detids++;
    else if (entries == 1)  valid_ok_detids++;
    else if (entries >  1){

      cell_ieta   = (int) hcalDetId.ieta();
      cell_iphi   = (int) hcalDetId.iphi();
      cell_depth  = (int) hcalDetId.depth();
      cell_subdet = (int) hcalDetId.subdet();

      h_doublecount_occ[cell_depth] -> Fill(cell_ieta, cell_iphi);

      valid_multicount_detids++;
      
      if (verbose){
	cout << "This HCAL Det ID showed up " << entries << " times" <<  endl
	     << "  subdet = " << cell_subdet 
	     << ", ieta  = "  << cell_ieta 
	     << ", iphi = "   << cell_iphi 
	     << ", depth = "  << cell_depth  << endl;
      }
    }
    else cout << "Entries = " << entries << endl;
  }  
  
  
}

void HcalTrigTowerGeometryChecker::fillCountingMap(std::vector<HcalTrigTowerDetId> trigTowerDetIds){

  int entries;

  std::vector<HcalTrigTowerDetId>::iterator trigTowerDetId_iter = trigTowerDetIds.begin();

  for(; trigTowerDetId_iter != trigTowerDetIds.end(); ++trigTowerDetId_iter){

    std::vector<HcalDetId> hcalDetIds_mapped = theTrigTowerGeometry.detIds(*trigTowerDetId_iter);
    
    std::vector<HcalDetId>::iterator hcalDetId_mapped_iter = hcalDetIds_mapped.begin();
    
    for (; hcalDetId_mapped_iter != hcalDetIds_mapped.end(); ++hcalDetId_mapped_iter){
      
      entries = hcalDetIdCountMap.count(*hcalDetId_mapped_iter);
      
      hcalDetIdCountMap.insert(HcalDetIdCount(*hcalDetId_mapped_iter, entries + 1));

    }
    
  }

}

void HcalTrigTowerGeometryChecker::checkForMapDisagreement(std::vector<DetId> detIds_original){
  
  bool mapping_successful;

  std::vector<DetId>::iterator detId_original_iter = detIds_original.begin();

  for (; detId_original_iter != detIds_original.end(); ++detId_original_iter){

    const HcalDetId hcalDetId_original = HcalDetId(*detId_original_iter);

    std::vector<HcalTrigTowerDetId> trigTowerDetIds = theTrigTowerGeometry.towerIds(hcalDetId_original);
    std::vector<HcalTrigTowerDetId>::iterator trigTowerDetId_iter = trigTowerDetIds.begin();

    for (; trigTowerDetId_iter != trigTowerDetIds.end(); ++trigTowerDetId_iter){

      std::vector<HcalDetId> hcalDetIds_mapped = theTrigTowerGeometry.detIds(*trigTowerDetId_iter);

      std::vector<HcalDetId>::iterator hcalDetId_mapped_iter = hcalDetIds_mapped.begin();

      mapping_successful = false;
 
      for (;  hcalDetId_mapped_iter != hcalDetIds_mapped.end(); ++hcalDetId_mapped_iter)
	mapping_successful |= (*hcalDetId_mapped_iter == hcalDetId_original);

      if (!mapping_successful) {
	cout << "-----------------------------------" << endl;
	cout << "Original HCAL Det ID:" << endl 
	     << "  subdet = " << hcalDetId_original.subdet() 
	     << ", ieta  = "  << hcalDetId_original.ieta() 
	     << ", iphi = "   << hcalDetId_original.iphi() 
	     << ", depth = "  << hcalDetId_original.depth() << endl
	     << "Mapped to CaloTower at: " << endl
	     << "  ieta = " << (*trigTowerDetId_iter).ieta() 
	     << ", iphi = " << (*trigTowerDetId_iter).iphi() << endl
	     << "But we were unable to map from the CaloTower back to the HcalDetId" << endl
	     << "Mapped HCAL digis included: " << endl;
	
	hcalDetId_mapped_iter = hcalDetIds_mapped.begin();
	for (;  hcalDetId_mapped_iter != hcalDetIds_mapped.end(); ++hcalDetId_mapped_iter){
	  cout << "  subdet = " << (*hcalDetId_mapped_iter).subdet() 
	       << ", ieta  = "  << (*hcalDetId_mapped_iter).ieta() 
	       << ", iphi = "   << (*hcalDetId_mapped_iter).iphi() 
	       << ", depth = "  << (*hcalDetId_mapped_iter).depth() << endl;
	}
      }
    }

  }

}

void HcalTrigTowerGeometryChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;

  setup();

  //-----------------------------------------------------
  // Get vectors of all valid HcalDetId's
  //-----------------------------------------------------
  
  edm::ESHandle<CaloGeometry> geometry;
  iSetup.get<CaloGeometryRecord>().get(geometry);     
  
  std::vector<DetId> hbCells    = geometry->getValidDetIds(DetId::Hcal, HcalBarrel      );
  std::vector<DetId> heCells    = geometry->getValidDetIds(DetId::Hcal, HcalEndcap      );
  std::vector<DetId> hoCells    = geometry->getValidDetIds(DetId::Hcal, HcalOuter       );
  std::vector<DetId> hfCells    = geometry->getValidDetIds(DetId::Hcal, HcalForward     );
  
  //-----------------------------------------------------
  // You can't do a similar method for the 
  // HcalTrigTowerDetId's.
  // Get them yourself from the HcalDetId's by using
  // HcalTrigTowerGeometry::towerIds
  //-----------------------------------------------------

  // This doesn't work:
  // std::vector<DetId> trigTowerDetIds = geometry->getValidDetIds(DetId::Hcal, HcalTriggerTower);   
  std::vector<HcalTrigTowerDetId> trigTowerDetIds = getTrigTowerDetIds(hbCells, heCells, hoCells, hfCells);
  
  //-----------------------------------------------------
  // Make sure that for every HcalDetId that gets mapped
  // to a trigger tower with HcalTrigTowerGeometry::towerIds,
  // your method can map back from that trigger tower to 
  // the original HcalDetId
  //-----------------------------------------------------
  
  checkForMapDisagreement(hbCells);
  checkForMapDisagreement(heCells);
  checkForMapDisagreement(hoCells);
  checkForMapDisagreement(hfCells);

  //-----------------------------------------------------
  // Check for double-counting.  Make sure that no 
  // single HcalDetId is a constituent of more than
  // one TriggerTower.
  //
  // EXCEPT: HE cells with 21 <= |ieta| <= 29, where this
  // double-counting is expected.
  //-----------------------------------------------------
  
  fillCountingMap (trigTowerDetIds);

  if (check_hb_doublecount) checkForDoubleCounting (hbCells);
  if (check_he_doublecount) checkForDoubleCounting (heCells);
  if (check_ho_doublecount) checkForDoubleCounting (hoCells);
  if (check_hf_doublecount) checkForDoubleCounting (hfCells);
  
  //-----------------------------------------------------
  // Print out our findings
  //-----------------------------------------------------
  cout << endl;
  cout << "There were a total of " << trigTowerDetIds.size() << " unique trigger towers " << endl;   
  cout << "Out of " << valid_total_detids       << " total valid HcalDetId's" << endl;
  cout << "  "      << valid_unassigned_detids  << " are not assigned to an HcalTrigTower" << endl;
  cout << "  "      << valid_ok_detids          << " are assigned to one HcalTrigTower" << endl;
  cout << "  "      << valid_multicount_detids  << " are assigned to multiple HcalTrigTowers" << endl;
  if (check_he_doublecount){
    cout << "We expect " << 1440 << " instances of double-counting."  << endl;
    cout << "This is due to HE cells with |ieta| between 21 and 29 "  << endl;
  }
}

void HcalTrigTowerGeometryChecker::beginJob(const edm::EventSetup&) {}

void HcalTrigTowerGeometryChecker::endJob() {

  plotFile -> cd();
  
  for (int depth = 1; depth <= 4; depth++)
    h_doublecount_occ[depth] -> Write(h_doublecount_occ[depth] -> GetName());

  plotFile -> Close();

}

DEFINE_FWK_MODULE(HcalTrigTowerGeometryChecker);
