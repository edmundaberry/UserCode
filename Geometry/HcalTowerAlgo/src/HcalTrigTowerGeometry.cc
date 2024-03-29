#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "Geometry/HcalTowerAlgo/src/HcalHardcodeGeometryData.h"

#include <iostream>
#include <cassert>

HcalTrigTowerGeometry::HcalTrigTowerGeometry() {
  useShortFibers_=true;
  useHFQuadPhiRings_=true;
}

void HcalTrigTowerGeometry::setupHF(bool useShortFibers, bool useQuadRings) {
  useShortFibers_=useShortFibers;
  useHFQuadPhiRings_=useQuadRings;
}

std::vector<HcalDetId>
HcalTrigTowerGeometry::detIds(const HcalTrigTowerDetId & hcalTrigTowerDetId) const {
  
  std::vector<HcalDetId> retval;

  int tower_ieta = hcalTrigTowerDetId.ieta();
  int tower_iphi = hcalTrigTowerDetId.iphi();

  int cell_ieta = tower_ieta;
  int cell_iphi = tower_iphi;

  int min_depth, n_depths;

  // HB
  
  if (abs(cell_ieta) <= theTopology.lastHBRing()){
    theTopology.depthBinInformation(HcalBarrel, abs(tower_ieta), n_depths, min_depth);
    for (int cell_depth = min_depth; cell_depth <= min_depth + n_depths - 1; cell_depth++)
      retval.push_back(HcalDetId(HcalBarrel,cell_ieta,cell_iphi,cell_depth));
  }

  // HO
  
  if (abs(cell_ieta) <= theTopology.lastHORing()){ 
    theTopology.depthBinInformation(HcalOuter , abs(tower_ieta), n_depths, min_depth);  
    for (int ho_depth = min_depth; ho_depth <= min_depth + n_depths - 1; ho_depth++)
      retval.push_back(HcalDetId(HcalOuter, cell_ieta,cell_iphi,ho_depth));
  }

  // HE 

  if (abs(cell_ieta) >= theTopology.firstHERing() && 
      abs(cell_ieta) <  theTopology.lastHERing()){   

    theTopology.depthBinInformation(HcalEndcap, abs(tower_ieta), n_depths, min_depth);
    
    // Special for double-phi cells
    if (abs(cell_ieta) >= theTopology.firstHEDoublePhiRing())
      if (tower_iphi%2 == 0) cell_iphi = tower_iphi - 1;
    
    for (int cell_depth = min_depth; cell_depth <= min_depth + n_depths - 1; cell_depth++)
      retval.push_back(HcalDetId(HcalEndcap, cell_ieta, cell_iphi, cell_depth));
    
    // Special for split-eta cells
    if (abs(tower_ieta) == 28){
      theTopology.depthBinInformation(HcalEndcap, abs(tower_ieta)+1, n_depths, min_depth);
      for (int cell_depth = min_depth; cell_depth <= min_depth + n_depths - 1; cell_depth++){
	if (tower_ieta < 0) retval.push_back(HcalDetId(HcalEndcap, tower_ieta - 1, cell_iphi, cell_depth));
	if (tower_ieta > 0) retval.push_back(HcalDetId(HcalEndcap, tower_ieta + 1, cell_iphi, cell_depth));
      }
    }
    
  }
    
  // HF 
  
  if (abs(cell_ieta) >= theTopology.firstHFRing()){  
    
    int HfTowerPhiSize     = 72 / nPhiBins(tower_ieta);
    int HfTowerEtaSize     = hfTowerEtaSize(tower_ieta);
    int FirstHFRingInTower = firstHFRingInTower(abs(tower_ieta));
    
    for (int iHFTowerPhiSegment = 0; iHFTowerPhiSegment < HfTowerPhiSize; iHFTowerPhiSegment++){      
            
      cell_iphi =  (tower_iphi / HfTowerPhiSize) * HfTowerPhiSize; // Find the minimum phi segment
      cell_iphi -= 2;                        // The first trigger tower starts at HCAL iphi = 71, not HCAL iphi = 1
      cell_iphi += iHFTowerPhiSegment;       // Get all of the HCAL iphi values in this trigger tower
      cell_iphi += 72;                       // Don't want to take the mod of a negative number
      cell_iphi =  cell_iphi % 72;           // There are, at most, 72 cells.
      cell_iphi += 1;                        // There is no cell at iphi = 0
      
      if (cell_iphi%2 == 0) continue;        // These cells don't exist.

      for (int iHFTowerEtaSegment = 0; iHFTowerEtaSegment < HfTowerEtaSize; iHFTowerEtaSegment++){
		
	cell_ieta = FirstHFRingInTower + iHFTowerEtaSegment;

	if (cell_ieta >= 40 && cell_iphi%4 == 1) continue;  // These cells don't exist.

	theTopology.depthBinInformation(HcalForward, cell_ieta, n_depths, min_depth);  

	// Negative tower_ieta -> negative cell_ieta
	if (tower_ieta < 0) cell_ieta *= -1;	       

	for (int cell_depth = min_depth; cell_depth <= min_depth + n_depths - 1; cell_depth++)
	  retval.push_back(HcalDetId(HcalForward, cell_ieta, cell_iphi, cell_depth));
	
      }    
    }
  }

  return retval;
}

std::vector<HcalTrigTowerDetId>
HcalTrigTowerGeometry::towerIds(const HcalDetId & cellId) const {

  std::vector<HcalTrigTowerDetId> results;

  if(cellId.subdet() == HcalForward) {
    // short fibers don't count
    if(cellId.depth() == 1 || useShortFibers_) {
      // first do eta
      int hfRing = cellId.ietaAbs();
      int ieta = firstHFTower(); 
      // find the tower that contains this ring
      while(hfRing >= firstHFRingInTower(ieta+1)) {
        ++ieta;
      }

      ieta *= cellId.zside();

      // now for phi
      // HF towers are quad, 18 in phi.
      // go two cells per trigger tower.
      int iphi = (((cellId.iphi()+1)/4) * 4 + 1)%72; // 71+1 --> 1, 3+5 --> 5
      if (useHFQuadPhiRings_ || cellId.ietaAbs() < theTopology.firstHFQuadPhiRing())
        results.push_back( HcalTrigTowerDetId(ieta, iphi) );
    }
      
  } else {
    // the first twenty rings are one-to-one
    if(cellId.ietaAbs() < theTopology.firstHEDoublePhiRing()) {    
      results.push_back( HcalTrigTowerDetId(cellId.ieta(), cellId.iphi()) );
    } else {
      // the remaining rings are two-to-one in phi
      int iphi1 = cellId.iphi();
      int ieta = cellId.ieta();
      // the last eta ring in HE is split.  Recombine.
      if(ieta == theTopology.lastHERing()) --ieta;
      if(ieta == -theTopology.lastHERing()) ++ieta;

      results.push_back( HcalTrigTowerDetId(ieta, iphi1) );
      results.push_back( HcalTrigTowerDetId(ieta, iphi1+1) );
    }
  }

  return results;
}

int HcalTrigTowerGeometry::hfTowerEtaSize(int ieta) const {
  int ietaAbs = abs(ieta); 
  assert(ietaAbs >= firstHFTower() && ietaAbs <= nTowers());
  // the first three come from rings 29-31, 32-34, 35-37. The last has 4 rings: 38-41
  return (ietaAbs == nTowers()) ? 4 : 3;
}


int HcalTrigTowerGeometry::firstHFRingInTower(int ietaTower) const {
  // count up to the correct HF ring
  int inputTower = abs(ietaTower);
  int result = theTopology.firstHFRing();
  for(int iTower = firstHFTower(); iTower != inputTower; ++iTower) {
    result += hfTowerEtaSize(iTower);
  }
  
  // negative in, negative out.
  if(ietaTower < 0) result *= -1;
  return result;
}


void HcalTrigTowerGeometry::towerEtaBounds(int ieta, double & eta1, double & eta2) const {
  int ietaAbs = abs(ieta);
  if(ietaAbs < firstHFTower()) {
    eta1 = theHBHEEtaBounds[ietaAbs-1];
    eta2 = theHBHEEtaBounds[ietaAbs];
    // the last tower is split, so get tower 29, too
    if(ietaAbs == theTopology.lastHERing()-1) {
      eta2 = theHBHEEtaBounds[ietaAbs+1];
    } 
  } else {
    // count from 0
    int hfIndex = firstHFRingInTower(ietaAbs) - theTopology.firstHFRing();
    eta1 = theHFEtaBounds[hfIndex];
    eta2 = theHFEtaBounds[hfIndex + hfTowerEtaSize(ieta)];
  }

  // get the signs and order right
  if(ieta < 0) {
    double tmp = eta1;
    eta1 = -eta2;
    eta2 = -tmp;
  }
}
