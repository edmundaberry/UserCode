#include <iostream>
#include "TFile.h"
#include "TH1F.h"

const char* FILE_NAME = "/afs/cern.ch/user/e/eberry/MuonMomentumStudy/CMSSW_3_2_4/src/Analyzers/CSCTFMuonPtAnalyzer/test/test.root";

const int MIN_PTBIN = 0;
const int MAX_PTBIN = 10;

const int MIN_COMBO_ID = 0;
const int MAX_COMBO_ID = 5;

void find_filled_hists(){
  
  TFile *file = new TFile (FILE_NAME);

  char hist_name[100];

  for (int iPtBin = MIN_PTBIN; iPtBin <= MAX_PTBIN; iPtBin++){
    for (int iComboId = MIN_COMBO_ID; iComboId <= MAX_COMBO_ID; iComboId++){
      
      sprintf(hist_name,"Combo%d_PtBin%d",iComboId, iPtBin);

      TH1F* hist = (TH1F*) file -> Get(hist_name);

      int entries = (int) hist -> GetEntries();

      if ( entries != 0 ){
	std::cout << "Hist with PtBin " << iPtBin 
		  << " and ComboId " << iComboId 
		  << " has " << entries << " entries " << std::endl;
      }

    }    
  }

}
