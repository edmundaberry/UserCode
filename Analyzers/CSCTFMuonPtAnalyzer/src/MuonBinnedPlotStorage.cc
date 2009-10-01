#include "Analyzers/CSCTFMuonPtAnalyzer/interface/MuonBinnedPlotStorage.h"
#include "TMath.h"
#include <iostream>

MuonBinnedPlotStorage::MuonBinnedPlotStorage(int nHitCombos, int nFRBCombos, int nPtBins, int nEtaBins, const std::string& fileName):
  m_nHitCombos(nHitCombos),
  m_nFRBCombos(nFRBCombos),
  m_nPtBins(nPtBins),  
  m_nEtaBins(nEtaBins),
  m_fileName(fileName),
  m_hists ( std::vector < std::vector <std::vector <std::vector< std::vector<TH1F*> > > > > (nHitCombos+1, std::vector <std::vector< std::vector < std::vector <TH1F*> > > > (nFRBCombos+1, std::vector< std::vector< std::vector <TH1F*> > > (nPtBins +1 , std::vector< std::vector <TH1F*> > ( nEtaBins+1, std::vector < TH1F* > ( 2, new TH1F()))))))
{

  char histName[200];

  for (int iEtaBin = 0; iEtaBin <= m_nEtaBins; iEtaBin++){
    for (int iPtBin = 0; iPtBin <= m_nPtBins; ++iPtBin){
      for (int iHitCombo = 0; iHitCombo <= m_nHitCombos; ++iHitCombo){
	for (int iFRBCombo = 0; iFRBCombo <= m_nFRBCombos; ++iFRBCombo){
	  for (int isFirstCombo = 0; isFirstCombo <= 1; ++isFirstCombo){
	    
	    sprintf(histName,"HitCombo%d_FDBCombo%d_PtBin%d_EtaBin%d_IsFirstCombo%d",
		    iHitCombo,iFRBCombo,iPtBin, iEtaBin, isFirstCombo);            

	    m_hists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][isFirstCombo] = new TH1F(histName,"",210,-525., 525. );
	    
	  }
	}
      }
    }
  }
}

MuonBinnedPlotStorage::~MuonBinnedPlotStorage(){}

void MuonBinnedPlotStorage::fill (int iHitCombo, int iFRBCombo, int iPtBin, int iEtaBin, int isFirstCombo, int deltaPhi) {
  m_hists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][isFirstCombo] -> Fill(deltaPhi);
}

void MuonBinnedPlotStorage::save (){

  TFile *file = new TFile(m_fileName.c_str(),"RECREATE");

  file -> cd();

  for (int iHitCombo = 0; iHitCombo <= m_nHitCombos - 1; ++iHitCombo){
    for (int iFRBCombo = 0; iFRBCombo <= m_nFRBCombos; ++iFRBCombo){
      for (int iPtBin = 0; iPtBin <= m_nPtBins; ++iPtBin){
	for (int iEtaBin = 0; iEtaBin <= m_nEtaBins; iEtaBin++){
	  for (int isFirstCombo = 0; isFirstCombo <= 1; ++isFirstCombo){

	    Double_t entries =  m_hists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][isFirstCombo] -> GetEntries();
	    
	    if (entries > 0.0) m_hists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][isFirstCombo]-> Write();
	    
	    delete m_hists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][isFirstCombo];
	  }
	}
      }
    }
  }

  file -> Close();

}
