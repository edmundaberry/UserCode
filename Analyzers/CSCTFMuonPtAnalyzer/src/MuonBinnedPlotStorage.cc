#include "Analyzers/CSCTFMuonPtAnalyzer/interface/MuonBinnedPlotStorage.h"
#include "TMath.h"

MuonBinnedPlotStorage::MuonBinnedPlotStorage(int nCombos, int nPtBins, const char* fileName):
  m_nCombos(nCombos),
  m_nPtBins(nPtBins),  
  m_fileName(fileName),
  m_hists ( std::vector< std::vector<TH1F*> > (nCombos+1, std::vector<TH1F*>(nPtBins, new TH1F())))
{
  
  char histName[100];

  for (int iPtBin = 0; iPtBin <= m_nPtBins; ++iPtBin){
    for (int iCombo = 0; iCombo <= m_nCombos; ++iCombo){
     
      sprintf(histName,"Combo%d_PtBin%d",iCombo,iPtBin);            
      m_hists[iCombo][iPtBin] = new TH1F(histName,"",100,-TMath::Pi(),TMath::Pi());
      
    }
  }
}

MuonBinnedPlotStorage::~MuonBinnedPlotStorage(){}

void MuonBinnedPlotStorage::fill (int iCombo, int iPtBin, float deltaPhi) {
  m_hists[iCombo][iPtBin] -> Fill(deltaPhi);
}

void MuonBinnedPlotStorage::save (){

  TFile *file = new TFile(m_fileName,"RECREATE");

  file -> cd();

  for (int iCombo = 0; iCombo <= m_nCombos - 1; ++iCombo){
    for (int iPtBin = 0; iPtBin <= m_nPtBins; ++iPtBin){
      m_hists[iCombo][iPtBin] -> Write();
      delete m_hists[iCombo][iPtBin];
    }
  }

  file -> Close();

}
