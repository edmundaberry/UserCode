#ifndef ANALYZERS_CSCTFMUONPTANALYZER_MUONBINNEDPLOTSTORAGE_H
#define ANALYZERS_CSCTFMUONPTANALYZER_MUONBINNEDPLOTSTORAGE_H

#include <vector>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"

class MuonBinnedPlotStorage {
 public:
  MuonBinnedPlotStorage(int nHitCombos, int nFRBCombos  , int nPtBins, 
			int nEtaBins  , int nQualityBin, 
			const std::string & dphiFileName,
			const std::string & detaFileName,
			const std::string & etaFileName);
  
  ~MuonBinnedPlotStorage();

  void fillDPhi   ( int hitCombo, int frbCombo, int ptBin, int etaBin, int s1ring, int quality, int isFirstCombo, int deltaPhi);
  void fillDEta   ( int hitCombo, int frbCombo, int ptBin, int etaBin, int s1ring, int quality, int isFirstCombo, int deltaEta);
  void fillEtaBin ( int hitCombo, int quality, int ptBin, int isFirstCombo, int etaBin );
  void save();

 private:
  
  const int m_nHitCombos;
  const int m_nFRBCombos;
  const int m_nPtBins;
  const int m_nEtaBins;
  const int m_nQualityBins;
  const std::string m_dphiFileName;
  const std::string m_detaFileName;
  const std::string m_etaFileName;
  
  std::vector < std::vector <std::vector <std::vector< std::vector < std::vector< std::vector <TH1F*> > > > > > > m_dphiHists;
  std::vector < std::vector <std::vector <std::vector< std::vector < std::vector< std::vector <TH1F*> > > > > > > m_detaHists;
  std::vector <std::vector <std::vector< std::vector < TH1F* > > > > m_etaHists;
};

#endif
