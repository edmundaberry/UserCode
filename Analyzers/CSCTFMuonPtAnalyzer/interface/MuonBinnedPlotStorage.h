#ifndef ANALYZERS_CSCTFMUONPTANALYZER_MUONBINNEDPLOTSTORAGE_H
#define ANALYZERS_CSCTFMUONPTANALYZER_MUONBINNEDPLOTSTORAGE_H

#include <vector>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"

class MuonBinnedPlotStorage {
 public:
  MuonBinnedPlotStorage(int nHitCombos, int nFRBCombos, int nPtBins, int nEtaBins, const std::string & fileName);
  ~MuonBinnedPlotStorage();

  void fill(int hitCombo, int frbCombo, int ptBin, int etaBin, int isFirstCombo, int deltaPhi);
  void save();

 private:
  
  const int m_nHitCombos;
  const int m_nFRBCombos;
  const int m_nPtBins;
  const int m_nEtaBins;
  const std::string m_fileName;

  std::vector <std::vector< std::vector < std::vector< std::vector <TH1F*> > > > > m_hists;
};

#endif
