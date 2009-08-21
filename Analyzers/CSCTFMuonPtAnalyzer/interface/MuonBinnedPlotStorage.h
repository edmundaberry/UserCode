#ifndef ANALYZERS_CSCTFMUONPTANALYZER_MUONBINNEDPLOTSTORAGE_H
#define ANALYZERS_CSCTFMUONPTANALYZER_MUONBINNEDPLOTSTORAGE_H

#include <vector>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"

class MuonBinnedPlotStorage {
 public:
  MuonBinnedPlotStorage(int nCombos, int nPtBins, const char* fileName);
  ~MuonBinnedPlotStorage();

  void fill(int combo, int ptBin, float deltaPhi);
  void save();

 private:
  
  const int m_nCombos;
  const int m_nPtBins;
  const char* m_fileName;
  
  std::vector< std::vector <TH1F*> > m_hists;
};

#endif
