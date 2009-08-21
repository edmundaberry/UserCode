#ifndef ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONPTANALYZER_H
#define ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONPTANALYZER_H

#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "Analyzers/CSCTFMuonPtAnalyzer/interface/MuonBinnedPlotStorage.h"

class CSCTFMuonPtAnalyzer : public edm::EDAnalyzer {
public:
  explicit CSCTFMuonPtAnalyzer(const edm::ParameterSet&);
  ~CSCTFMuonPtAnalyzer();
  
  
private:

  //--------------------------------------------------
  // Standard analysis functions
  //--------------------------------------------------
  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  //--------------------------------------------------
  // My analysis functions
  //--------------------------------------------------
  
  void analyzeGenMuons();
  void analyzeDataMuons();
  void analyzeCSCTrigPrims();

  //--------------------------------------------------
  // Event and setup pointers
  //--------------------------------------------------

  const edm::Event* m_event;
  const edm::EventSetup* m_setup;

  //--------------------------------------------------
  // Muon information storage struct
  //--------------------------------------------------
  
  struct RealMuon {
    
    RealMuon(double temp_eta, double temp_phi, double temp_pt, int temp_ptBin );
    double eta, phi, pt;
    int ptBin;

  };
  
  void fillHist ( const std::vector<float> & station_average_phi,
		  const std::vector<float> & station_average_eta,
		  const int comboId, const int ptBin );
  
  int reducedComboId ( const int stationId1, const int stationId2);

  void int2bin (uint32_t val, char* string);
  
  //--------------------------------------------------
  // Build objects in the constructor
  //--------------------------------------------------
  
  void buildSRLUTs();
  void buildPtBins();

  //--------------------------------------------------
  // Binning functions
  //--------------------------------------------------

  int getPtBin( double pt );

  //--------------------------------------------------
  // Get the global eta and phi coordinates of a digi
  //--------------------------------------------------

  void getHitCoordinates( const CSCDetId& digi,  const CSCCorrelatedLCTDigi& digi, float& eta, float& phi );
  
  //--------------------------------------------------
  // List of "real" muons
  //--------------------------------------------------

  std::vector<RealMuon> m_realMuonList;
      
  //--------------------------------------------------
  // edm::InputTag's
  //--------------------------------------------------

  edm::InputTag m_cscCorrLCTDigiTag;
  edm::InputTag m_genParticlesTag;

  //--------------------------------------------------
  // Which muons to analyze? (User-set)
  //--------------------------------------------------
  
  const bool m_analyzeGenMuons;
  const bool m_analyzeDataMuons;

  //--------------------------------------------------
  // Binning info
  //--------------------------------------------------

  const int m_nPtBins;
  const double m_minBinnablePt;
  const double m_maxBinnablePt;
  std::vector<double> m_ptBins;

  //--------------------------------------------------
  // Plot Storage
  //--------------------------------------------------
  
  MuonBinnedPlotStorage m_plotStorage;

  //--------------------------------------------------
  // CSC Sector-receiver LUT's
  //--------------------------------------------------

  CSCSectorReceiverLUT * m_srLUTs[5][2];

  //--------------------------------------------------
  // Assume all sectors are identical,
  // So we only declare one sector variable
  //--------------------------------------------------

  const int m_sector;

  //--------------------------------------------------
  // Trigger Mother Board (TMB) firmware variable
  //--------------------------------------------------
  
  const bool m_TMB07;

};

#endif
