#ifndef ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONPTANALYZER_H
#define ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONPTANALYZER_H

// Framework
#include "FWCore/Framework/interface/ESHandle.h"

// LUTs
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"

// My headers
#include "Analyzers/CSCTFMuonPtAnalyzer/interface/MuonBinnedPlotStorage.h"
#include "Analyzers/CSCTFMuonPtAnalyzer/interface/FillCSCTFMuonTree.h"
#include "Analyzers/CSCTFMuonPtAnalyzer/interface/CSCTFMuonTree.h"

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
  // Conditions getter function
  //--------------------------------------------------

  void getConditions();

  //--------------------------------------------------
  // My analysis functions
  //--------------------------------------------------
  
  void analyzeGenMuons();
  void analyzeDataMuons();
  void analyzeCSCTrigPrims();
  void analyzeCSCSimHits();

  //--------------------------------------------------
  // Helper functions
  //--------------------------------------------------
  
  int reducedComboHitId ( const int stationId1, const int stationId2);
  int reducedComboFRBId ( const int frbit1    , const int frbit2    );

  void int2bin (uint32_t val, char* string);
  
  //--------------------------------------------------
  // Build objects in the constructor
  //--------------------------------------------------
  
  void buildSRLUTs();

  //--------------------------------------------------
  // Binning functions
  //--------------------------------------------------

  int getPtBin ( float pt  );
  int getEtaBin( float eta );

  //--------------------------------------------------
  // Get the global eta and phi coordinates of a digi
  //--------------------------------------------------

  void getHitCoordinates (const CSCDetId& detId , const CSCCorrelatedLCTDigi& digi,
			  int& lclPhi_phi_bend, int& lclPhi_phi, int& gblPhi_phi, int& gblEta_eta,
			  float& cms_eta, float& cms_phi, bool& bad_phi );
  
  //--------------------------------------------------
  // trig prim information storage
  //--------------------------------------------------  

  struct CSCTPInfo {
    int station;
    int ring;
    int sector;
    int frBit;    
    int gblPhi;
    float eta;
  };

  void analyzeCombos( const std::vector < std::vector < CSCTPInfo > > & v_tpInfo );

  //--------------------------------------------------
  // Event and setup pointers
  //--------------------------------------------------

  const edm::Event* m_event;
  const edm::EventSetup* m_setup;

  //--------------------------------------------------
  // edm::InputTag's
  //--------------------------------------------------

  const edm::InputTag m_cscCorrLCTDigiTag;
  const edm::InputTag m_genParticlesTag;
  const edm::InputTag m_simHitsTag;
  const edm::InputTag m_simTracksTag;

  //--------------------------------------------------
  // Conditions information
  //--------------------------------------------------

  edm::ESHandle<CSCGeometry> m_geometry;

  //--------------------------------------------------
  // GEN muon info
  //--------------------------------------------------
  
  float m_gen_muon_pt;

  //---------------------------------------------
  // ROOT tree objects
  //---------------------------------------------
  
  FillCSCTFMuonTree m_fillTree;
  CSCTFMuonTree     m_tree;
  
  //--------------------------------------------------
  // Which muons to analyze? (User-set)
  //--------------------------------------------------
  
  const bool m_analyzeGenMuons;
  const bool m_analyzeDataMuons;

  //--------------------------------------------------
  // Verbosity bool
  //--------------------------------------------------

  const bool m_verbose;

  //--------------------------------------------------
  // Binning info
  //--------------------------------------------------
  
  std::vector<double> m_ptBins;
  std::vector<double> m_etaBins;

  //--------------------------------------------------
  // File names
  //--------------------------------------------------

  std::string m_fileName;
  std::string m_histName;

  //--------------------------------------------------
  // Plot Storage
  //--------------------------------------------------
  
  MuonBinnedPlotStorage* m_plotStorage;

  //--------------------------------------------------
  // Various LUTs
  //--------------------------------------------------
  
  CSCSectorReceiverLUT * m_srLUTs[6][2];

};

#endif
