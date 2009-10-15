#ifndef ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONPTANALYZER_H
#define ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONPTANALYZER_H

// Framework
#include "FWCore/Framework/interface/ESHandle.h"

// Sector receivers and LUTs
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "L1Trigger/CSCTrackFinder/src/CSCTFDTReceiver.h"

// Data formats
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/CSCTriggerContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/TrackStub.h"

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
  // trig prim information storage
  //--------------------------------------------------  

  struct CSCTPInfo {
    int station;
    int ring;
    int sector;
    int frBit;    
    int gblPhi;
    int gblEta;
    float eta;
  };
  
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
  void analyzeMuonTrigPrims();

  //--------------------------------------------------
  // Conditions getter function
  //--------------------------------------------------

  void getConditions();

  //--------------------------------------------------
  // Build objects in the constructor
  //--------------------------------------------------
  
  void buildSRLUTs();
    
  //--------------------------------------------------
  // Fill a list of stubs with DT and CSC information
  //--------------------------------------------------

  void fillStubList ( std::vector<csctf::TrackStub> & v_stubs );

  //--------------------------------------------------
  // Get the global eta and phi coordinates of a stub
  //--------------------------------------------------
  
  void getStubCoordinates( csctf::TrackStub & stub, 
			   int& gblEta_eta, int& gblPhi_phi, 
			   float& cms_eta, float& cms_phi );

  //--------------------------------------------------
  // Analyze combinations of stubs
  //--------------------------------------------------
  
  void analyzeStubCombos    ( const std::vector < std::vector < std::vector< CSCTPInfo > > > & v_tpInfo );
  void analyzeStubCombos_new( const std::vector < std::vector < std::vector< CSCTPInfo > > > & v_tpInfo );

  //--------------------------------------------------
  // Get the quality of a track
  //--------------------------------------------------
  
  void getTrackQuality( const std::vector < std::vector < std::vector < CSCTPInfo > > > & v_tpInfo,
			std::vector <int>&  trackQuality_bySector);
  
  //--------------------------------------------------
  // Binning functions
  //--------------------------------------------------

  int getPtBin ( float pt  );
  int getEtaBin( float eta );

  //--------------------------------------------------
  // Helper functions
  //--------------------------------------------------
  
  int reducedComboHitId ( const int stationId1, const int stationId2);
  int reducedComboFRBId ( const int frbit1    , const int frbit2    );
  void int2bin (uint32_t val, char* string);

  //--------------------------------------------------
  // Event and setup pointers
  //--------------------------------------------------

  const edm::Event* m_event;
  const edm::EventSetup* m_setup;

  //--------------------------------------------------
  // edm::InputTag's
  //--------------------------------------------------

  const edm::InputTag m_cscCorrLCTDigiTag;
  const edm::InputTag m_dttpDigiTag;
  const edm::InputTag m_genParticlesTag;

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
  std::string m_dphiHistName;
  std::string m_detaHistName;
  std::string m_etaHistName;

  //--------------------------------------------------
  // Limits of stations and sectors
  //--------------------------------------------------

  const int m_firstStation, m_lastStation;
  const int m_firstSector , m_lastSector ;
  
  const int m_maxHitCombo ;
  const int m_maxFrbCombo ;
  const int m_nQualityBins;
  
  //--------------------------------------------------
  // Plot Storage
  //--------------------------------------------------
  
  MuonBinnedPlotStorage* m_plotStorage;

  //--------------------------------------------------
  // Various LUTs
  //--------------------------------------------------
  
  CSCSectorReceiverLUT * m_srLUTs[5][6][2];
  CSCTFDTReceiver* m_dtrc;

};

#endif
