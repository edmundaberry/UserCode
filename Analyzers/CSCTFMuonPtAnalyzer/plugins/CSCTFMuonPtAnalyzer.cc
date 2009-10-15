#include <memory>
#include <stdio.h>
#include <iostream>

// ROOT
#include "TMath.h"

// Framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Data formats
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCStatusDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// LUTs
#include "L1Trigger/CSCCommonTrigger/interface/CSCFrontRearLUT.h"

// Geometry/Conditions info
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

// My headers
#include "Analyzers/CSCTFMuonPtAnalyzer/interface/CSCTFMuonPtAnalyzer.h"

CSCTFMuonPtAnalyzer::CSCTFMuonPtAnalyzer(const edm::ParameterSet& iConfig):
  m_cscCorrLCTDigiTag ( iConfig.getParameter<edm::InputTag>("CSCCorrLCTDigiTag")),
  m_dttpDigiTag       ( iConfig.getParameter<edm::InputTag>("DTTPDigiTag")),
  m_genParticlesTag   ( iConfig.getParameter<edm::InputTag>("GenParticlesTag")),
  m_analyzeGenMuons   ( iConfig.getParameter<bool>         ("AnalyzeGenMuons")),
  m_analyzeDataMuons  ( iConfig.getParameter<bool>         ("AnalyzeDataMuons")),
  m_verbose           ( iConfig.getParameter<bool>         ("Verbose")),
  m_ptBins            ( iConfig.getParameter < std::vector <double> > ("ScalePt")),
  m_etaBins           ( iConfig.getParameter < std::vector <double> > ("ScaleEta")),
  m_fileName          ( iConfig.getParameter<std::string>  ("FileName")),
  m_dphiHistName      ( iConfig.getParameter<std::string>  ("DPhiHistName")),
  m_detaHistName      ( iConfig.getParameter<std::string>  ("DEtaHistName")),
  m_etaHistName       ( iConfig.getParameter<std::string>  ("EtaHistName")),
  m_firstStation      ( 0 ), // includes the DT's
  m_lastStation       ( 4 ),
  m_firstSector       ( 1 ),
  m_lastSector        ( 6 ), 
  m_maxHitCombo       ( 10 ),
  m_maxFrbCombo       ( 4 ),
  m_nQualityBins      ( 4 )
{
  
  buildSRLUTs ();

  m_dtrc = new CSCTFDTReceiver();
  
  m_fillTree.init(m_fileName,&m_tree);

  m_plotStorage = new MuonBinnedPlotStorage (m_maxHitCombo, m_maxFrbCombo, m_ptBins.size() , 
					     m_etaBins.size() , m_nQualityBins, 
					     m_dphiHistName, m_detaHistName, m_etaHistName );

}

CSCTFMuonPtAnalyzer::~CSCTFMuonPtAnalyzer() {
  if (m_plotStorage) delete m_plotStorage;
  if (m_dtrc) delete m_dtrc;
}

void CSCTFMuonPtAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //--------------------------------------------------
  // Assign event and setup pointers
  //--------------------------------------------------

  m_event = &iEvent;
  m_setup = &iSetup;

  //--------------------------------------------------
  // Initialize the ROOT tree
  //--------------------------------------------------  
  
  m_tree.init();

  //--------------------------------------------------
  // Get run and event
  //--------------------------------------------------

  int run   = (int) m_event -> id().run();
  int event = (int) m_event -> id().event();

  m_tree.run   = run;
  m_tree.event = event;
  
  //--------------------------------------------------
  // Analyze "real" muons, whether from data or MC
  //--------------------------------------------------

  if (m_analyzeGenMuons ) analyzeGenMuons ();
  if (m_analyzeDataMuons) analyzeDataMuons();

  //--------------------------------------------------
  // Analyze the trigger primitives
  //--------------------------------------------------

  analyzeMuonTrigPrims();

  //--------------------------------------------------
  // Fill ROOT tree
  //--------------------------------------------------

  m_fillTree.fill();

}

//--------------------------------------------------
// Analyze generator muons from MC, if asked to
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::analyzeGenMuons(){

  edm::Handle< reco::GenParticleCollection > GenParticles;
  
  bool gotGenParticles = m_event -> getByLabel (m_genParticlesTag, GenParticles);  
  
  if (!gotGenParticles){
    edm::LogWarning("CSCTFMuonPtAnalyzer") << "Could not extract: " << m_genParticlesTag;
    return;
  }

  reco::GenParticleCollection::const_iterator genParticle     = GenParticles -> begin();
  reco::GenParticleCollection::const_iterator genParticle_end = GenParticles -> end();
  
  int ngmu = 0;

  for (; genParticle != genParticle_end; ++genParticle){

    if ( ngmu >= CSCTFMuonTree::MAXNGMU ) continue;

    int pdg = (int) genParticle -> pdgId();
    
    if ( abs(pdg) != 13 ) continue;
    
    reco::Particle::LorentzVector p4 = genParticle -> p4();
    
    float gen_pt  = (float) p4.pt();
    float gen_phi = (float) p4.phi();
    float gen_eta = (float) p4.eta();

    m_tree.gmu_pt [ngmu] = gen_pt;
    m_tree.gmu_phi[ngmu] = gen_phi;
    m_tree.gmu_eta[ngmu] = gen_eta;

    if (m_verbose) 
      std::cout << " GEN pt = " << gen_pt  << ", " 
		<< "eta = "     << gen_eta << ", " 
		<< "phi = "     << gen_phi << std::endl;
    
    m_gen_muon_pt = gen_pt;

    ++ngmu;
    
  }

  int ptBin = getPtBin( m_gen_muon_pt );

  m_tree.ptBin = ptBin;
  m_tree.ngmu  = ngmu;

}

//--------------------------------------------------
// Analyze muons from data, if asked to
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::analyzeDataMuons(){

  edm::LogWarning("CSCTFMuonPtAnalyzer") << "Analyzing muons from data is not implemented!";
  return;

}

//--------------------------------------------------
// Once you have your "real" muons, analyze CSC TP's
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::analyzeMuonTrigPrims(){

  //--------------------------------------------------  
  // Combine CSC and DT info into a vector of stubs
  //--------------------------------------------------  
  
  std::vector<csctf::TrackStub> v_stubs;
  fillStubList ( v_stubs );

  //--------------------------------------------------  
  // Now loop over the stubs
  //--------------------------------------------------  

  std::vector<csctf::TrackStub>::iterator stub     = v_stubs.begin();
  std::vector<csctf::TrackStub>::iterator stub_end = v_stubs.end();
  
  std::vector< std::vector< std::vector <CSCTPInfo> > > 
    v_comboInfo ( m_lastSector + 1, std::vector< std::vector< CSCTPInfo> >(m_lastStation + 1, std::vector<CSCTPInfo> () ) );

  int nstub = 0;

  for (; stub != stub_end; ++stub){

    //--------------------------------------------------
    // These values exist for CSC and DT stubs
    //--------------------------------------------------

    int  station  = stub -> station();
    int  sector   = stub -> sector ();
    int  endcap   = stub -> endcap();
    int  cscid    = stub -> cscid();

    //--------------------------------------------------
    // These values won't be filled for the DT stubs
    //--------------------------------------------------

    // Digi values
    int  strip    = -2;
    int  pattern  = -2;
    int  quality  = -2;
    int  bend     = -2;
    int  keyWG    = -2;

    // DetId values
    int  ring     = -2;
    int  chamber  = -2;
    int  frBit    = 0;
    bool bad_phi  = false;

    //--------------------------------------------------
    // Is this in the CSC?
    //--------------------------------------------------

    if ( station != 5 ){

      //--------------------------------------------------
      // Get the CSCDetId and the digi
      //--------------------------------------------------
      
      const CSCDetId * cscDetId = new CSCDetId ( stub -> getDetId() );
      const CSCCorrelatedLCTDigi * digi = stub -> getDigi();

      //--------------------------------------------------
      // Get the digi values
      //--------------------------------------------------
            
      strip   = digi -> getStrip();
      pattern = digi -> getPattern();
      quality = digi -> getQuality();
      bend    = digi -> getBend();
      keyWG   = digi -> getKeyWG();

      //--------------------------------------------------
      // Get the detector values
      //--------------------------------------------------

      ring    = cscDetId -> ring();
      chamber = cscDetId -> chamber();      
      frBit   = CSCFrontRearLUT::getFRBit ( cscDetId -> triggerSector(),
					    CSCTriggerNumbering::triggerSubSectorFromLabels(*cscDetId),
					    cscDetId -> station(),
					    cscDetId -> triggerCscId() );

      //--------------------------------------------------
      // Is this a bad phi chamber?
      //--------------------------------------------------
      
      bad_phi = (chamber > 8 && station == 4 && sector == 2);

      //--------------------------------------------------
      // Memory clean-up, if necessary 
      //--------------------------------------------------
      
      if (cscDetId) delete cscDetId;

    }

    //--------------------------------------------------
    // Get stub coordinates
    //--------------------------------------------------
    
    int stub_gblEta, stub_gblPhi;
    float cms_eta, cms_phi;

    getStubCoordinates( *stub, stub_gblEta,  stub_gblPhi, cms_eta,  cms_phi );

    //--------------------------------------------------
    // Store the digi values
    //--------------------------------------------------

    m_tree.stub_digi_eta    [nstub] = cms_eta;
    m_tree.stub_digi_phi    [nstub] = cms_phi;
    m_tree.stub_digi_strip  [nstub] = strip;
    m_tree.stub_digi_keyWG  [nstub] = keyWG;
    m_tree.stub_digi_quality[nstub] = quality;
    m_tree.stub_digi_badphi [nstub] = (bad_phi ? 1 : 0);
    m_tree.stub_digi_pattern[nstub] = pattern;
    m_tree.stub_digi_bend   [nstub] = bend;    
    m_tree.stub_digi_gblPhi [nstub] = stub_gblEta;
    m_tree.stub_digi_gblEta [nstub] = stub_gblPhi;       			       
    
    //--------------------------------------------------
    // Store the detector values
    //--------------------------------------------------

    m_tree.stub_frBit[nstub] = frBit  ;
    m_tree.stub_stat [nstub] = station;
    m_tree.stub_ring [nstub] = ring   ;
    m_tree.stub_cham [nstub] = chamber;
    m_tree.stub_endc [nstub] = endcap ;
    m_tree.stub_sect [nstub] = sector ;
    m_tree.stub_cscid[nstub] = cscid  ;
        
    //--------------------------------------------------
    // Is this phi OK?
    //--------------------------------------------------

    if (!bad_phi){

      //--------------------------------------------------
      // If it is OK, print out info on the stub
      //--------------------------------------------------

      if (m_verbose) 
	std::cout << "   " << endcap << ":" << station << ":" << ring << ":" << chamber << ", "
		  << "sector = " << sector  << ", " 
		  << "eta = "    << cms_eta << ", "
		  << "phi = "    << cms_phi << ", "
		  << "gblPhi = " << stub_gblPhi << std::endl;

      //--------------------------------------------------
      // Also pass on info about the stub to be made into
      // a combo for momentum assignment
      //--------------------------------------------------

      int comboStation = station;
      if (station == 5) comboStation = 0;

      CSCTPInfo tpInfo;
      
      tpInfo.station = comboStation;
      tpInfo.sector  = sector;
      tpInfo.ring    = ring;
      tpInfo.frBit   = frBit;
      tpInfo.gblPhi  = stub_gblPhi;
      tpInfo.gblEta  = stub_gblEta;
      tpInfo.eta     = cms_eta;      

      v_comboInfo[sector][comboStation].push_back(tpInfo);

    }

    //--------------------------------------------------
    // Say how many stubs we found
    //--------------------------------------------------

    ++nstub;

  }

  m_tree.nstub = nstub;

  analyzeStubCombos ( v_comboInfo );
  
}

//--------------------------------------------------
// Method for getting the track quality across 
// all sectors in the CSC's
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::getTrackQuality( const std::vector < std::vector < std::vector< CSCTPInfo > > > & v_tpInfo,
					   std::vector <int>& v_trackQuality_bySector) {

  //--------------------------------------------------
  // Resize the quality vector to fit all sectors
  //--------------------------------------------------

  v_trackQuality_bySector.resize( m_lastSector + 1 );

  //--------------------------------------------------
  // Loop over the sectors available
  //--------------------------------------------------

  for (int iSector = m_firstSector; iSector <= m_lastSector; ++iSector){
    
    //--------------------------------------------------
    // Allot values in memory for quality assigment
    //--------------------------------------------------

    int  quality;
    int  n_stations_with_hits_in_this_sector = 0;
    bool hit_station_1 = false;
    bool hit_station_2 = false; 
    bool hit_station_3 = false;

    //--------------------------------------------------
    // Get the hits in this sector as ordered by station
    //--------------------------------------------------

    const std::vector < std::vector <CSCTPInfo> > * v_hits_in_this_sector = & v_tpInfo[iSector];

    //--------------------------------------------------
    // Loop over the stations
    //--------------------------------------------------

    for (int iStation = m_firstStation; iStation <= m_lastStation; ++iStation){
      
      //--------------------------------------------------
      // Do not count the DT's in the quality assigment
      //--------------------------------------------------
      
      if ( iStation == 0 ) continue;

      //--------------------------------------------------
      // Get the hits in this sector and station
      //--------------------------------------------------

      const std::vector < CSCTPInfo > * v_hits_in_this_sector_and_station = & v_hits_in_this_sector -> at( iStation );
      
      //--------------------------------------------------
      // Get values you need to assign quality
      //--------------------------------------------------
      
      int n_hits_in_this_sector_and_station = (int) v_hits_in_this_sector_and_station -> size() ;
      
      if ( n_hits_in_this_sector_and_station > 0 ) {

	n_stations_with_hits_in_this_sector++;	
	
	if      ( iStation == 1 ) hit_station_1 = true;
	else if ( iStation == 2 ) hit_station_2 = true;
	else if ( iStation == 3 ) hit_station_3 = true;      
	
      }

    }

    //--------------------------------------------------
    // Assign quality
    //--------------------------------------------------
     
    if      (  n_stations_with_hits_in_this_sector < 2                   ) quality = -2; // can't make a combo
    else if ( !hit_station_2 && !hit_station_3                           ) quality = -1; // no key stations
    else if ( !hit_station_1 && n_stations_with_hits_in_this_sector == 2 ) quality =  1;
    else if (  hit_station_1 && n_stations_with_hits_in_this_sector == 2 ) quality =  2;
    else if (  hit_station_1 && n_stations_with_hits_in_this_sector >= 3 ) quality =  3;
    else                                                                   quality =  0;

    v_trackQuality_bySector [ iSector ] = quality;
    
  }

}

//--------------------------------------------------
// Analyze combinations of stubs 
// and make tracks from them
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::analyzeStubCombos( const std::vector < std::vector < std::vector <CSCTPInfo> > > & v_tpInfo ) {

  //--------------------------------------------------
  // Loop over possible first stations
  //--------------------------------------------------
  
  std::vector<int> v_trackQuality;
  getTrackQuality ( v_tpInfo, v_trackQuality );

  if (m_verbose){
    std::cout << "   " << "I am predicting tracks in sectors: " << std::endl;
    for (int iSector = m_firstSector; iSector <= m_lastSector ; ++iSector){
      int quality = v_trackQuality[iSector];
      if (quality > 0)	std::cout << "          " << iSector << " with quality = " << quality << std::endl;
    }
  }

  std::vector<int> sect_ncombo;
    
  int nsector = 0;

  //--------------------------------------------------    
  // Loop over the sectors
  //--------------------------------------------------    

  for ( int iSector = m_firstSector; iSector <= m_lastSector; ++iSector ){

    //--------------------------------------------------
    // Don't bother with sectors that are of poor quality
    //--------------------------------------------------    

    int quality = v_trackQuality[iSector];    
    if (quality < 0) continue;

    //--------------------------------------------------    
    // Get all hits in this sector
    //--------------------------------------------------    

    const std::vector < std::vector < CSCTPInfo > > * v_hits_in_this_sector = & v_tpInfo[iSector];
       
    //--------------------------------------------------    
    // Initialize values for this sector
    //--------------------------------------------------    

    int  first_combo_station1 = -999;
    int  first_combo_station2 = -999;
    bool noMoreStationsWithHits = false; 
    int  station1 = m_firstStation;
    int  ncombo = 0;

    //--------------------------------------------------    
    // Loop over stations
    //--------------------------------------------------    

    while ( station1 <= ( m_lastStation - 1) && !noMoreStationsWithHits ){
      
      //--------------------------------------------------
      // How many hits in this station? 
      // If none, go to the next station.
      //--------------------------------------------------
            
      const std::vector<CSCTPInfo> * v_station1_hits = & v_hits_in_this_sector -> at ( station1 );
    
      int n_station1_hits = (int) v_station1_hits -> size();
      
      if ( n_station1_hits == 0) {      
	station1++;
	continue;
      }
    
      if (m_verbose)
	std::cout << "   Trying 1st Station is Station " << station1 
		  << " with " << n_station1_hits << " hits..." << std::endl;
      
      //--------------------------------------------------
      // Loop over candidates for the second station, 
      // until you find one that has a hit.
      //
      // Be careful not to go past the last station.
      // 
      // If you never found another non-empty station,
      // leave the loop over station1 candidates
      //--------------------------------------------------
    
      int station2 = station1 + 1;
      int n_station2_hits = 0;
      bool foundNextNonEmptyStation = false;
      
      while ( !foundNextNonEmptyStation && station2 <= m_lastStation){
      
	n_station2_hits = (int) v_hits_in_this_sector -> at ( station2 ).size();
      
	if ( n_station2_hits > 0 ) foundNextNonEmptyStation = true;
	else station2++;
      
      }
    
      if (!foundNextNonEmptyStation) {
	if (m_verbose) std::cout << "          2nd Station cannot be found!" << std::endl;
	noMoreStationsWithHits = true;
	continue;
      }

      if (m_verbose) 
	std::cout << "          " << "2nd Station is Station " << station2 
		  << " with " << n_station2_hits << " hits " << std::endl;
      
      //--------------------------------------------------
      // Get your vector of station 2 hits
      //--------------------------------------------------
      
      const std::vector<CSCTPInfo> * v_station2_hits = & v_hits_in_this_sector -> at ( station2 );
      
      //--------------------------------------------------
      // Now that you have both stations, loop over first
      // station hits
      //--------------------------------------------------
      
      for ( int iHit1 = 0; iHit1 < n_station1_hits; ++iHit1 ) {
      
	const CSCTPInfo * station1_hit = & v_station1_hits -> at ( iHit1 );
	
	//--------------------------------------------------
	// ... and second station hits
	//--------------------------------------------------
      
	for (int iHit2 = 0; iHit2 < n_station2_hits; ++iHit2 ) {
	  
	  const CSCTPInfo * station2_hit = & v_station2_hits -> at ( iHit2 );
	  
	  //--------------------------------------------------
	  // Get the dphi and deta values
	  //--------------------------------------------------
	  
	  int dphi = (station1_hit -> gblPhi) - (station2_hit -> gblPhi);
	  int deta = (station1_hit -> gblEta) - (station2_hit -> gblEta);
	
	  //--------------------------------------------------
	  // Don't bother with hits in different sectors
	  //--------------------------------------------------
	  
	  if ( station1_hit -> sector != station2_hit -> sector ) {
	    if (m_verbose) std::cout << "          " << "- Combo FAIL: Same sector cut " << std::endl;
	    continue;
	  }
	  
	  //--------------------------------------------------
	  // Don't bother with combos outside of our dphi cut
	  //
	  // Is this the first combination of stations?
	  // If it isn't, our dphi cut gets stricter
	  //--------------------------------------------------

	  int isFirstCombo = 0;
	  
	  if ( first_combo_station1 == -999 || 
	       first_combo_station2 == -999  ){
	    first_combo_station1 = station1;
	    first_combo_station2 = station2;
	  }
	  
	  if ( station1 == first_combo_station1 && 
	       station2 == first_combo_station2 ){	  
	    
	    isFirstCombo = 1;	 
	    
	    if ( abs(dphi) > 512 ) {
	      if (m_verbose) std::cout <<"          " << "- Combo FAIL: 1st dphi cut " << std::endl;
	      continue;
	    }
	  }
	  
	  else {
	    if ( abs(dphi) > 256 ) {
	      if (m_verbose) std::cout << "          " << "- Combo FAIL: 2nd dphi cut " << std::endl;
	      continue;
	    }
	  }

	  //--------------------------------------------------
	  // Get information on the Station 1 (of 4) ring
	  //--------------------------------------------------
	  
	  int s1ring = 0;
	  int s1ring_bit = 0;
	  
	  if (station1 == 1) s1ring = station1_hit -> ring;
	  if (station2 == 1) s1ring = station2_hit -> ring;	  

	  if (s1ring == 1) s1ring_bit = 1;

	  //--------------------------------------------------
	  // Get binning information
	  //--------------------------------------------------
	  
	  int combo_hitId  = reducedComboHitId ( station1 , station2 );
	  int combo_frbId  = reducedComboFRBId ( station1_hit -> frBit , station2_hit -> frBit );	  
	  int etaBin       = getEtaBin         ( station2_hit -> eta );
	  
	  //--------------------------------------------------
	  // Store info about this combo in the tree
	  //--------------------------------------------------
	  
	  m_tree.sect_combo_s1ring [nsector][ncombo] = s1ring_bit;
	  m_tree.sect_combo_isFirst[nsector][ncombo] = isFirstCombo;
	  m_tree.sect_combo_etaBin [nsector][ncombo] = etaBin;
	  m_tree.sect_combo_dphi   [nsector][ncombo] = dphi;
	  m_tree.sect_combo_deta   [nsector][ncombo] = deta;
	  m_tree.sect_combo_hitId  [nsector][ncombo] = combo_hitId;
	  m_tree.sect_combo_frbId  [nsector][ncombo] = combo_frbId;	  
	  m_tree.sect_combo_stat1  [nsector][ncombo] = station1;
	  m_tree.sect_combo_stat2  [nsector][ncombo] = station2;
	  m_tree.sect_combo_eta1   [nsector][ncombo] = station1_hit -> gblEta;
	  m_tree.sect_combo_eta2   [nsector][ncombo] = station2_hit -> gblEta;
	  m_tree.sect_combo_phi1   [nsector][ncombo] = station1_hit -> gblPhi;	  
	  m_tree.sect_combo_phi2   [nsector][ncombo] = station2_hit -> gblPhi;
	  
	  if (m_verbose) std::cout << "          " << "+ Combo #" << ncombo + 1 << " in Sector #" << nsector + 1 << " PASS: dphi = " << dphi << ", isFirst = " << isFirstCombo <<  std::endl;      
	  
	  //--------------------------------------------------
	  // Fill the histograms
	  //--------------------------------------------------
	  
	  m_plotStorage -> fillEtaBin ( combo_hitId, quality, m_tree.ptBin, isFirstCombo, etaBin );
	  m_plotStorage -> fillDPhi   ( combo_hitId, combo_frbId, m_tree.ptBin, etaBin, s1ring_bit, quality, isFirstCombo, dphi );
	  m_plotStorage -> fillDEta   ( combo_hitId, combo_frbId, m_tree.ptBin, etaBin, s1ring_bit, quality, isFirstCombo, deta );

	  //--------------------------------------------------
	  // Count the number of combos
	  //--------------------------------------------------
	  
	  ncombo++;
	  
	}
      }
      
      //--------------------------------------------------
      // Get ready to look at the next set of combinations
      //--------------------------------------------------
      
      station1 = station2;
      
    }

    //--------------------------------------------------
    // Store the total sector information
    //--------------------------------------------------

    m_tree.sect_quality[nsector] = quality;
    m_tree.sect_num    [nsector] = iSector;
    m_tree.sect_ncombo [nsector] = ncombo;
    
    //--------------------------------------------------
    // Count sectors
    //--------------------------------------------------

    nsector++;
    
  }

  //--------------------------------------------------
  // Store the number of sectors
  //--------------------------------------------------

  m_tree.nsector = nsector;
  
}


int CSCTFMuonPtAnalyzer::reducedComboFRBId ( const int frbit1, const int frbit2 ){

  int id;
  
  id  = frbit2 << 1;
  id += frbit1;

  return id;

}

int CSCTFMuonPtAnalyzer::reducedComboHitId ( const int stationId1, const int stationId2){
  
  assert( stationId1 != stationId2 );

  int binary_id;
  binary_id  = (1<<stationId1);
  binary_id += (1<<stationId2);

  if      (binary_id == 3 ) return 0;
  else if (binary_id == 5 ) return 1; // * 
  else if (binary_id == 6 ) return 2;
  else if (binary_id == 9 ) return 3;
  else if (binary_id == 10) return 4;
  else if (binary_id == 12) return 5; // *
  else if (binary_id == 17) return 6;
  else if (binary_id == 18) return 7;
  else if (binary_id == 20) return 8;
  else if (binary_id == 24) return 9;
  else {
    edm::LogWarning("CSCTFMuonPtAnalyzer") << "Unbinnable combo id: " << binary_id;
    return -999;
  }

  return binary_id;
}

//--------------------------------------------------
// Get stub loocation
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::getStubCoordinates( csctf::TrackStub & stub, 
					      int& gblEta_eta, int& gblPhi_phi, 
					      float& cms_eta, float& cms_phi ){
    
  //--------------------------------------------------
  // Get stub detector information
  //--------------------------------------------------

  unsigned int endcap    = stub.endcap();
  unsigned int station   = stub.station();
  unsigned int sector    = stub.sector();
  unsigned int subsector = stub.subsector();
  unsigned int cscid     = stub.cscid();
  
  //--------------------------------------------------
  // Stations 1-4 (CSC stations) get their coordinates
  // from the SR LUTs
  //--------------------------------------------------

  if( station != 5) {

    //--------------------------------------------------
    // Get global and local phi/eta structs
    //--------------------------------------------------

    unsigned int fpga = ( subsector ? subsector-1 : station );
    lclphidat lclPhi = m_srLUTs[fpga][sector-1][endcap-1] -> localPhi( stub.getStrip(), 
								       stub.getPattern(), 
								       stub.getQuality(), 
								       stub.getBend());
    
    gblphidat gblPhi = m_srLUTs[fpga][sector-1][endcap-1] -> globalPhiME( lclPhi.phi_local, 
									  stub.getKeyWG(), cscid);
    
    gbletadat gblEta = m_srLUTs[fpga][sector-1][endcap-1] -> globalEtaME(lclPhi.phi_bend_local, 
									 lclPhi.phi_local, 
									 stub.getKeyWG(), cscid);
    
    //--------------------------------------------------
    // Set the packed values for the stub
    //--------------------------------------------------

    stub.setEtaPacked(gblEta.global_eta);
    stub.setPhiPacked(gblPhi.global_phi);
    
  }  
  
  //--------------------------------------------------
  // Get the global eta and phi values
  //--------------------------------------------------

  gblEta_eta = stub.etaPacked();
  gblPhi_phi = stub.phiPacked();

  //--------------------------------------------------
  // Get the "real" eta value 
  //--------------------------------------------------
  
  cms_eta = gblEta_eta/127. * 1.5 + 0.9;

  //--------------------------------------------------
  // Get the "real" phi value in radians
  //--------------------------------------------------
  
  float cmsPhiInDegrees = (gblPhi_phi/4096.*62.) - 1. + 15.+60.*sector;  
  if(cmsPhiInDegrees > 360) cmsPhiInDegrees -= 360;
  
  cms_phi = cmsPhiInDegrees * TMath::Pi() / 180.0;
  if (cms_phi > TMath::Pi()) cms_phi -= (2.0 * TMath::Pi());
    
}

//--------------------------------------------------
// Process for building LUTs within the constructor
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::buildSRLUTs(){
  
  bzero(m_srLUTs,sizeof(m_srLUTs));
  bool TMB07=true;
  edm::ParameterSet srLUTset;
  srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
  srLUTset.addUntrackedParameter<bool>("Binary",   false);
  srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");
  for(int endcap = 1; endcap<=2; endcap++) {
    for(int sector=1; sector<=6; sector++) {
      for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
	if(station==1)
	  for(int subSector=0; subSector<2; subSector++)
	    m_srLUTs[fpga++][sector-1][endcap-1] = new CSCSectorReceiverLUT(endcap, sector, subSector+1, station, srLUTset, TMB07);
	else
	  m_srLUTs[fpga++][sector-1][endcap-1] = new CSCSectorReceiverLUT(endcap, sector, 0, station, srLUTset, TMB07);
      }
    }
  }
}

int CSCTFMuonPtAnalyzer::getEtaBin ( float eta ){

  if ( eta > m_etaBins.back() ) return m_etaBins.size();
  if ( eta < m_etaBins[0]     ) return -1;
  
  int size = m_etaBins.size();

  for (int iEtaBin = 0; iEtaBin < size - 1; ++iEtaBin){
    double bin_min_eta = m_etaBins[iEtaBin];
    double bin_max_eta = m_etaBins[iEtaBin+1];
    
    if (eta >= bin_min_eta && eta < bin_max_eta) return iEtaBin;
  }

  edm::LogWarning("CSCTFMuonPtAnalyzer") << "Unbinnable eta value: " << eta;
  return -2;
}

int CSCTFMuonPtAnalyzer::getPtBin( float pt ){

  if ( pt > m_ptBins.back()) return m_ptBins.size();
  if ( pt < m_ptBins[0]    ) return -1;

  int size = m_ptBins.size();

  for (int iPtBin = 0; iPtBin < size - 1; ++iPtBin){

    double bin_min_pt = m_ptBins[iPtBin];
    double bin_max_pt = m_ptBins[iPtBin+1];
    
    if (pt >= bin_min_pt && pt < bin_max_pt) return iPtBin;

  }
  
  edm::LogWarning("CSCTFMuonPtAnalyzer") << "Unbinnable pt value: " << pt;

  return -2;

}

//--------------------------------------------------
// Beginning and ending job functions
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::beginJob(){}

void CSCTFMuonPtAnalyzer::endJob()  {
  m_fillTree.finalize();
  m_plotStorage -> save();
}

//--------------------------------------------------
// int2bin function, useful for checking bit logic
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::int2bin (uint32_t val, char *string){
  int i, j;
  unsigned long tmpval;
  char *stringPtr = &string[0];
  int cnt;
  unsigned long tval = val;

  for (cnt=0; tval>0; tval=tval/2,cnt++) ;

  for (j=0, i = cnt-1; val>0; i-- ,j++){
    
    tmpval = val%2;
    val = val/2;

    if (tmpval) *(stringPtr + i) = '0' + 1;
    else        *(stringPtr + i) = '0' + 0;

  }
  
  *(stringPtr + cnt) = '\0';
  return;
}

//--------------------------------------------------
// Fill a list of stubs from the CSC's and the DT's
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::fillStubList ( std::vector<csctf::TrackStub> & v_stubs ) {

  //--------------------------------------------------
  // Get the handles
  //--------------------------------------------------
  
  edm::Handle< L1MuDTChambPhContainer > DTTPs;
  edm::Handle< CSCCorrelatedLCTDigiCollection > CorrLCTs;
  bool gotCSCCorrLCTDigis = m_event -> getByLabel (m_cscCorrLCTDigiTag, CorrLCTs);
  bool gotDTTPDigis       = m_event -> getByLabel (m_dttpDigiTag,       DTTPs   );
  
  if (!gotDTTPDigis){
    edm::LogWarning("CSCTFMuonPtAnalyzer") << "Could not extract: " << m_dttpDigiTag;
    return;
  }

  if (!gotCSCCorrLCTDigis){
    edm::LogWarning("CSCTFMuonPtAnalyzer") << "Could not extract: " << m_cscCorrLCTDigiTag;
    return;
  }  

  //--------------------------------------------------
  // Make a special container of stubs to fill
  //--------------------------------------------------

  CSCTriggerContainer<csctf::TrackStub> stub_list;

  //--------------------------------------------------
  // Process the DT trig prims to get the DT stubs
  //--------------------------------------------------
  
  CSCTriggerContainer<csctf::TrackStub> dtStubs = m_dtrc -> process(DTTPs.product());
  stub_list.push_many(dtStubs);
  
  //--------------------------------------------------
  // Loop over the CSC trig prims to get the CSC stubs
  //--------------------------------------------------
  
  CSCCorrelatedLCTDigiCollection::DigiRangeIterator chamber     = CorrLCTs.product() -> begin();
  CSCCorrelatedLCTDigiCollection::DigiRangeIterator chamber_end = CorrLCTs.product() -> end();
  
  for (; chamber != chamber_end; ++chamber) {
    
    CSCDetId * cscDetId = &((*chamber).first);

    CSCCorrelatedLCTDigiCollection::Range corrLCTDigi_range = CorrLCTs.product() -> get( *cscDetId  );    
    CSCCorrelatedLCTDigiCollection::const_iterator corrLCTDigi     = corrLCTDigi_range.first;
    CSCCorrelatedLCTDigiCollection::const_iterator corrLCTDigi_end = corrLCTDigi_range.second;
    CSCCorrelatedLCTDigiCollection::const_iterator bestCorrLCTDigi;
    
    //--------------------------------------------------
    // Find the best-quality digi in this chamber
    //--------------------------------------------------

    int max_quality = -999;

    for (; corrLCTDigi != corrLCTDigi_end; ++corrLCTDigi ){    

      int this_digi_quality = (int) corrLCTDigi -> getQuality();

      if (this_digi_quality > max_quality) {
	max_quality = this_digi_quality;
	bestCorrLCTDigi = corrLCTDigi;
      }
    }

    
    csctf::TrackStub theStub(*bestCorrLCTDigi, *cscDetId );
    stub_list.push_back(theStub);

  }

  //--------------------------------------------------
  // Pass the filled list back to be analyzed
  //--------------------------------------------------

  v_stubs = stub_list.get();

}


//--------------------------------------------------
// Define framework module
//--------------------------------------------------

DEFINE_FWK_MODULE(CSCTFMuonPtAnalyzer);
