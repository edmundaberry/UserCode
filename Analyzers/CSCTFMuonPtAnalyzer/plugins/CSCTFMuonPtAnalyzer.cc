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
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
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
  m_genParticlesTag   ( iConfig.getParameter<edm::InputTag>("GenParticlesTag")),
  m_simHitsTag        ( iConfig.getParameter<edm::InputTag>("SimHitsTag")),
  m_simTracksTag      ( iConfig.getParameter<edm::InputTag>("SimTracksTag")),
  m_analyzeGenMuons   ( iConfig.getParameter<bool>         ("AnalyzeGenMuons")),
  m_analyzeDataMuons  ( iConfig.getParameter<bool>         ("AnalyzeDataMuons")),
  m_verbose           ( iConfig.getParameter<bool>         ("Verbose")),
  m_ptBins            ( iConfig.getParameter < std::vector <double> > ("ScalePt")),
  m_etaBins           ( iConfig.getParameter < std::vector <double> > ("ScaleEta")),
  m_fileName          ( iConfig.getParameter<std::string>  ("FileName")),
  m_histName          ( iConfig.getParameter<std::string>  ("HistName"))
{
  
  buildSRLUTs ();

  m_fillTree.init(m_fileName,&m_tree);

  m_plotStorage = new MuonBinnedPlotStorage (6, 4, m_ptBins.size()  , m_etaBins.size() , m_histName);

}

CSCTFMuonPtAnalyzer::~CSCTFMuonPtAnalyzer() {
  if (m_plotStorage) delete m_plotStorage;
}

void CSCTFMuonPtAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //--------------------------------------------------
  // Assign event and setup pointers
  //--------------------------------------------------

  m_event = &iEvent;
  m_setup = &iSetup;

  // iSetup.get<MuonGeometryRecord>().get( m_geometry );

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
  // Get conditions
  //--------------------------------------------------

  getConditions();

  //--------------------------------------------------
  // Analyze "real" muons, whether from data or MC
  //--------------------------------------------------

  if (m_analyzeGenMuons ) analyzeGenMuons ();
  if (m_analyzeDataMuons) analyzeDataMuons();

  //--------------------------------------------------
  // Analyze simulation muon hits
  //--------------------------------------------------

  // analyzeCSCSimHits();

  //--------------------------------------------------
  // Analyze the trigger primitives
  //--------------------------------------------------

  analyzeCSCTrigPrims();

  //--------------------------------------------------
  // Fill ROOT tree
  //--------------------------------------------------

  m_fillTree.fill();

}

//--------------------------------------------------
// Get conditions information
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::getConditions(){
  
  
  
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
// Analyze SIM hits 
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::analyzeCSCSimHits(){

  //--------------------------------------------------
  // Get the SimHits
  //--------------------------------------------------
  
  edm::Handle< edm::PSimHitContainer > SimHits;
  edm::Handle< edm::SimTrackContainer > SimTracks;

  bool gotSimHits = m_event -> getByLabel (m_simHitsTag, SimHits);
  if (!gotSimHits){
    edm::LogWarning("CSCTFMuonPtAnalyzer") << "Could not extract: " << m_simHitsTag;
    return;
  }

  //--------------------------------------------------
  // Store the last CSCDetId you ran over
  //--------------------------------------------------

  CSCDetId lastDetId;
  int lastPartId =-999;

  //--------------------------------------------------
  // Loop over the sim hits
  //--------------------------------------------------  

  edm::PSimHitContainer::const_iterator simHit     = SimHits.product() -> begin();
  edm::PSimHitContainer::const_iterator simHit_end = SimHits.product() -> end();

  int nsh = 0;

  for (; simHit != simHit_end; ++simHit){

    if ( nsh >= CSCTFMuonTree::MAXNSIMHIT ) continue;
    
    //--------------------------------------------------
    // Sim hits will contain lots of electrons from
    // scattering.  Is this a muon or an electron?
    //--------------------------------------------------

    int pdg = simHit -> particleType();
    
    if (abs(pdg) != 13) continue;

    //--------------------------------------------------
    // Next make sure this is within the CSC
    //--------------------------------------------------

    DetId * detId = new DetId( simHit -> detUnitId() );
    
    if ( detId -> det()      != DetId::Muon ||
    	 detId -> subdetId() != MuonSubdetId::CSC ){
      if (detId) delete detId;
      continue;
    }
    
    //--------------------------------------------------
    // If this is within the CSC, and we haven't looked
    // at this particular particle / chamber before,
    // store the information about the sim hit
    //--------------------------------------------------
    CSCDetId cscDetId = (CSCDetId) simHit -> detUnitId();
    const CSCLayer* cscLayer = m_geometry -> layer ( cscDetId );
    
    LocalPoint  hitLP = simHit -> localPosition();
    GlobalPoint hitGP = cscLayer -> toGlobal(hitLP);
    
    float phi = hitGP.phi();
    float eta = hitGP.eta();
    
    int frBit = CSCFrontRearLUT::getFRBit ( cscDetId.triggerSector(),
					    CSCTriggerNumbering::triggerSubSectorFromLabels(cscDetId),
					    cscDetId.station(),
					    cscDetId.triggerCscId() );
    
    m_tree.sh_pdg  [nsh] = pdg;
    m_tree.sh_eta  [nsh] = eta;
    m_tree.sh_phi  [nsh] = phi;
    m_tree.sh_frBit[nsh] = frBit;
    m_tree.sh_stat [nsh] = cscDetId.station();
    m_tree.sh_ring [nsh] = cscDetId.ring();
    m_tree.sh_cham [nsh] = cscDetId.chamber();
    m_tree.sh_layr [nsh] = cscDetId.layer();
    m_tree.sh_endc [nsh] = cscDetId.endcap();
    m_tree.sh_sect [nsh] = cscDetId.triggerSector();
    m_tree.sh_cscid[nsh] = cscDetId.triggerCscId();
    
    // SimHitInfo simHitInfo = { cscDetId.chamberId(), pdg, eta, phi };
    // m_simHits.push_back(simHitInfo);

    //--------------------------------------------------
    // Update the last particle and chamber information
    //--------------------------------------------------
    
    lastDetId  = cscDetId;
    lastPartId = pdg;

    //--------------------------------------------------
    // Update the number of sim hits
    //--------------------------------------------------

    ++nsh;
  }
  
  m_tree.nsh = nsh;
  
}


//--------------------------------------------------
// Once you have your "real" and sim muons, analyze CSC TP's
//--------------------------------------------------


void CSCTFMuonPtAnalyzer::analyzeCSCTrigPrims(){

  //--------------------------------------------------
  // Get the handle
  //--------------------------------------------------
  
  edm::Handle< CSCCorrelatedLCTDigiCollection > CorrLCTs;

  bool gotCSCCorrLCTDigis = m_event -> getByLabel (m_cscCorrLCTDigiTag, CorrLCTs);
  
  if (!gotCSCCorrLCTDigis){
    edm::LogWarning("CSCTFMuonPtAnalyzer") << "Could not extract: " << m_cscCorrLCTDigiTag;
    return;
  }  

  //--------------------------------------------------  
  // Loop over the chambers
  // Note that the DigiRangeIterator represents chamber
  //    -- First  element is the chamber DetId 
  //    -- Second element is a **collection** of 
  //       associated digis.  May be up to 2.
  //--------------------------------------------------

  CSCCorrelatedLCTDigiCollection::DigiRangeIterator chamber     = CorrLCTs.product() -> begin();
  CSCCorrelatedLCTDigiCollection::DigiRangeIterator chamber_end = CorrLCTs.product() -> end();

  int nl1detid = 0;
  
  std::vector< std::vector<CSCTPInfo> > v_comboInfo ( 4, std::vector<CSCTPInfo>() );

  for (; chamber != chamber_end; ++chamber){   

    //--------------------------------------------------
    // Don't look at more chambers than we can store
    //--------------------------------------------------

    if ( nl1detid >= CSCTFMuonTree::MAXNL1DETID ) continue;

    //--------------------------------------------------
    // Get a pointer to the det id that identifies this
    // chamber
    //--------------------------------------------------

    CSCDetId * cscDetId = &((*chamber).first);

    //--------------------------------------------------
    // Cache information about the CSCDetId
    //--------------------------------------------------

    int station = cscDetId -> station();	  
    int ring    = cscDetId -> ring();	  
    int chamber = cscDetId -> chamber();	  
    int layer   = cscDetId -> layer();	  
    int endcap  = cscDetId -> endcap();	  
    int sector  = cscDetId -> triggerSector(); 
    int cscid   = cscDetId -> triggerCscId();  
    int frBit   = CSCFrontRearLUT::getFRBit ( sector, CSCTriggerNumbering::triggerSubSectorFromLabels(*cscDetId), station, cscid );
    
    //--------------------------------------------------
    // Store information about the CSCDetId
    //--------------------------------------------------

    m_tree.l1detid_frBit[nl1detid] = frBit  ;
    m_tree.l1detid_stat [nl1detid] = station;
    m_tree.l1detid_ring [nl1detid] = ring   ;
    m_tree.l1detid_cham [nl1detid] = chamber;
    m_tree.l1detid_layr [nl1detid] = layer  ;
    m_tree.l1detid_endc [nl1detid] = endcap ;
    m_tree.l1detid_sect [nl1detid] = sector ;
    m_tree.l1detid_cscid[nl1detid] = cscid  ;

    //--------------------------------------------------
    // For each chamber DetId, get the range of LCT digis
    //--------------------------------------------------
    
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

    //--------------------------------------------------
    // Cache information about this best digi
    //--------------------------------------------------

    int strip   = (int) bestCorrLCTDigi -> getStrip();
    int pattern = (int) bestCorrLCTDigi -> getPattern();
    int quality = (int) bestCorrLCTDigi -> getQuality();
    int bend    = (int) bestCorrLCTDigi -> getBend();
    int keyWG   = (int) bestCorrLCTDigi -> getKeyWG();

    //--------------------------------------------------
    // Store information about this best digi
    //--------------------------------------------------

    m_tree.l1detid_digi_strip      [nl1detid] = strip;
    m_tree.l1detid_digi_pattern    [nl1detid] = pattern;
    m_tree.l1detid_digi_quality    [nl1detid] = quality;
    m_tree.l1detid_digi_bend       [nl1detid] = bend;
    m_tree.l1detid_digi_keyWG      [nl1detid] = keyWG;

    //--------------------------------------------------
    // Cache information about the location of the hits
    //--------------------------------------------------
    
    int lclPhi_phi_bend, lclPhi_phi, gblPhi_phi, gblEta_eta;
    float cms_eta, cms_phi;  
    bool bad_phi;
      
    getHitCoordinates (*cscDetId , *bestCorrLCTDigi,
		       lclPhi_phi_bend, lclPhi_phi, gblPhi_phi, gblEta_eta,
		       cms_eta, cms_phi, bad_phi );

    if ( bad_phi )  m_tree.l1detid_digi_badphi[nl1detid] = 1;
    else {
      m_tree.l1detid_digi_badphi[nl1detid] = 0;
      if (m_verbose) 
	std::cout << *cscDetId << ", " 
		  << "sector = " << sector  << ", " 
		  << "eta = "    << cms_eta << ", "
		  << "phi = "    << cms_phi << ", "
		  << "gblPhi = " << gblPhi_phi << std::endl;
    }

    //--------------------------------------------------
    // Store information about the location of the hits
    //--------------------------------------------------

    m_tree.l1detid_digi_eta        [nl1detid] = cms_eta;
    m_tree.l1detid_digi_phi        [nl1detid] = cms_phi;    
    m_tree.l1detid_digi_lclPhi     [nl1detid] = lclPhi_phi;
    m_tree.l1detid_digi_lclPhiBend [nl1detid] = lclPhi_phi_bend;
    m_tree.l1detid_digi_gblPhi     [nl1detid] = gblPhi_phi;
    m_tree.l1detid_digi_gblEta     [nl1detid] = gblEta_eta;    

    //--------------------------------------------------
    // Save important values about the collection of 
    // digis to the tree, if the phi wasn't bad
    //--------------------------------------------------
    
    if (!bad_phi){

      CSCTPInfo tpInfo;
      
      tpInfo.station = cscDetId -> station();
      tpInfo.sector  = cscDetId -> triggerSector();
      tpInfo.ring    = cscDetId -> ring();
      tpInfo.frBit   = frBit;
      tpInfo.gblPhi  = gblPhi_phi;
      tpInfo.eta     = cms_eta;

      // Store trigger primitives to make combos

      v_comboInfo[station - 1].push_back(tpInfo);

    }

    ++nl1detid;

  }

  m_tree.nl1detid = nl1detid;

  analyzeCombos ( v_comboInfo );

}

void CSCTFMuonPtAnalyzer::analyzeCombos( const std::vector < std::vector < CSCTPInfo > > & v_tpInfo ) {

  //--------------------------------------------------
  // Loop over possible first stations
  //--------------------------------------------------
  
  int first_combo_station1 = -999;
  int first_combo_station2 = -999;
  
  int ncombo = 0;
  int station1 = 0;
  
  bool noMoreStationsWithHits = false;
  
  while ( station1 <= 2 && !noMoreStationsWithHits ){
    
    //--------------------------------------------------
    // How many hits in this station? 
    // If none, go to the next station.
    //--------------------------------------------------
    
    const std::vector<CSCTPInfo> * station1_hits = & v_tpInfo [station1];
    
    int nStation1Hits = (int) station1_hits -> size();
    
    if ( nStation1Hits == 0) {      
      station1++;
      continue;
    }
    
    if (m_verbose)
      std::cout << "  Trying 1st Station is Station " << station1 + 1 
		<< " with " << nStation1Hits << " hits..." << std::endl;
    
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
    int nStation2Hits = 0;
    bool foundNextNonEmptyStation = false;
    
    while ( !foundNextNonEmptyStation && station2 <= 3){
      
      nStation2Hits = (int) v_tpInfo[station2].size();
      
      if ( nStation2Hits > 0 ) foundNextNonEmptyStation = true;
      else station2++;
      
    }
    
    if (!foundNextNonEmptyStation) {
      if (m_verbose) std::cout << "         2nd Station cannot be found!" << std::endl;
      noMoreStationsWithHits = true;
      continue;
    }

    if (m_verbose) 
      std::cout << "         " << "2nd Station is Station " << station2 + 1 
		<< " with " << nStation2Hits << " hits " << std::endl;

    //--------------------------------------------------
    // Get your vector of station 2 hits
    //--------------------------------------------------
    
    const std::vector<CSCTPInfo> * station2_hits = & v_tpInfo[station2];       
    
    //--------------------------------------------------
    // Now that you have both stations, loop over first
    // station hits
    //--------------------------------------------------
	
    for ( int iHit1 = 0; iHit1 < nStation1Hits; ++iHit1 ) {
      
      const CSCTPInfo * station1_hit = & station1_hits -> at ( iHit1 );
      
      //--------------------------------------------------
      // ... and second station hits
      //--------------------------------------------------
      
      for (int iHit2 = 0; iHit2 < nStation2Hits; ++iHit2 ) {
	
	const CSCTPInfo * station2_hit = & station2_hits -> at ( iHit2 );
	  	
	//--------------------------------------------------
	// Get the dphi value
	//--------------------------------------------------

	int dphi = (station1_hit -> gblPhi) - (station2_hit -> gblPhi);
	
	//--------------------------------------------------
	// Don't bother with hits that aren't in the same
	// sector.  Count how many combos fail this cut.
	//--------------------------------------------------

	if ( station1_hit -> sector != station2_hit -> sector ) {
	  m_tree.nskip++;
	  if (m_verbose) std::cout << "         " << "- Combo FAIL: Same sector cut " << std::endl;
	  continue;
	}
	
	//--------------------------------------------------
	// Don't bother with combos outside of our dphi cut
	//
	// Is this the first combination of stations?
	// If it isn't, our dphi cut gets stricter
	//--------------------------------------------------

	int isFirstCombo = 0;

	if ( first_combo_station1 == -999 || first_combo_station2 == -999  ){
	  first_combo_station1 = station1;
	  first_combo_station2 = station2;
	}
	
	if ( station1 == first_combo_station1 || station2 == first_combo_station2 ){	  

	  isFirstCombo = 1;

	  if ( abs(dphi) > 512 ) {
	    if (m_verbose) std::cout <<"         " << "- Combo FAIL: 1st dphi cut " << std::endl;
	    continue;
	  }
	}
	
	else {
	  if ( abs(dphi) > 256 ) {
	    if (m_verbose) std::cout << "         " << "- Combo FAIL: 2nd dphi cut " << std::endl;
	    continue;
	  }
	}
		
	int combo_hitId  = reducedComboHitId ( station1 , station2 );
	int combo_frbId  = reducedComboFRBId ( station1_hit -> frBit , station2_hit -> frBit );	  
	int etaBin       = getEtaBin         ( station1_hit -> eta );
	
	//--------------------------------------------------
	// Store info about this combo in the tree
	//--------------------------------------------------

	m_tree.combo_isFirst[ncombo]    = isFirstCombo;
	m_tree.combo_etaBin [ncombo]    = etaBin;
	m_tree.combo_dphi   [ncombo]    = dphi;
	m_tree.combo_hitId  [ncombo]    = combo_hitId;
	m_tree.combo_frbId  [ncombo]    = combo_frbId;	  
	m_tree.combo_sector [ncombo]    = station1_hit -> sector;
	m_tree.combo_eta    [ncombo][0] = station1_hit -> eta;
	m_tree.combo_eta    [ncombo][1] = station2_hit -> eta;
	m_tree.combo_phi    [ncombo][0] = station1_hit -> gblPhi;	  
	m_tree.combo_phi    [ncombo][1] = station2_hit -> gblPhi;

	if (m_verbose) std::cout << "         " << "+ Combo PASS: dphi = " << dphi << std::endl;
	
	//--------------------------------------------------
	// Fill the histogram
	//--------------------------------------------------
	
	m_plotStorage -> fill( combo_hitId, combo_frbId, m_tree.ptBin, etaBin, isFirstCombo, dphi );
	
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

  m_tree.ncombo = ncombo;

}

int CSCTFMuonPtAnalyzer::reducedComboFRBId ( const int frbit1, const int frbit2 ){

  int id = frbit2 << 1;
  id += frbit1;

  return id;

}

int CSCTFMuonPtAnalyzer::reducedComboHitId ( const int stationId1, const int stationId2){
  
  assert( stationId1 != stationId2 );

  int binary_id = (1<<stationId1);
  binary_id += (1<<stationId2);

  if      (binary_id == 3 ) return 0;
  else if (binary_id == 5 ) return 1;
  else if (binary_id == 6 ) return 2;
  else if (binary_id == 9 ) return 3;
  else if (binary_id == 10) return 4;
  else if (binary_id == 12) return 5;
  else {
    return -999;
  }
}



//--------------------------------------------------
// From a given Correlated LCT digi at a given DetId,
// get the hit eta and phi
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::getHitCoordinates( const CSCDetId& detId, const CSCCorrelatedLCTDigi& digi, 
					     int& lclPhi_phi_bend, int& lclPhi_phi, int& gblPhi_phi, int& gblEta_eta,
					     float& cms_eta, float& cms_phi, bool& bad_phi ){
  
  //--------------------------------------------------
  // Get information about where we are in the detector
  //--------------------------------------------------

  int endcap  = detId.endcap()        - 1; 
  int station = detId.station()       - 1;
  int sector  = detId.triggerSector() - 1;
  int cscId   = detId.triggerCscId()  - 1;
  int chamber = detId.chamber();
  
  bad_phi = (chamber > 8 && station == 3 && sector == 1);

  int subSector = CSCTriggerNumbering::triggerSubSectorFromLabels(detId);
  int fpga      = ( subSector ? subSector-1 : station+1 );
  int endarg    = endcap + 1;

  //--------------------------------------------------
  // Initialize the return values to avoid seg faults
  //--------------------------------------------------

  lclPhi_phi_bend = -990;
  lclPhi_phi      = -990;
  gblPhi_phi      = -990;
  gblEta_eta      = -990;
  cms_eta         = -990.;
  cms_phi         = -990.;

  //--------------------------------------------------
  // Make sure our values are sensible.
  //--------------------------------------------------

  if( endcap<0||endcap>1 || sector<0||sector>6 || station<0||station>3 || cscId<0||cscId>8 || fpga<0||fpga>4) {
    edm::LogError("L1CSCTF: CSC TP are out of range: ") <<"  endcap: "<<(endcap+1)<<"  station: "<<(station+1) <<"  sector: "<<(sector+1)<<"  subSector: "<<subSector<<"  fpga: "<<fpga<<"  cscId: "<<(cscId+1);
      return;
  }

  //--------------------------------------------------
  // Get structs for local and global position vars
  //--------------------------------------------------

  lclphidat lclPhi; // phi within the chamber  
  gblphidat gblPhi; // phi within the sector
  gbletadat gblEta; // eta within the sector

  float cmsPhi, cmsPhiInDegrees;  // consistent with generator phi (CMS phi)
  float cmsEta;                   // consistent with generator eta (CMS eta)
  
  lclPhi = m_srLUTs[fpga][endarg]->localPhi(digi.getStrip(), digi.getPattern(), digi.getQuality(), digi.getBend() );
  gblPhi = m_srLUTs[fpga][endarg]->globalPhiME(lclPhi.phi_local ,digi.getKeyWG(), cscId+1);
  gblEta = m_srLUTs[fpga][endarg]->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local, digi.getKeyWG(), cscId+1);

  //--------------------------------------------------
  // Get the CMS eta value (easy)
  //--------------------------------------------------

  cmsEta = gblEta.global_eta/127. * 1.5 + 0.9;
  
  //--------------------------------------------------
  // The CSC sector phi goes from 0 to 360 degrees
  // Convert to CMS phi, which goes from -pi to pi
  // in radians.  
  // 
  // Phi = 0 is the same for both.
  // Both systems have increasing phi in the same
  // direction.
  //--------------------------------------------------
  
  cmsPhiInDegrees = (gblPhi.global_phi/4096.*62.) - 1. + 15.+60.*sector;  
  if(cmsPhiInDegrees > 360) cmsPhiInDegrees -= 360;

  cmsPhi = cmsPhiInDegrees * TMath::Pi() / 180.0;
  if (cmsPhi > TMath::Pi()) cmsPhi -= (2.0 * TMath::Pi());

  //--------------------------------------------------
  // Pass the values to the output
  //--------------------------------------------------

  cms_phi         = cmsPhi;
  cms_eta         = cmsEta;
  lclPhi_phi_bend = lclPhi.phi_bend_local;
  lclPhi_phi      = lclPhi.phi_local;
  gblPhi_phi      = gblPhi.global_phi;
  gblEta_eta      = gblEta.global_eta;

}

//--------------------------------------------------
// Process for building LUTs within the constructor
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::buildSRLUTs(){
  
  bzero(m_srLUTs,sizeof(m_srLUTs));
  int sector=1; // assume SR LUTs are all same for every sector
  bool TMB07=true; // specific TMB firmware
  // Create a dumy pset for SR LUTs
  edm::ParameterSet srLUTset;
  srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
  srLUTset.addUntrackedParameter<bool>("Binary",   false);
  srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");
  for(int endcap = 1; endcap<=2; endcap++)
    {
      for(int station=1,fpga=0; station<=4 && fpga<5; station++)
	{
	  if(station==1)
	    for(int subSector=0; subSector<2 && fpga<5; subSector++)
	      m_srLUTs[fpga++][endcap] = new CSCSectorReceiverLUT(endcap,
								 sector, subSector+1, station, srLUTset, TMB07);
	  else
	    m_srLUTs[fpga++][endcap] = new CSCSectorReceiverLUT(endcap, sector,
							       0, station, srLUTset, TMB07);
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
// Define framework module
//--------------------------------------------------

DEFINE_FWK_MODULE(CSCTFMuonPtAnalyzer);
