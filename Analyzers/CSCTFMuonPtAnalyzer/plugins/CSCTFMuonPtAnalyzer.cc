#include <memory>
#include <stdio.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCStatusDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "Analyzers/CSCTFMuonPtAnalyzer/interface/CSCTFMuonPtAnalyzer.h"

CSCTFMuonPtAnalyzer::CSCTFMuonPtAnalyzer(const edm::ParameterSet& iConfig):
  m_cscCorrLCTDigiTag ( iConfig.getParameter<edm::InputTag>("CSCCorrLCTDigiTag")),
  m_genParticlesTag   ( iConfig.getParameter<edm::InputTag>("GenParticlesTag")),
  m_analyzeGenMuons   ( iConfig.getParameter<bool>("AnalyzeGenMuons")),
  m_analyzeDataMuons  ( iConfig.getParameter<bool>("AnalyzeDataMuons")),
  m_nPtBins           ( iConfig.getParameter<int> ("NPtBins")),
  m_minBinnablePt     ( iConfig.getParameter<double> ("MinBinnablePt")),
  m_maxBinnablePt     ( iConfig.getParameter<double> ("MaxBinnablePt")),
  m_plotStorage ( MuonBinnedPlotStorage(6,m_nPtBins, "test.root" )),
  m_sector(1),
  m_TMB07(true)
{
  
  buildSRLUTs();
  buildPtBins();
  
}

CSCTFMuonPtAnalyzer::~CSCTFMuonPtAnalyzer() {}

void CSCTFMuonPtAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //--------------------------------------------------
  // Assign event and setup pointers
  //--------------------------------------------------

  m_event = &iEvent;
  m_setup = &iSetup;
  
  //--------------------------------------------------
  // Analyze "real" muons, whether from data or MC
  //--------------------------------------------------

  if (m_analyzeGenMuons ) analyzeGenMuons ();
  if (m_analyzeDataMuons) analyzeDataMuons();

  //--------------------------------------------------
  // Analyze the trigger primitives
  //--------------------------------------------------

  analyzeCSCTrigPrims();

  //--------------------------------------------------
  // Free up memory
  //--------------------------------------------------

  m_realMuonList.clear();
  
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
  
  for (; genParticle != genParticle_end; ++genParticle){

    int pdg = (int) genParticle -> pdgId();

    if ( abs(pdg) != 13 ) continue;
    
    reco::Particle::LorentzVector p4 = genParticle -> p4();
    
    double gen_pt  = p4.pt();
    double gen_phi = p4.phi();
    double gen_eta = p4.eta();

    int gen_ptBin = getPtBin ( gen_pt );

    if ( gen_ptBin < 0 || gen_ptBin == m_nPtBins ) continue;

    m_realMuonList.push_back(RealMuon (gen_eta, gen_phi, gen_pt, gen_ptBin));

  }
  
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
  // Loop over the LCT digis
  // Note that the DigiRangeIterator represents a pair
  //    -- First  element is the DetId
  //    -- Second element is a **collection** of 
  //       associated digis.  May be more than one!
  //--------------------------------------------------

  CSCCorrelatedLCTDigiCollection::DigiRangeIterator corrLCT     = CorrLCTs.product() -> begin();
  CSCCorrelatedLCTDigiCollection::DigiRangeIterator corrLCT_end = CorrLCTs.product() -> end();

  std::vector<float> station_average_phi(5,0);
  std::vector<float> station_average_eta(5,0);
  std::vector<int  > station_num_hits   (5,0);

  for (; corrLCT != corrLCT_end; ++corrLCT){
    
    //--------------------------------------------------
    // For each DetId, get the range of digis
    //--------------------------------------------------

    CSCDetId * cscDetId = &((*corrLCT).first);

    int station = cscDetId -> station() - 1;

    CSCCorrelatedLCTDigiCollection::Range corrLCTDigi_range = CorrLCTs.product() -> get( *cscDetId  );    
    CSCCorrelatedLCTDigiCollection::const_iterator corrLCTDigi     = corrLCTDigi_range.first;
    CSCCorrelatedLCTDigiCollection::const_iterator corrLCTDigi_end = corrLCTDigi_range.second;
    
    //--------------------------------------------------
    // Now loop over the LCT's
    //--------------------------------------------------
    
    float detId_average_phi = 0;
    float detId_average_eta = 0;
    int   ndigi = 0;

    for (; corrLCTDigi != corrLCTDigi_end; ++corrLCTDigi ){    
      
      float hit_eta, hit_phi;
      
      getHitCoordinates( *cscDetId , *corrLCTDigi, hit_eta, hit_phi );

      detId_average_phi += hit_phi;
      detId_average_eta += hit_eta;
      ++ndigi;
    }

    detId_average_phi /= ((float) ndigi );
    detId_average_eta /= ((float) ndigi );

    station_average_phi[station] += detId_average_phi;
    station_average_eta[station] += detId_average_eta;
    station_num_hits   [station] ++;
    
  }
  
  int combo_binary_id = 0;
  int num_station_hits = 0;

  for (int iStation = 0; iStation <= 3; iStation++){
    if (station_num_hits[iStation] != 0){
      num_station_hits++;
      combo_binary_id += (1<<iStation);
      station_average_phi[iStation] /= ((float) station_num_hits[iStation]);
      station_average_eta[iStation] /= ((float) station_num_hits[iStation]);
    }
  }  
  
  //--------------------------------------------------
  // Fill the histograms appropriately
  //--------------------------------------------------

  int ptBin   = m_realMuonList[0].ptBin;

  if (num_station_hits >= 2) fillHist ( station_average_phi, station_average_eta, combo_binary_id, ptBin);
  
}

void CSCTFMuonPtAnalyzer::fillHist ( const std::vector<float> & station_average_phi,
				     const std::vector<float> & station_average_eta,
				     const int combo_binary_id,
				     const int ptBin ){

  for (int iHit1 = 0; iHit1 < 4; iHit1++){
    for (int iHit2 = iHit1 + 1; iHit2 < 4; iHit2++){
      
      bool hit_1st_station = (combo_binary_id == (combo_binary_id | (1 << iHit1)));
      bool hit_2nd_station = (combo_binary_id == (combo_binary_id | (1 << iHit2)));

      if (!hit_1st_station) continue;
      if (!hit_2nd_station) continue;      

      double deltaPhi = reco::deltaPhi( station_average_phi[iHit1],station_average_phi[iHit2]);
      int thisReducedComboId = reducedComboId (iHit1, iHit2);

//std::cout << "phi1 = " << station_average_phi[iHit1] << ", " 
//		<< "phi2 = " << station_average_phi[iHit2] << ", " 
//		<< "deltaPhi = " << deltaPhi << std::endl;
      
      m_plotStorage.fill( thisReducedComboId , ptBin, deltaPhi );
             
    }
  }
  
}

int CSCTFMuonPtAnalyzer::reducedComboId ( const int stationId1, const int stationId2){
  
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
					     float& eta, float& phi ){
  
  int endcap  = detId.endcap()        - 1; 
  int station = detId.station()       - 1;
  int sector  = detId.triggerSector() - 1;
  int cscId   = detId.triggerCscId()  - 1;
  
  int subSector = CSCTriggerNumbering::triggerSubSectorFromLabels(detId);
  int fpga      = ( subSector ? subSector-1 : station+1 );
  int endarg    = endcap + 1;
  
  lclphidat lclPhi; // phi within the chamber
  
  gblphidat gblPhi; // phi within the sector
  gbletadat gblEta; // eta within the sector
  
  float sectorPhi;  // consistent with generator phi (CMS phi)
  float sectorEta;  // consistent with generator eta (CMS eta)
  
  lclPhi = m_srLUTs[fpga][endarg]->localPhi(digi.getStrip(), digi.getPattern(), digi.getQuality(), digi.getBend() );
  
  gblPhi = m_srLUTs[fpga][endarg]->globalPhiME(lclPhi.phi_local ,digi.getKeyWG(), cscId+1);
  gblEta = m_srLUTs[fpga][endarg]->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local, digi.getKeyWG(), cscId+1);
  
  sectorEta = gblEta.global_eta/127. * 1.5 + 0.9;
  sectorPhi = (gblPhi.global_phi/4096.*62.) - 1. + 15.+60.*sector;
  if(sectorPhi > 360) sectorPhi -= 360;
  
  phi = sectorPhi;
  eta = sectorEta;
  
}

//--------------------------------------------------
// Process for building LUTs within the constructor
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::buildSRLUTs(){
  
  //--------------------------------------------------
  // Create a dummy parameter set for the SR LUTs
  //--------------------------------------------------
  
  edm::ParameterSet srLUTset;

  srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
  srLUTset.addUntrackedParameter<bool>("Binary",   false);
  srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");

  //--------------------------------------------------
  // Loop over all components of CSC & assign LUTs
  //--------------------------------------------------

  for(int endcap = 1; endcap<=2; ++endcap) {
    
    for(int station=1,fpga=0; station<=4 && fpga<5; ++station) {
      if(station==1)
	for(int subSector=0; subSector<2 && fpga<5; subSector++)
	  m_srLUTs[fpga++][endcap] = new CSCSectorReceiverLUT(endcap, m_sector, subSector+1, station, srLUTset, m_TMB07);
      else
	m_srLUTs[fpga++][endcap] = new CSCSectorReceiverLUT(endcap, m_sector,0, station, srLUTset, m_TMB07);
    }
  }
}

//--------------------------------------------------
// Pt bins
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::buildPtBins(){

  double nPtBins = (double) m_nPtBins;
  double step = (m_maxBinnablePt - m_minBinnablePt) / nPtBins;

  for (int iPtBin = 0; iPtBin <= nPtBins; ++iPtBin){
    double boundary = ((double) iPtBin) * step;
    m_ptBins.push_back(boundary);
  }

}

int CSCTFMuonPtAnalyzer::getPtBin( double pt ){

  if ( pt > m_maxBinnablePt ) return m_nPtBins + 1;
  if ( pt < m_minBinnablePt ) return -1;
  
  int size = m_ptBins.size();
  
  for (int iPtBin = 0; iPtBin < size - 1; ++iPtBin){
    double bin_min_pt = m_ptBins[iPtBin];
    double bin_max_pt = m_ptBins[iPtBin+1];
    
    if (pt > bin_min_pt && pt < bin_max_pt) return iPtBin;
  }

  return -2;

}

//--------------------------------------------------
// Muon storage struct initializer
//--------------------------------------------------

CSCTFMuonPtAnalyzer::RealMuon::RealMuon(double temp_eta, double temp_phi, double temp_pt, int temp_ptBin):
  eta(temp_eta), phi(temp_phi), pt(temp_pt), ptBin(temp_ptBin){}

//--------------------------------------------------
// Beginning and ending job functions
//--------------------------------------------------

void CSCTFMuonPtAnalyzer::beginJob(){}
void CSCTFMuonPtAnalyzer::endJob()  {
  m_plotStorage.save();
}

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
