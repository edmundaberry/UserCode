#ifndef HCALBEAMHALO_HCALBEAMHALOPRODUCER_HCALBEAMHALOALGO_H
#define HCALBEAMHALO_HCALBEAMHALOPRODUCER_HCALBEAMHALOALGO_H

#include <string>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "HcalBeamHalo/HcalBeamHaloProducer/interface/HcalBeamHaloCollection.h"

class HcalBeamHaloAlgo {

 public:
  HcalBeamHaloAlgo( const std::vector<int> & rechit_blacklist, double minRecHitEnergy, double minBandEnergy, int minCounts, int width, int maxNWindows, bool verbose );
  ~HcalBeamHaloAlgo();
  
  void process ( const HBHERecHitCollection & hbheRecHits, const reco::TrackCollection & saMuons );
  void setGeom ( const CaloGeometry * calo_geometry );
  void getHalo ( std::auto_ptr <HcalBeamHaloCollection> & coll );

 private:

  //--------------------------------------------------------------
  // Functions for blacklisting a phi strip, window, or rechit
  //--------------------------------------------------------------
  
  void initialize_blacklists   ();
  void add_to_iphi_blacklist   ( int iphi   );
  void add_to_window_blacklist ( int window );
  bool check_iphi_blacklist    ( int iphi   );
  bool check_window_blacklist  ( int window );
  bool check_rechit_blacklist  ( const HcalDetId & id );

  //--------------------------------------------------------------
  // Functions for building the halo tracks
  //--------------------------------------------------------------

  void initialize_iphi_strips();
  void initialize_windows();
  void initialize_halo ();
  void initialize_rechits();
  void sum_iphi_strips ( const HBHERecHitCollection & hbheRecHits );
  void build_windows();
  void build_halo();
  bool match_to_sa_muons ( const reco::TrackCollection & saMuons );
  bool match ( float mu_phi_plus, float mu_phi_minus , float window_phi_plus, float window_phi_minus );

  //--------------------------------------------------------------
  // Helper functions (fiduciality, etc)
  //--------------------------------------------------------------

  bool isHBFiducial ( const reco::Track & track, float & phi_plus, float & phi_minus );

  //--------------------------------------------------------------
  // User-set constant global variables
  //--------------------------------------------------------------

  std::vector <int> m_rechit_blacklist_raw;
  std::vector <HcalDetId> m_rechit_blacklist;

  const int    n_rechit_blacklist;
  const double m_rechit_minEnergy;
  const double m_window_minEnergy;
  const int    m_window_minCounts;
  const int    m_window_width    ;
  const int    m_min_iphi        ;
  const int    m_max_iphi        ;
  const int    m_max_n_windows   ;

  //--------------------------------------------------------------
  // Verbosity
  //--------------------------------------------------------------

  const bool m_verbose;
  
  //--------------------------------------------------------------
  // Geometry 
  //--------------------------------------------------------------

  const CaloGeometry * m_calo_geometry;
  const CaloSubdetectorGeometry * m_hb_geometry;
  
  //--------------------------------------------------------------
  // HcalBeamHalo vector
  //--------------------------------------------------------------

  std::vector < HcalBeamHalo > m_beam_halo; 

  //--------------------------------------------------------------
  // Window vectors
  //--------------------------------------------------------------

  std::vector <int> m_best_window_min_iphi_indices;
  std::vector <int> m_best_window_max_iphi_indices;
  std::vector <int> m_best_window_energies        ;

  //--------------------------------------------------------------
  // Blacklist vectors
  //--------------------------------------------------------------

  std::vector <int> m_blacklist_iphi;
  std::vector <int> m_blacklist_window;

  //--------------------------------------------------------------
  // iphi strip vectors
  //--------------------------------------------------------------

  std::vector<float> m_iphi_strip_energy;
  std::vector < std::vector < HcalDetId > > m_iphi_strip_cells;

   
};

#endif
