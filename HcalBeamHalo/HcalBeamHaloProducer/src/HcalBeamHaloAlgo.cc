#include "HcalBeamHalo/HcalBeamHaloProducer/interface/HcalBeamHaloAlgo.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include <cmath>

HcalBeamHaloAlgo::HcalBeamHaloAlgo ( const std::vector<int> & rechit_blacklist, double minRecHitEnergy, double minBandEnergy, int minCounts, int width, int maxNWindows, bool verbose ):
  m_rechit_blacklist_raw (rechit_blacklist),
  n_rechit_blacklist (rechit_blacklist.size()),
  m_rechit_minEnergy (minRecHitEnergy),
  m_window_minEnergy (minBandEnergy),
  m_window_minCounts (minCounts),
  m_window_width     (width),
  m_min_iphi         ( 1  ), 
  m_max_iphi         ( 72 ),
  m_max_n_windows    ( maxNWindows ),
  m_verbose          ( verbose )
{
  initialize_rechits();
}

HcalBeamHaloAlgo::~HcalBeamHaloAlgo(){}

//--------------------------------------------------------------
// Geometry setting function
//--------------------------------------------------------------

void HcalBeamHaloAlgo::setGeom ( const CaloGeometry * calo_geometry ){
  m_calo_geometry = calo_geometry;
  m_hb_geometry = calo_geometry -> getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
}

//--------------------------------------------------------------
// Main workhorse function
//--------------------------------------------------------------

void HcalBeamHaloAlgo::process ( const HBHERecHitCollection & hbheRecHits,
				 const reco::TrackCollection & saMuons ){

  //--------------------------------------------------------------
  // Print out status if asked
  //--------------------------------------------------------------
  
  if (m_verbose) {
    std::cout << "Require RecHit HadE >= " << m_rechit_minEnergy << " GeV " << std::endl;
    std::cout << "Require Window HadE >= " << m_window_minEnergy << " GeV " << std::endl;
    std::cout << "Require Window Width = " << m_window_width << " cells " << std::endl;
    std::cout << "Require " << m_window_minCounts << " RecHits for a Window " << std::endl;
    std::cout << "Allow no more than " << m_max_n_windows << " Windows " << std::endl;
  }

  //--------------------------------------------------------------
  // Initialization
  //--------------------------------------------------------------
  
  initialize_blacklists();
  initialize_iphi_strips();
  initialize_windows();
  initialize_halo();
  
  //--------------------------------------------------------------
  // "Track" building
  //--------------------------------------------------------------

  sum_iphi_strips ( hbheRecHits );
  build_windows ();
  build_halo ();
  bool found_fiducial_muon = match_to_sa_muons ( saMuons ) ;
  
}

//--------------------------------------------------------------
// Sum together strips of HCAL cells that are constant in iphi
//--------------------------------------------------------------

void HcalBeamHaloAlgo::sum_iphi_strips ( const HBHERecHitCollection & hbheRecHits ){
  
  if (m_verbose) std::cout << " Summing iphi strips " << std::endl;

  HBHERecHitCollection::const_iterator hbheRecHit     = hbheRecHits.begin();
  HBHERecHitCollection::const_iterator hbheRecHit_end = hbheRecHits.end();

  int nrechit = 0;
  
  for (; hbheRecHit != hbheRecHit_end; ++hbheRecHit ){
    float           energy = hbheRecHit -> energy();
    HcalDetId       id     = hbheRecHit -> id();
    HcalSubdetector subdet = id.subdet();
    
    bool blacklisted = check_rechit_blacklist ( id ) ;

    if (subdet != HcalBarrel       ) continue;
    if (energy < m_rechit_minEnergy) continue;
    if (blacklisted                ) continue;
    
    int iphi = id.iphi();

    m_iphi_strip_cells [iphi].push_back(id);
    m_iphi_strip_energy[iphi] += energy;

    nrechit++;

  }

  if ( !m_verbose ) return;

  for ( int iphi = m_min_iphi; iphi <= m_max_iphi; ++iphi ){
    if ( m_iphi_strip_cells[iphi].size() < 1 ) continue;
    std::cout << "  Found " << m_iphi_strip_cells[iphi].size() << " RecHits in iphi = " << iphi << " with combined energy = " << m_iphi_strip_energy[iphi] << " GeV " << std::endl;
  }

  if (m_verbose) std::cout << " Summed  iphi strips " << std::endl;

}

//--------------------------------------------------------------
// Select and build sliding windows from strips of HCAL cells
// that are constant in iphi
//--------------------------------------------------------------

void HcalBeamHaloAlgo::build_windows(){

  if (m_verbose) std::cout << " Building sliding windows " << std::endl;

  //--------------------------------------------------------------
  // Look at prescribed number of sliding windows
  //--------------------------------------------------------------

  for ( int iwindow = 0; iwindow < m_max_n_windows ; ++iwindow ){

    //--------------------------------------------------------------
    // There are 72 possibilities for each sliding window
    //--------------------------------------------------------------

    float best_window_energy = -1.0;
    int   best_window_counts = -1;
    int   best_window_min_iphi_index_withBL  =  999;
    int   best_window_max_iphi_index_withBL  = -999;
    int   best_window_min_iphi_index         =  999;
    int   best_window_max_iphi_index         = -999;

    for ( int window_index = m_min_iphi ; window_index <= m_max_iphi ; ++window_index) {

      //--------------------------------------------------------------
      // If this window is blacklisted, skip it
      //--------------------------------------------------------------

      bool blacklisted_window = check_window_blacklist ( window_index ) ;
      if ( blacklisted_window ) continue;

      //--------------------------------------------------------------
      // Each window has a minimum and a maximum iphi
      //--------------------------------------------------------------

      int min_iphi_index = window_index;
      int max_iphi_index = window_index + m_window_width - 1;
      
      //--------------------------------------------------------------
      // Loop over the iphi strips in each window and sum the energies
      //--------------------------------------------------------------

      float window_energy = 0.0;
      int   window_counts = 0;

      for (int iphi_index = min_iphi_index; iphi_index <= max_iphi_index; ++iphi_index ){

	int iphi = iphi_index;
	if (iphi > m_max_iphi) iphi -= m_max_iphi;

	bool is_blacklisted_iphi = check_iphi_blacklist ( iphi );
	if ( is_blacklisted_iphi ) continue;

	window_energy += m_iphi_strip_energy [ iphi ];
	window_counts += m_iphi_strip_cells  [ iphi ].size();
	
      }	

      //--------------------------------------------------------------
      // If this window doesn't pass our cuts, skip it
      //--------------------------------------------------------------

      if ( window_energy < m_window_minEnergy ||
	   window_counts < m_window_minCounts ) continue;
      
      //--------------------------------------------------------------
      // If this is our best window, save it
      //--------------------------------------------------------------
      
      if ( window_energy > best_window_energy ){
	best_window_energy          = window_energy;
	best_window_counts          = window_counts;
	best_window_min_iphi_index  = min_iphi_index;
	best_window_max_iphi_index  = max_iphi_index;
      }
    }

    //--------------------------------------------------------------
    // Now that we know the best window, add its iphi strips to the
    // iphi blacklist
    //--------------------------------------------------------------

    for (int iphi_index = best_window_min_iphi_index; iphi_index <= best_window_max_iphi_index; ++iphi_index ){
      int iphi = iphi_index;
      if (iphi > m_max_iphi) iphi -= m_max_iphi;
      
      bool is_in_iphi_blacklist = check_iphi_blacklist ( iphi ) ;
      if ( is_in_iphi_blacklist ) continue;

      if ( iphi < best_window_min_iphi_index_withBL ) best_window_min_iphi_index_withBL = iphi_index;
      if ( iphi > best_window_max_iphi_index_withBL ) best_window_max_iphi_index_withBL = iphi_index;
      
      add_to_iphi_blacklist (iphi);
    }

    //--------------------------------------------------------------
    // Save the properties of the best window
    //--------------------------------------------------------------
    
     if (best_window_energy > 0) {
       m_best_window_min_iphi_indices.push_back ( best_window_min_iphi_index_withBL );
       m_best_window_max_iphi_indices.push_back ( best_window_max_iphi_index_withBL );
       m_best_window_energies.push_back         ( best_window_energy         );
       if (m_verbose) 
	 std::cout << "  Found " << best_window_counts << " RecHits in window with range = ( " << best_window_min_iphi_index_withBL << ", " << best_window_max_iphi_index_withBL << "), with combined energy = " << best_window_energy << " GeV " << std::endl;
     }

    //--------------------------------------------------------------
    // Add the window to the window blacklsit
    //--------------------------------------------------------------

    add_to_window_blacklist ( best_window_min_iphi_index );

  }

  if (m_verbose) std::cout << " Built    sliding windows " << std::endl;
  
}

//--------------------------------------------------------------
// Build the actual halo objects from the sliding windows
//--------------------------------------------------------------

void HcalBeamHaloAlgo::build_halo (){
  
  int n_windows = m_best_window_energies.size();
  
  if (m_verbose) {
    std::cout << " Building halo " << std::endl;
    std::cout << " Found " << n_windows << " halo candidates " << std::endl;
  }
  
  std::vector <HcalDetId> ids;

  for ( int iwindow = 0; iwindow < n_windows; ++iwindow ){
    
    ids.clear();

    int window_min_iphi_index = m_best_window_min_iphi_indices[iwindow];
    int window_max_iphi_index = m_best_window_max_iphi_indices[iwindow];
    
    int window_ieta   = 1;
    int window_depth  = 1;
    int window_iphi_1 = window_min_iphi_index;
    int window_iphi_2 = window_max_iphi_index;
    float window_energy = m_best_window_energies [iwindow];

    if (window_iphi_1 > m_max_iphi) window_iphi_1 -= m_max_iphi;
    if (window_iphi_2 > m_max_iphi) window_iphi_2 -= m_max_iphi;
    
    for ( int iphi_index = window_min_iphi_index ; iphi_index <= window_max_iphi_index; ++iphi_index ){
      
      int iphi = iphi_index;
      if (iphi > m_max_iphi ) iphi -= m_max_iphi;

      ids.insert ( ids.end(), m_iphi_strip_cells[iphi].begin(), m_iphi_strip_cells[iphi].end() );

    }

    HcalDetId cell_1 = HcalDetId ( HcalBarrel, window_ieta, window_iphi_1 , window_depth );
    HcalDetId cell_2 = HcalDetId ( HcalBarrel, window_ieta, window_iphi_2 , window_depth );
	  
    const CaloCellGeometry* cell_1_geometry = m_hb_geometry -> getGeometry(cell_1);
    const CaloCellGeometry* cell_2_geometry = m_hb_geometry -> getGeometry(cell_2);

    float phi_1 = cell_1_geometry -> getPosition().phi();
    float phi_2 = cell_2_geometry -> getPosition().phi();

    float half_cell_width = 2.5 * 3.141592 / 180.0;

    phi_1 -= half_cell_width;
    phi_2 += half_cell_width;

    if (m_verbose){
      std::cout << "  Inserting " << ids.size() << " RecHits into halo candidate #" << iwindow + 1 << std::endl;
      std::cout << "    phi_1 = " << phi_1 << std::endl;
      std::cout << "    phi_2 = " << phi_2 << std::endl;
    }

    HcalBeamHalo halo;    
    halo.setConstituents ( ids );
    halo.setPhiBoundaries ( window_iphi_1, window_iphi_2, phi_1, phi_2 );
    halo.setEnergy ( window_energy );

    m_beam_halo.push_back(halo);

  }  
  
  if (m_verbose) std::cout << " Built    halo " << std::endl;
  
}

//--------------------------------------------------------------
// Try to match the halo objects to StandAlone muons
//--------------------------------------------------------------

bool HcalBeamHaloAlgo::match_to_sa_muons ( const reco::TrackCollection & saMuons ){

  reco::TrackCollection::const_iterator saMuon     = saMuons.begin();
  reco::TrackCollection::const_iterator saMuon_end = saMuons.end();

  int nSaMuon = 0;
  int nHaloWindows = m_beam_halo.size();

  bool fiducial_muon = false;
  
  for (; saMuon != saMuon_end; ++saMuon){
    
    float mu_phi_plus, mu_phi_minus;
    
    if (m_verbose) std::cout << " Found StandAlone muon #" << nSaMuon + 1 << std::endl;

    bool isFiducial = isHBFiducial ( * saMuon, mu_phi_plus, mu_phi_minus );

    if (!isFiducial) continue;

    fiducial_muon = true;

    for ( int iHalo = 0; iHalo < nHaloWindows; ++iHalo ){

      bool matched = match ( mu_phi_plus, mu_phi_minus, m_beam_halo[iHalo].maxPhi(), m_beam_halo[iHalo].minPhi());
      m_beam_halo[iHalo].setMatchedStatus ( matched ) ;

    }
      
    nSaMuon++;

  }

  return fiducial_muon;

}

//--------------------------------------------------------------
// saMuons functions
//--------------------------------------------------------------

void HcalBeamHaloAlgo::add_to_window_blacklist (int window){

  bool is_in_window_blacklist = check_window_blacklist ( window ) ;
  if (!is_in_window_blacklist ) m_blacklist_window.push_back (window) ;

}

void HcalBeamHaloAlgo::add_to_iphi_blacklist ( int iphi ){

  bool is_in_iphi_blacklist = check_iphi_blacklist ( iphi ) ;
  if (!is_in_iphi_blacklist) m_blacklist_iphi.push_back ( iphi );

}

bool HcalBeamHaloAlgo::check_window_blacklist (int window ){
  
  std::vector<int>::iterator iter     = m_blacklist_window.begin();
  std::vector<int>::iterator iter_end = m_blacklist_window.end();
  
  bool is_in_blacklist = false;
  
  for (; iter != iter_end; ++iter)
    is_in_blacklist |= ((*iter) == window );
  
  return is_in_blacklist;

}

bool HcalBeamHaloAlgo::check_iphi_blacklist  ( int iphi ){

  std::vector<int>::iterator iter     = m_blacklist_iphi.begin();
  std::vector<int>::iterator iter_end = m_blacklist_iphi.end();

  bool is_in_blacklist = false;

  for (; iter != iter_end; ++iter)
    is_in_blacklist |= ((*iter) == iphi );
  
  return is_in_blacklist;

}

bool HcalBeamHaloAlgo::check_rechit_blacklist ( const HcalDetId & id ) {

  bool blacklisted = false;

  int test_ieta  = id.ieta();
  int test_iphi  = id.iphi();
  int test_depth = id.depth();

  for (int iRechit = 0 ; iRechit < n_rechit_blacklist; ++iRechit ) {
        
    int bad_ieta   = m_rechit_blacklist[iRechit].ieta();
    int bad_iphi   = m_rechit_blacklist[iRechit].iphi();
    int bad_depth  = m_rechit_blacklist[iRechit].depth();    

    bool match = ( test_ieta  == bad_ieta &&
		   test_iphi  == bad_iphi &&
		   test_depth == bad_depth );
    
    blacklisted |= match;

  }

  return blacklisted;

}

//--------------------------------------------------------------
// Initialization functions
//--------------------------------------------------------------


void HcalBeamHaloAlgo::initialize_windows(){

  if (m_verbose) std::cout << " Initializing windows " << std::endl;

  m_best_window_min_iphi_indices.clear();
  m_best_window_max_iphi_indices.clear();
  m_best_window_energies.clear();

  if (m_verbose) std::cout << " Initialized  windows " << std::endl;

}

void HcalBeamHaloAlgo::initialize_iphi_strips(){

  if (m_verbose) std::cout << " Initialazing iphi strips " << std::endl;

  m_iphi_strip_energy = std::vector<float> (m_max_iphi + 1,0.0);
  m_iphi_strip_cells  = std::vector < std::vector < HcalDetId > > ( m_max_iphi + 1 );

  if (m_verbose) std::cout << " Initialazed iphi strips " << std::endl;

}

void HcalBeamHaloAlgo::initialize_blacklists(){

  if (m_verbose) std::cout << " Initializing blacklists " << std::endl;
  m_blacklist_iphi.clear();
  m_blacklist_window = std::vector< int> ( m_max_n_windows, 0 );
  if (m_verbose) std::cout << " Initialized  blacklists " << std::endl;
}

bool HcalBeamHaloAlgo::isHBFiducial ( const reco::Track & track, float & phi_plus, float & phi_minus ){

  float hb_z_min = -570.0;
  float hb_z_max =  570.0;
  float hb_r_min =  177.0;
  float hb_r_max =  295.0;

  float x1 = track.innerPosition().x();
  float y1 = track.innerPosition().y();
  float z1 = track.innerPosition().z();
  float r1 = sqrt ( ( x1 * x1 ) + ( y1 * y1 ) );
  float phi1 = atan2 ( y1, x1 );
  
  float x2 = track.outerPosition().x();
  float y2 = track.outerPosition().y();
  float z2 = track.outerPosition().z(); 
  float r2 = sqrt ( ( x2 * x2 ) + ( y2 * y2 ) );
  float phi2 = atan2 ( y2, x2 );
  
  float mx = ( x1 - x2 ) / ( z1 - z2 );
  float my = ( y1 - y2 ) / ( z1 - z2 );

  float dz_plus  = hb_z_max - z1;
  float dz_minus = hb_z_min - z1;

  float bx = x1;
  float by = y1;

  float x_hb_minus_edge = mx * dz_minus + bx;
  float x_hb_plus_edge  = mx * dz_plus  + bx; 
  float y_hb_minus_edge = my * dz_minus + by;
  float y_hb_plus_edge  = my * dz_plus  + by; 

  float phi_hb_minus_edge = atan2( y_hb_minus_edge, x_hb_minus_edge );
  float phi_hb_plus_edge  = atan2( y_hb_plus_edge , x_hb_plus_edge  );
  
  float rad_hb_minus_edge = sqrt ( ( x_hb_minus_edge * x_hb_minus_edge ) +
				   ( y_hb_minus_edge * y_hb_minus_edge ) );
  
  float rad_hb_plus_edge  = sqrt ( ( x_hb_plus_edge  * x_hb_plus_edge  ) +
				   ( y_hb_plus_edge  * y_hb_plus_edge  ) );
  
  bool hits_minus_edge = ( rad_hb_minus_edge >= hb_r_min && rad_hb_minus_edge <= hb_r_max );
  bool hits_plus_edge  = ( rad_hb_plus_edge  >= hb_r_min && rad_hb_plus_edge  <= hb_r_max );

  bool hits_outer_edge = ( ( rad_hb_minus_edge > hb_r_max &&      // Over the outer edge on the - side
			     rad_hb_plus_edge  < hb_r_max ) ||    //    ... under the outer edge on the + side
			   ( rad_hb_minus_edge < hb_r_max &&      // Under the outer edge on the - side
			     rad_hb_plus_edge  > hb_r_max ) ) ;   //    ... over the outer edge on the + side

  bool hits_inner_edge = ( ( rad_hb_minus_edge > hb_r_min &&      // Over the inner edge on the - side
			     rad_hb_plus_edge  < hb_r_min ) ||    //    ... under the inner edge on the + side
			   ( rad_hb_minus_edge < hb_r_min  &&     // Under the inner edge on the - side
			     rad_hb_plus_edge  > hb_r_min ) ) ;   //    ... over the inner edge on the + side 

  phi_minus = phi_hb_minus_edge;
  phi_plus  = phi_hb_plus_edge ;
  
  if (m_verbose){
    std::cout << "   InnerP: x = " << x1 << ", y = " << y1 << ", z = " << z1 << ", r = " << r1 << ", phi = " << phi1 << std::endl;
    std::cout << "   OuterP: x = " << x2 << ", y = " << y2 << ", z = " << z2 << ", r = " << r2 << ", phi = " << phi2 << std::endl;
    std::cout << "   HB +  : x = " << x_hb_plus_edge << ", y = " << y_hb_plus_edge << ", z = " << hb_z_max << ", r = " << rad_hb_plus_edge << ", phi = " << phi_hb_plus_edge << std::endl;
    std::cout << "   HB -  : x = " << x_hb_minus_edge << ", y = " << y_hb_minus_edge << ", z = " << hb_z_min << ", r = " << rad_hb_minus_edge << ", phi = " << phi_hb_minus_edge << std::endl;
    
    if (hits_minus_edge) std::cout << "    Hits minus edge!" << std::endl;
    if (hits_plus_edge ) std::cout << "    Hits plus  edge!" << std::endl;
    if (hits_inner_edge) std::cout << "    Hits inner edge!" << std::endl;
    if (hits_outer_edge) std::cout << "    Hits outer edge!" << std::endl;
  }

  bool retval = ( hits_minus_edge || 
		  hits_plus_edge  || 
		  hits_inner_edge || 
		  hits_outer_edge  );

  if (m_verbose) if (!retval) std::cout << "    NO HITS!" << std::endl;

  return retval;
			   


}

bool HcalBeamHaloAlgo::match ( float mu_phi_plus, float mu_phi_minus , float window_phi_plus, float window_phi_minus ){

  bool retval = false;

  return retval;

}

void HcalBeamHaloAlgo::getHalo ( std::auto_ptr <HcalBeamHaloCollection> & coll ){
  
  (*coll).insert((*coll).begin(), m_beam_halo.begin() , m_beam_halo.end());

}

void HcalBeamHaloAlgo::initialize_halo () {
  m_beam_halo.clear();
}

void HcalBeamHaloAlgo::initialize_rechits(){

  std::vector<int>::iterator raw_id_ptr     = m_rechit_blacklist_raw.begin();
  std::vector<int>::iterator raw_id_ptr_end = m_rechit_blacklist_raw.end();

  for (; raw_id_ptr != raw_id_ptr_end; ++raw_id_ptr){
    
    int raw_id = * raw_id_ptr;

    int depth = raw_id - (raw_id / 100) * 100;
    int iphi  = ( raw_id / 100 ) - ( raw_id / 10000 ) * 100;
    int ieta  = ( raw_id / 10000 );
    
    HcalDetId id = HcalDetId ( HcalBarrel, ieta, iphi , depth ) ;

    m_rechit_blacklist.push_back ( id );

  }

}
