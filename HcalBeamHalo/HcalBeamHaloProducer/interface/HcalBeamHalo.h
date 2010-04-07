#ifndef HCALBEAMHALO_HCALBEAMHALOPRODUCER_HCALBEAMHALO_H
#define HCALBEAMHALO_HCALBEAMHALOPRODUCER_HCALBEAMHALO_H

#include <vector>

#include "DataFormats/HcalDetId/interface/HcalDetId.h"

class HcalBeamHalo {

 public:

  //------------------------------------------------------
  // Constructor/Destructor
  //------------------------------------------------------
  
  HcalBeamHalo();
  HcalBeamHalo ( const HcalBeamHalo& rhs );
  //  HcalBeamHalo( std::vector < HcalDetId > v ); 
  ~HcalBeamHalo();

  //------------------------------------------------------
  // User functions
  //------------------------------------------------------
  
  void setPhiBoundaries ( int min_iphi, int max_iphi,
			  float min_phi, float max_phi );

  float minPhi  () { return m_min_phi  ; }
  float maxPhi  () { return m_max_phi  ; }
  int   minIPhi () { return m_min_iphi ; }
  int   maxIPhi () { return m_max_iphi ; }

  const std::vector <HcalDetId> & getConstituents() { return m_constituents; }
  void addConstituent  ( HcalDetId id ) { m_constituents.push_back ( id ); }
  void setConstituents ( std::vector <HcalDetId> v ) { m_constituents = v; }

  bool getMatchedStatus () { return m_matched; }
  void setMatchedStatus (bool m) { m_matched = m; }

  float getEnergy() { return m_energy; }
  void  setEnergy( float e ) { m_energy = e ; } 

  //------------------------------------------------------
  // Private data
  //------------------------------------------------------
  
 private:
  
  std::vector <HcalDetId> m_constituents;
  
  int   m_min_iphi;
  int   m_max_iphi;
  float m_min_phi;
  float m_max_phi;
  float m_energy ; 
  bool  m_matched;
  
};

#endif
