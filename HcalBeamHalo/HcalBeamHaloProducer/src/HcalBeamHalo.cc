#include "HcalBeamHalo/HcalBeamHaloProducer/interface/HcalBeamHalo.h"

HcalBeamHalo::HcalBeamHalo ():
  m_min_iphi (-999),
  m_max_iphi (-999),
  m_min_phi  (-999.),
  m_max_phi  (-999.),
  m_energy   (-999.),
  m_matched  ( false )
{}

HcalBeamHalo::HcalBeamHalo ( const HcalBeamHalo& rhs ) :
  m_constituents (rhs.m_constituents),
  m_min_iphi     (rhs.m_min_iphi    ), 
  m_max_iphi     (rhs.m_max_iphi    ), 
  m_min_phi      (rhs.m_min_phi     ), 
  m_max_phi      (rhs.m_max_phi     ),
  m_energy       (rhs.m_energy      ),
  m_matched      (rhs.m_matched     )
{}

HcalBeamHalo::~HcalBeamHalo(){}

void HcalBeamHalo::setPhiBoundaries ( int min_iphi, int max_iphi,
				      float min_phi, float max_phi ){

  m_min_iphi = min_iphi;
  m_max_iphi = max_iphi;
  m_min_phi  = min_phi ;
  m_max_phi  = max_phi ;

}
