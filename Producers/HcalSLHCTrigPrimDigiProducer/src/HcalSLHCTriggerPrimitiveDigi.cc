#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTriggerPrimitiveDigi.h"

//------------------------------------------------------
// Constructors/destructor
//------------------------------------------------------

HcalSLHCTriggerPrimitiveDigi::HcalSLHCTriggerPrimitiveDigi(){}
HcalSLHCTriggerPrimitiveDigi::~HcalSLHCTriggerPrimitiveDigi(){}
HcalSLHCTriggerPrimitiveDigi::HcalSLHCTriggerPrimitiveDigi(const HcalTrigTowerDetId& id):
  m_id  (id),
  m_size(0),
  m_hcalPresamples(0)
{}

//------------------------------------------------------
// Set methods
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveDigi::setSize(int size) {
  if (size<0) m_size=0;
  else if (size>MAXSAMPLES) m_size=MAXSAMPLES;
  else m_size=size;
}

void HcalSLHCTriggerPrimitiveDigi::setPresamples(int ps) {
  if (ps<0) m_hcalPresamples&=0xFFFFFF0;
  else m_hcalPresamples|=ps&0xF;
}

void HcalSLHCTriggerPrimitiveDigi::setZSInfo(bool unsuppressed, bool markAndPass) {
  if (markAndPass) m_hcalPresamples|=0x10;
  if (unsuppressed) m_hcalPresamples|=0x20;
}

//------------------------------------------------------
// Print out to std::cout, etc
//------------------------------------------------------

std::ostream& operator<<(std::ostream& s, const HcalSLHCTriggerPrimitiveDigi& digi) {
  s << digi.id() << " " << digi.size() << " samples "<< digi.presamples() << " presamples";
  if (digi.zsUnsuppressed()) s << " zsUS";
  if (digi.zsMarkAndPass()) s << " zsM&P";
  s << std::endl;
  for (int i=0; i<digi.size(); i++) 
    s << "  " << digi.sample(i) << std::endl;
  return s;
}
