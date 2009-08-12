#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTriggerPrimitiveDigi.h"

HcalSLHCTriggerPrimitiveDigi::HcalSLHCTriggerPrimitiveDigi(){}
HcalSLHCTriggerPrimitiveDigi::~HcalSLHCTriggerPrimitiveDigi(){}
HcalSLHCTriggerPrimitiveDigi::HcalSLHCTriggerPrimitiveDigi(const HcalTrigTowerDetId& id):
  m_id  (id),
  m_size(0),
  m_hcalPresamples(0)
{}
