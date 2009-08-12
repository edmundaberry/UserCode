#ifndef PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLCHTRIGGERPRIMITIVEDIGI_H
#define PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLCHTRIGGERPRIMITIVEDIGI_H

#include <ostream>
#include <vector>
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveSample.h"

class HcalSLHCTriggerPrimitiveDigi {

 public:
  typedef HcalTrigTowerDetId key_type;

  explicit HcalSLHCTriggerPrimitiveDigi(const HcalTrigTowerDetId& id);
  
  HcalSLHCTriggerPrimitiveDigi(); 
  ~HcalSLHCTriggerPrimitiveDigi();

  const HcalTrigTowerDetId& id() const { return m_id; }

 private:
  
  HcalTrigTowerDetId m_id;
  int m_size;
  int m_hcalPresamples;

};

#endif
