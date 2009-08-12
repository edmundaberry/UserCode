#ifndef PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLCHTRIGGERPRIMITIVEDIGI_H
#define PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLCHTRIGGERPRIMITIVEDIGI_H

#include <ostream>
#include <vector>
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTriggerPrimitiveSample.h"

class HcalSLHCTriggerPrimitiveDigi {

 public:
  
  //------------------------------------------------------
  // Needed for consistency with SortedCollection
  //------------------------------------------------------
  
  typedef HcalTrigTowerDetId key_type;
  explicit HcalSLHCTriggerPrimitiveDigi(const HcalTrigTowerDetId& id);
  const HcalTrigTowerDetId& id() const { return m_id; }

  //------------------------------------------------------
  // Constructor/Destructor
  //------------------------------------------------------
  
  HcalSLHCTriggerPrimitiveDigi(); 
  ~HcalSLHCTriggerPrimitiveDigi();

  //------------------------------------------------------
  // Set information.  MAXSAMPLES sets size limit.
  //------------------------------------------------------

  void setSize      ( int  size );
  void setPresamples( int  presamples );
  void setZSInfo    ( bool unsuppressed, bool markAndPass);
  void setSample    ( int i, const HcalSLHCTriggerPrimitiveSample& sample ) { m_data[i] = sample; }

  static const int MAXSAMPLES = 10;

  //------------------------------------------------------
  // Get the number of samples / presamples
  //------------------------------------------------------
  
  int size      () const { return (m_size           & 0xF); } // 
  int presamples() const { return (m_hcalPresamples & 0xF); } 
  
  //------------------------------------------------------
  // Get nformation about the ZS
  //------------------------------------------------------

  bool zsMarkAndPass () const { return m_hcalPresamples & 0x10; }
  bool zsUnsuppressed() const { return m_hcalPresamples & 0x20; }

  //------------------------------------------------------
  // Get information about individual samples
  //------------------------------------------------------
  
  // Access all stored samples

  const HcalSLHCTriggerPrimitiveSample& operator[](int i) const { return m_data[i]; }
  const HcalSLHCTriggerPrimitiveSample& sample    (int i) const { return m_data[i]; }
  
  // Access "sample of interest" directly
  
  const HcalSLHCTriggerPrimitiveSample& t0() const { return m_data[presamples()]; }
  int SOI_fineGrain      () const { return t0().fineGrain      (); }
  int SOI_compressedEt   () const { return t0().compressedEt   (); }
  int SOI_compressedIsoEt() const { return t0().compressedIsoEt(); }
  
 private:
  
  HcalTrigTowerDetId m_id;
  int m_size;
  int m_hcalPresamples;
  HcalSLHCTriggerPrimitiveSample m_data [MAXSAMPLES];


};

#endif
