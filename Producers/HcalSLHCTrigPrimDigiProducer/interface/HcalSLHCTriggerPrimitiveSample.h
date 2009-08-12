#ifndef PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLHCTRIGGERPRIMITIVESAMPLE_H
#define PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLHCTRIGGERPRIMITIVESAMPLE_H

#include <boost/cstdint.hpp>
#include <ostream>

class HcalSLHCTriggerPrimitiveSample {
  
 public:  
  
  HcalSLHCTriggerPrimitiveSample();
  HcalSLHCTriggerPrimitiveSample(uint32_t data);
  HcalSLHCTriggerPrimitiveSample(int encodedEt, int encodedIsoEt, int fineGrain, 
				 int slb, int slbchan);
  
  uint32_t raw() const {return theSample; }
  
  int compressedEt   () const { return ((theSample >> 0 ) & 0xFF); }
  int compressedIsoEt() const { return ((theSample >> 8 ) & 0xFF); }
  int fineGrain      () const { return ((theSample >> 16) & 0xFF); }
  int slbChan        () const { return ((theSample >> 27) & 0x3 ); }
  int slb            () const { return ((theSample >> 29) & 0x7 ); }
  int slbAndChan     () const { return ((theSample >> 27) & 0x1F); }
  
 private:
  uint32_t theSample;

};

std::ostream& operator<<(std::ostream& s, const HcalSLHCTriggerPrimitiveSample& samp);

#endif
