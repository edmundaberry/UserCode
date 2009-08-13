#ifndef PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLHCTRIGGERPRIMITIVESAMPLE_H
#define PRODUCERS_HCALSLHCTRIGPRIMDIGIPRODUCER_HCALSLHCTRIGGERPRIMITIVESAMPLE_H

#include <boost/cstdint.hpp>
#include <ostream>

class HcalSLHCTriggerPrimitiveSample {
  
 public:  
  
  HcalSLHCTriggerPrimitiveSample();
  HcalSLHCTriggerPrimitiveSample(uint32_t data);
  HcalSLHCTriggerPrimitiveSample(int encodedEt, int encodedIsoEt, int fineGrain, int slb, int slbchan);
  
  uint32_t raw() const {return theSample; }  
  int slb            () const { return ((theSample >> 29) & 0x7 ); } // slb          has 3 bits, (30-32)
  int slbChan        () const { return ((theSample >> 27) & 0x3 ); } // slbchan      has 2 bits, (28-29)
  int slbAndChan     () const { return ((theSample >> 27) & 0x1F); }
  int fineGrain      () const { return ((theSample >> 16) & 0xFF); } // fine grain   has 8 bits, (17-24)
  int compressedIsoEt() const { return ((theSample >> 8 ) & 0xFF); } // isolation et has 8 bits, (9 -16)
  int compressedEt   () const { return ((theSample >> 0 ) & 0xFF); } // et           has 8 bits, (1 -8 )
  
 private:
  uint32_t theSample;

};

std::ostream& operator<<(std::ostream& s, const HcalSLHCTriggerPrimitiveSample& samp);

#endif
