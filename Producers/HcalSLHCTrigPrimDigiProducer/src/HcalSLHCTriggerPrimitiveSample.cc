#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTriggerPrimitiveSample.h"

//------------------------------------------------------
// Null constructor
//------------------------------------------------------

HcalSLHCTriggerPrimitiveSample::HcalSLHCTriggerPrimitiveSample(): theSample(0){}

//------------------------------------------------------
// Constructor where we already know the data word
//------------------------------------------------------

HcalSLHCTriggerPrimitiveSample::HcalSLHCTriggerPrimitiveSample(uint32_t data): theSample(data){}

//------------------------------------------------------
// Constructor where we make the data word ourselves
//------------------------------------------------------

HcalSLHCTriggerPrimitiveSample::HcalSLHCTriggerPrimitiveSample(int encodedEt, int encodedIsoEt, int encodedFG, int slb, int slbchan){
  
  // bits are numbered starting from 1 and going to 32, inclusive.
  // there are 3 unused bits: from 25-27, inclusive.
  
  theSample = 
    (((slb         )&0x7 )<<29) | // slb          has 3 bits, (30-32)
    (((slbchan     )&0x3 )<<27) | // slbchan      has 2 bits, (28-29)
    (((encodedFG   )&0xFF)<<16) | // fine grain   has 8 bits, (17-24)
    (((encodedIsoEt)&0xFF)<<8 ) | // isolation et has 8 bits, (9 -16)
    (((encodedEt   )&0xFF)<<0 );  // et           has 8 bits, (1 -8 )
  
}

//------------------------------------------------------
// Print function for std::cout, etc
//------------------------------------------------------

std::ostream& operator<<(std::ostream& s, const HcalSLHCTriggerPrimitiveSample& samp) {
  return s << "Et="      << samp.compressedEt()    << ", " 
	   << "IsoEt = " << samp.compressedIsoEt() << ", " 
	   << "FG="      << samp.fineGrain() ;
}

