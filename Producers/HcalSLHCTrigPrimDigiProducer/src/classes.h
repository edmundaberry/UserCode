#include <vector>
#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiCollection.h"

namespace {

  struct dictionary {
    
    std::vector<HcalSLHCTriggerPrimitiveSample> vTPS_;
    HcalSLHCTrigPrimDigiCollection theHTP_; 
    std::vector<HcalSLHCTriggerPrimitiveDigi> theVecHTP_;
    edm::SortedCollection<HcalSLHCTriggerPrimitiveDigi> vHTP_;
    edm::Wrapper<edm::SortedCollection<HcalSLHCTriggerPrimitiveDigi> > anotherHTP_;
    edm::Wrapper<HcalSLHCTrigPrimDigiCollection> theHTPw_; 

  };

}
