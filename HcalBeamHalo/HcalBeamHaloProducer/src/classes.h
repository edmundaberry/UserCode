#ifndef HCALBEAMHALO_HCALBEAMHALOPRODUCER_classes_h
#define HCALBEAMHALO_HCALBEAMHALOPRODUCER_classes_h

#include "HcalBeamHalo/HcalBeamHaloProducer/interface/HcalBeamHaloCollection.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {

  struct dictionary {
    std::vector<HcalDetId> hcalDetId_;
    std::vector<HcalBeamHalo> hcalBeamHaloCollection_;


    HcalDetId theHcalDetId_;
    HcalBeamHaloCollection theHcalBeamHaloCollection_;
    
    edm::Wrapper<HcalBeamHaloCollection> anotherHcalBeamHaloCollection;
    edm::Wrapper<HcalDetId> anotherHcalDetIdCollection;
    
    edm::Wrapper< std::vector < HcalBeamHalo > > finalHcalBeamHaloCollection;
    edm::Wrapper< std::vector < HcalDetId > > finalHcalDetIdCollection;
  };

}

#endif
