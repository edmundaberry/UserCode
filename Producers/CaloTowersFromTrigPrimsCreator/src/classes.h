// #include "DataFormats/Common/interface/EDCollection.h"
// #include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowerNonSortedCollection.h"

namespace {

  struct dictionary {
    edm::EDCollection<CaloTower> caloTowerNonSortedCollection_;
    CaloTowerNonSortedCollection theCaloTowerNonSortedCollection_;
    edm::Wrapper<CaloTowerNonSortedCollection> anotherCaloTowerNonSortedCollection;
    edm::Wrapper< edm::EDCollection < CaloTower > > finalCaloTowerNonSortedCollection;
  };

}
