#include "Analyzers/HcalDigiAnalyzer/interface/HcalDigiAnalyzer.h"
#include "Analyzers/HcalDigiAnalyzer/interface/DigiTree.h"
#include "Analyzers/HcalDigiAnalyzer/interface/FillDigiTree.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "FWCore/Framework/interface/MakerMacros.h"

typedef HcalDigiAnalyzer<HBHEDigiCollection> HBHEDigiAnalyzer;
typedef HcalDigiAnalyzer<HODigiCollection  > HODigiAnalyzer;
typedef HcalDigiAnalyzer<HFDigiCollection  > HFDigiAnalyzer;

DEFINE_FWK_MODULE (HBHEDigiAnalyzer);
DEFINE_FWK_MODULE (HFDigiAnalyzer);
DEFINE_FWK_MODULE (HODigiAnalyzer);
