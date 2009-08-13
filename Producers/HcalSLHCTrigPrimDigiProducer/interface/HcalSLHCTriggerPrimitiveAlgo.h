#ifndef HcalSimAlgos_HcalTriggerPrimitiveAlgo_h
#define HcalSimAlgos_HcalTriggerPrimitiveAlgo_h

#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"

#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiCollection.h"

#include <map>
#include <vector>
class CaloGeometry;
class IntegerCaloSamples;

class HcalSLHCTriggerPrimitiveAlgo {
public:
  HcalSLHCTriggerPrimitiveAlgo(bool pf, const std::vector<double>& w, 
			       int latency, uint32_t FG_threshold, 
			       uint32_t ZS_threshold, int firstTPSample, int TPSize,
			       int minIsoDepth, int maxIsoDepth );
  ~HcalSLHCTriggerPrimitiveAlgo();

  void run(const HcalTPGCoder * incoder,
	   const HcalTPGCompressor * outcoder,
	   const HBHEDigiCollection & hbheDigis,
           const HFDigiCollection & hfDigis,
           HcalSLHCTrigPrimDigiCollection & result);
 private:

  typedef std::pair <IntegerCaloSamples, IntegerCaloSamples> IntegerCaloSamplesPair;
  
  typedef std::map<HcalTrigTowerDetId, IntegerCaloSamplesPair> SumMap;
  typedef std::map<uint32_t, IntegerCaloSamples> SumMapFG;
  typedef std::multimap<HcalTrigTowerDetId, IntegerCaloSamples> TowerMapFG;

  /// adds the signal to the map
  void addSignal  (const HBHEDataFrame      & frame  );
  void addSignal  (const HFDataFrame        & frame  );
  void addSignal  (const IntegerCaloSamples & samples, int depth );
  void addSignalFG(const IntegerCaloSamples & samples);
  /// adds the actual RecHits
  void analyze  (IntegerCaloSamplesPair & samplesPair, HcalSLHCTriggerPrimitiveDigi & result);
  void analyzeHF(IntegerCaloSamplesPair & samplesPair, HcalSLHCTriggerPrimitiveDigi & result);
  
  
  bool doSampleSum (const IntegerCaloSamples& inputSamples,
		    IntegerCaloSamples& summedSamples,
		    int outlength);

  void doSampleCollapse (const IntegerCaloSamples& originalSamples,
			 const IntegerCaloSamples& summedSamples,
			 IntegerCaloSamples& collapsedSamples,
			 int newprelength,
			 bool SOI_pegged);

  // Performs compression and fills HcalSLHCTriggerPrimitiveDigis
  void compress (const IntegerCaloSamples& etSamples,
		 const IntegerCaloSamples& isoSamples,
		 const std::vector<int>& fineGrainSamples,
		 HcalSLHCTriggerPrimitiveDigi & digi);
  
  std::vector<HcalTrigTowerDetId> towerIds(const HcalDetId & id) const;

  HcalTrigTowerGeometry theTrigTowerGeometry; // from event setup eventually?

  const HcalTPGCoder * incoder_;
  const HcalTPGCompressor * outcoder_;

  SumMap     theSumMap;  
  SumMap     theIsoSumMap;
  SumMapFG   theFGSumMap;
  TowerMapFG theTowerMapFG;

  double theThreshold;
  bool peakfind_;
  std::vector<double> weights_;
  int latency_;
  uint32_t FG_threshold_;
  uint32_t ZS_threshold_;
  int firstTPSample_;
  int TPSize_;
  
  int minIsoDepth_;
  int maxIsoDepth_;
};


#endif



























