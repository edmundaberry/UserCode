#ifndef HcalSimAlgos_HcalTriggerPrimitiveAlgo_h
#define HcalSimAlgos_HcalTriggerPrimitiveAlgo_h

#include <map>
#include <vector>
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"

#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTrigPrimDigiCollection.h"

class CaloGeometry;
class IntegerCaloSamples;

class HcalSLHCTriggerPrimitiveAlgo {
public:
  
  //------------------------------------------------------
  // Constructor/destructor
  //------------------------------------------------------
  
  HcalSLHCTriggerPrimitiveAlgo(bool pf, const std::vector<double>& w, 
			       int latency, 
			       int firstTPSample, int TPSize,
			       int minIsoDepth, int maxIsoDepth, bool excludeDepth5 );
  ~HcalSLHCTriggerPrimitiveAlgo();

  //------------------------------------------------------
  // Event-by-event algorithm run function
  //------------------------------------------------------

  void run(const HcalTPGCoder * incoder,
	   const HcalTPGCompressor * outcoder,
	   const HBHEDigiCollection & hbheDigis,
           const HFDigiCollection & hfDigis,
           HcalSLHCTrigPrimDigiCollection & result);
 private:

  //------------------------------------------------------
  // Mapping typedefs
  //------------------------------------------------------

  typedef std::pair <IntegerCaloSamples, IntegerCaloSamples> IntegerCaloSamplesPair;
  typedef std::map<HcalTrigTowerDetId, IntegerCaloSamplesPair> SumMap;

  //------------------------------------------------------
  /// These convert from HBHEDataFrame & HFDataFrame to 
  //  linear IntegerCaloSamples add the linear signal to the map
  //------------------------------------------------------
  
  void addSignal  (const HBHEDataFrame      & frame  );
  void addSignal  (const HFDataFrame        & frame  );
  void addSignal  (const IntegerCaloSamples & samples, int depth );

  //------------------------------------------------------
  /// These convert from linear IntegerCaloSamples to
  // compressed trigger primitive digis
  //------------------------------------------------------
  
  void analyze  (IntegerCaloSamplesPair & samplesPair, HcalSLHCTriggerPrimitiveDigi & result);
  void analyzeHF(IntegerCaloSamplesPair & samplesPair, HcalSLHCTriggerPrimitiveDigi & result);  
  
  //------------------------------------------------------
  // Conversion methods for HB & HE (for analyze function)
  //------------------------------------------------------

  // Weighted sum
  bool doSampleSum      (const IntegerCaloSamples& inputSamples,
			 IntegerCaloSamples& summedSamples,
			 int outlength);
  
  // Collapse/peakfinding
  void doSampleCollapse (const IntegerCaloSamples& originalSamples,
			 const IntegerCaloSamples& summedSamples,
			 IntegerCaloSamples& collapsedSamples,
			 int newprelength,
			 bool SOI_pegged);
  
  // Compression using LUT's
  void doSampleCompress (const IntegerCaloSamples& etSamples,
			 const IntegerCaloSamples& isoSamples,
			 const std::vector<int>& fineGrainSamples,
			 HcalSLHCTriggerPrimitiveDigi & digi);

  //------------------------------------------------------
  // Trig tower geometry
  //------------------------------------------------------

  HcalTrigTowerGeometry theTrigTowerGeometry;

  //------------------------------------------------------
  // Coders
  //------------------------------------------------------

  const HcalTPGCoder * incoder_;
  const HcalTPGCompressor * outcoder_;

  //------------------------------------------------------
  // Energy sum mapping
  //------------------------------------------------------

  SumMap theSumMap;  
  
  //------------------------------------------------------
  // Python file input
  //------------------------------------------------------

  double theThreshold;
  bool peakfind_;
  std::vector<double> weights_;
  int latency_;
  int firstTPSample_;
  int TPSize_;  
  int minIsoDepth_;
  int maxIsoDepth_;
  
  bool excludeDepth5_;
};


#endif



























