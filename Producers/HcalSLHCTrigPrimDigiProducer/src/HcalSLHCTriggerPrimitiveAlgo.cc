#include <iostream>
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTriggerPrimitiveAlgo.h"

HcalSLHCTriggerPrimitiveAlgo::HcalSLHCTriggerPrimitiveAlgo(bool pf, const std::vector<double>& w,
							   int latency, 
							   int firstTPSample, int TPSize,
							   int minIsoDepth, int maxIsoDepth, bool excludeDepth5, float lsb , bool slhcMode ):
  incoder_(0), 
  outcoder_(0), 
  theThreshold(0),
  peakfind_(pf), 
  weights_(w), 
  latency_(latency), 
  firstTPSample_(firstTPSample), 
  TPSize_(TPSize),
  minIsoDepth_(minIsoDepth),
  maxIsoDepth_(maxIsoDepth),
  excludeDepth5_(excludeDepth5),
  slhcMode_ ( slhcMode ),
  compressionLsb_( lsb )
{}

HcalSLHCTriggerPrimitiveAlgo::~HcalSLHCTriggerPrimitiveAlgo(){}

void HcalSLHCTriggerPrimitiveAlgo::run(const HcalTPGCoder * incoder,
				       const HcalTPGCompressor * outcoder,
				       const HBHEDigiCollection & hbheDigis,
				       const HFDigiCollection & hfDigis,
				       HcalSLHCTrigPrimDigiCollection & result){
  
  //------------------------------------------------------
  // Set coders and clear energy maps
  //------------------------------------------------------
  
  incoder_=incoder;
  outcoder_=outcoder;

  theSumMap.clear();      

  //------------------------------------------------------
  // Loop over the digi collections and fill energy maps
  //------------------------------------------------------

  HBHEDigiCollection::const_iterator hbheItr     = hbheDigis.begin();
  HBHEDigiCollection::const_iterator hbheItr_end = hbheDigis.end();

  HFDigiCollection::const_iterator hfItr     = hfDigis.begin();
  HFDigiCollection::const_iterator hfItr_end = hfDigis.end();
  
  for(; hbheItr != hbheItr_end; ++hbheItr)
    addSignal(*hbheItr);

  for(; hfItr != hfItr_end; ++hfItr)
    addSignal(*hfItr);

  //------------------------------------------------------
  // Once energy maps are filled,
  //   1 - Loop over entries in the energy sum map
  //   2 - Convert summed IntegerCaloSamples to TP's
  //------------------------------------------------------
  
  for(SumMap::iterator mapItr = theSumMap.begin(); mapItr != theSumMap.end(); ++mapItr){
      
    // Push back an empty TP digi that contains only the HcalTrigTowerDetId
    result.push_back(HcalSLHCTriggerPrimitiveDigi(mapItr->first));
    
    HcalTrigTowerDetId detId((mapItr->second).first.id());
    
    // Fill in the rest of the information for the TP digi
    if(detId.ietaAbs() >= theTrigTowerGeometry.firstHFTower())
      analyzeHF(mapItr->second, result.back());
    else
      analyze  (mapItr->second, result.back());
  }
  
  //------------------------------------------------------
  // Free up memory and return
  //------------------------------------------------------

  theSumMap.clear();  

  return;
}

//------------------------------------------------------
// LUT function
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveAlgo::adc2Linear(const HBHEDataFrame& frame, IntegerCaloSamples & sample ){
  
  int first_sample = 0;
  int last_sample  = frame.size() - 1;
  
  for (int iSample = first_sample; iSample <= last_sample; ++iSample){
    int value = (int) (frame[iSample].adc() * compressionLsb_);
    value &= 0x3FF;

    sample[iSample] = value;
  }
  
}

//------------------------------------------------------
// Add signals from the HBHE digis to the mapping
// 
// This should be the same as the pre-upgrade method.
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveAlgo::addSignal(const HBHEDataFrame & frame) {

  int depth = frame.id().depth();
  
  //------------------------------------------------------
  // "Hack for 300_pre10" from the original code.
  // User can turn this off in the python cfg file.
  //------------------------------------------------------
  
  if (excludeDepth5_ && frame.id().depth()==5) return;
  
  //------------------------------------------------------
  // Get the trigger tower id(s) for this digi
  //------------------------------------------------------
  
  std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry.towerIds(frame.id());
  assert(ids.size() == 1 || ids.size() == 2);
  
  IntegerCaloSamples samples1(ids[0], int(frame.size()));

  samples1.setPresamples(frame.presamples());
  
  //------------------------------------------------------
  // Compress HBHEDataFrame ADC -> IntegerCaloSamples
  // 
  //------------------------------------------------------

  adc2Linear ( frame, samples1 ) ;
  //incoder_->adc2Linear(frame, samples1);  
  
  //------------------------------------------------------
  // If there are two id's make a second trigprim for 
  // the other one, and split the energy
  //------------------------------------------------------

  if(ids.size() == 2) {
    
    IntegerCaloSamples samples2(ids[1], samples1.size());

    for(int i = 0; i < samples1.size(); ++i) {
      samples1[i] = uint32_t(samples1[i]*0.5);
      samples2[i] = samples1[i];
    }
    
    samples2.setPresamples(frame.presamples());
    addSignal(samples2, depth);
  }
  
  addSignal(samples1, depth);

}


void HcalSLHCTriggerPrimitiveAlgo::addSignal(const HFDataFrame & frame) {
  
  int depth = frame.id().depth();

  if( depth == 1 || depth == 2) {

    std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry.towerIds(frame.id());
    
    assert(ids.size() == 1);
    
    IntegerCaloSamples samples(ids[0], frame.size());
    
    samples.setPresamples(frame.presamples());

    incoder_->adc2Linear(frame, samples);
    
    addSignal(samples, depth);

  }
}
  
//------------------------------------------------------
// Add this IntegerCaloSamples to SumMap
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveAlgo::addSignal(const IntegerCaloSamples & samples, int depth ){

  //------------------------------------------------------
  // Is this from an isolated digi?
  //------------------------------------------------------

  bool hasIsoDepth = (depth >= minIsoDepth_ && depth <= maxIsoDepth_ );

  //------------------------------------------------------
  // If this ID isn't in the map already, add it
  //
  // Recall that SumMap maps from an HcalTrigTowerDetId
  //  to a pair of IntegerCaloSamples:
  //   First  IntegerCaloSamples is for everything
  //   Second IntegerCaloSamples is for isolation depths
  //------------------------------------------------------

  HcalTrigTowerDetId id(samples.id());
  SumMap::iterator itr = theSumMap.find(id);
  if(itr == theSumMap.end()) {
    if (hasIsoDepth) 
      theSumMap.insert(std::make_pair(id, std::make_pair(samples, samples)));
    else {

      IntegerCaloSamples emptySamples(id, samples.size());
      emptySamples.setPresamples (samples.presamples());

      theSumMap.insert(std::make_pair(id, std::make_pair(samples, emptySamples)));

    }
  } 

  //------------------------------------------------------
  // If this ID is in the map, add to the existing entry
  //------------------------------------------------------
  
  else {
    for(int i = 0; i < samples.size(); ++i)  {
      ((itr->second).first)[i] += samples[i];   
      if (hasIsoDepth) ((itr->second).second)[i] += samples[i];   
    }
  }
}

//------------------------------------------------------
// Go from an IntegerCaloSamples to a TP digi (for HBHE)
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveAlgo::analyze(IntegerCaloSamplesPair & samples, 
					   HcalSLHCTriggerPrimitiveDigi & result){

  //------------------------------------------------------
  // First find output samples length and new presamples
  //------------------------------------------------------
  
  int shrink       = weights_.size()-1; 
  int outlength    = samples.first.size() - shrink;
  int newprelength = ((samples.first.presamples()+1)-weights_.size())+latency_;
  
  //------------------------------------------------------
  // Declare all IntegerCaloSamples required
  // NOTE: this is a really inefficient way of doing this.  
  // Must be fixed.
  //------------------------------------------------------

  IntegerCaloSamples allSum(samples.first.id(), outlength);
  IntegerCaloSamples isoSum(samples.first.id(), outlength);

  IntegerCaloSamples allCollapsed (samples.first.id(), TPSize_  );  
  IntegerCaloSamples isoCollapsed (samples.first.id(), TPSize_  );

  std::vector<int> nullFineGrain ( allCollapsed.size(), 0 );

  //------------------------------------------------------
  // Sum over all samples and make a summed sample
  // 
  // doSampleSum returns a bool that says if the input
  //   IntegerCaloSample is saturated
  //------------------------------------------------------
  
  bool allSOI_pegged = doSampleSum (samples.first , allSum, outlength);
  bool isoSOI_pegged = doSampleSum (samples.second, isoSum, outlength);
  
  //------------------------------------------------------
  // Collapse the sample.  Peakfinder is here.
  //------------------------------------------------------
    
  doSampleCollapse (samples.first , allSum, allCollapsed, newprelength, allSOI_pegged);
  doSampleCollapse (samples.second, isoSum, isoCollapsed, newprelength, isoSOI_pegged);
  
  //------------------------------------------------------
  // Compress into an HcalSLHCTriggerPrimitiveDigi
  // FineGrain is left blank intentionally -- this should 
  // be filled in later.
  //------------------------------------------------------
  
  doSampleCompress ( allCollapsed, isoCollapsed, nullFineGrain, result );
  
}

//------------------------------------------------------
// Go from an IntegerCaloSamples to a TP digi (for HF)
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveAlgo::analyzeHF(IntegerCaloSamplesPair & samples, 
					     HcalSLHCTriggerPrimitiveDigi & result){
  
  //------------------------------------------------------
  // Declare needed IntegerCaloSamples 
  //------------------------------------------------------
  
  IntegerCaloSamples output(samples.first.id(),TPSize_);
  output.setPresamples(samples.first.presamples() - firstTPSample_);

  //------------------------------------------------------
  // "Collapse" algorithm for HF is just a scaling
  // Saturation is at 8 bits    
  //------------------------------------------------------
  
  for(int ibin2 = 0; ibin2 < TPSize_; ++ibin2) {    
    output[ibin2]=samples.first[ibin2+firstTPSample_]/4;  // Scaling
    if (output[ibin2] > 0x3FF) output[ibin2] = 0x3FF;     // Saturation
  }

  //------------------------------------------------------
  // Get null isolation sample and null fine grain
  //------------------------------------------------------

  std::vector<int>   nullFineGrain ( samples.first.size(), 0 );
  IntegerCaloSamples nullIsoSamples( samples.first.id(), output.size() );
  nullIsoSamples.setPresamples(samples.first.presamples());
  
  doSampleCompress ( output, nullIsoSamples, nullFineGrain, result );  

}

//------------------------------------------------------
// "Collapse" method for HBHE.  Peakfinder is here.
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveAlgo::doSampleCollapse (const IntegerCaloSamples& originalSamples,
						     const IntegerCaloSamples& summedSamples,
						     IntegerCaloSamples& collapsedSamples,
						     int newprelength, bool SOI_pegged ){
  
  //------------------------------------------------------
  // ...with peakfinding
  //------------------------------------------------------

  if(peakfind_) {    
    
    collapsedSamples.setPresamples(newprelength - firstTPSample_);
    
    for(int ibin2 = 0; ibin2 < TPSize_; ++ibin2) {
      
      int idx = firstTPSample_ + ibin2;
      
      // peak found
      if ( originalSamples[idx] >  originalSamples[idx-1] && 
	   originalSamples[idx] >= originalSamples[idx+1] && 
	   originalSamples[idx] >  theThreshold) 	 
	collapsedSamples[ibin2]=summedSamples[idx];
      
      // no peak found
      else collapsedSamples[ibin2]=0; 
    }
    
    // saturation
    if(SOI_pegged == true){
      collapsedSamples[collapsedSamples.presamples()] = 0x3FF;
    }        
  }
  
  //------------------------------------------------------
  // ... without peakfinding (just pass everything)
  //------------------------------------------------------
  
  else {
    
    collapsedSamples.setPresamples(newprelength - firstTPSample_ +1);
    
    for(int ibin2 = 0; ibin2 < TPSize_; ++ibin2) 
      collapsedSamples[ibin2]=summedSamples[ibin2+firstTPSample_];
    
  }   
  
}

//------------------------------------------------------
// Weighted sum method
//------------------------------------------------------

bool HcalSLHCTriggerPrimitiveAlgo::doSampleSum (const IntegerCaloSamples& inputSamples, 
						IntegerCaloSamples& summedSamples,
						int outlength ){
  
  
  bool SOI_pegged = (inputSamples[inputSamples.presamples()] > 0x3FF);
  
  //slide algo window
  for(int ibin = 0; ibin < outlength ; ++ibin) {
    
    int algosumvalue = 0;
    
    //add up value * scale factor      
    for(unsigned int i = 0; i < weights_.size(); i++) 
      algosumvalue += int(inputSamples[ibin+i] * weights_[i]);
    
    if      (algosumvalue < 0    ) summedSamples[ibin] = 0;           // low-side
    else if (algosumvalue > 0x3FF) summedSamples[ibin] = 0x3FF;       // high-side
    else                           summedSamples[ibin] = algosumvalue;// assign value to sum[]
  }
  
  return SOI_pegged;

}

//------------------------------------------------------
// Compression method (for HBHE and HF)
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveAlgo::doSampleCompress (const IntegerCaloSamples& etSamples,
						     const IntegerCaloSamples& isoSamples,
						     const std::vector<int> & fineGrainSamples,
						     HcalSLHCTriggerPrimitiveDigi & digi){
  
  //------------------------------------------------------
  // Looks like we have to go through 
  //   HcalTriggerPrimitiveSample in order to compress.
  //------------------------------------------------------

  for (int i = 0; i < etSamples.size(); ++i){
    
    int compressedEt  = outcoder_ -> compress ( etSamples.id() , etSamples [i], false ).compressedEt();
    int compressedIso = outcoder_ -> compress ( isoSamples.id(), isoSamples[i], false ).compressedEt();
    int fineGrain     = fineGrainSamples[i];

    digi.setSample( i, HcalSLHCTriggerPrimitiveSample ( compressedEt, compressedIso, fineGrain, 0, 0 ));

  }

  
  digi.setPresamples ( etSamples.presamples() );

  
}
