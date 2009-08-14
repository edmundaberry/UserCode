#include "Producers/HcalSLHCTrigPrimDigiProducer/interface/HcalSLHCTriggerPrimitiveAlgo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
using namespace std;

HcalSLHCTriggerPrimitiveAlgo::HcalSLHCTriggerPrimitiveAlgo(bool pf, const std::vector<double>& w,
							   int latency, uint32_t FG_threshold, 
							   uint32_t ZS_threshold, int firstTPSample, int TPSize,
							   int minIsoDepth, int maxIsoDepth ):
  incoder_(0), 
  outcoder_(0), 
  theThreshold(0),
  peakfind_(pf), 
  weights_(w), 
  latency_(latency), 
  FG_threshold_(FG_threshold), 
  ZS_threshold_(ZS_threshold), 
  firstTPSample_(firstTPSample), 
  TPSize_(TPSize),
  minIsoDepth_(minIsoDepth),
  maxIsoDepth_(maxIsoDepth)
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
  theFGSumMap.clear();    
  theTowerMapFG.clear();  

  //------------------------------------------------------
  // Loop over the digi collections and fill energy maps
  //------------------------------------------------------

  HBHEDigiCollection::const_iterator hbheItr     = hbheDigis.begin();
  HBHEDigiCollection::const_iterator hbheItr_end = hbheDigis.end();

  HFDigiCollection::const_iterator hfItr     = hfDigis.begin();
  HFDigiCollection::const_iterator hfItr_end = hfDigis.end();
  
  // this is ok
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
// Add signals from the HBHE digis to the mapping
// 
// This should be the same as the pre-upgrade method.
//------------------------------------------------------

void HcalSLHCTriggerPrimitiveAlgo::addSignal(const HBHEDataFrame & frame) {

  int depth = frame.id().depth();
  
  //Hack for 300_pre10, should be removed.
  if (frame.id().depth()==5) return;
  
  //------------------------------------------------------
  // Get the trigger tower id(s) for this digi
  //------------------------------------------------------
  
  std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry.towerIds(frame.id());
  assert(ids.size() == 1 || ids.size() == 2);
  
  IntegerCaloSamples samples1(ids[0], int(frame.size()));

  samples1.setPresamples(frame.presamples());
  incoder_->adc2Linear(frame, samples1);
  
  //------------------------------------------------------
  // If there are two id's make a second trigprim for 
  // the other one, and split the energy
  //------------------------------------------------------

  if(ids.size() == 2) {
    
    IntegerCaloSamples samples2(ids[1], samples1.size());
    for(int i = 0; i < samples1.size(); ++i) 
      {
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
        
    uint32_t fgid;

    // Mask off depths: fgid is the same for both depths
    
    fgid = (frame.id().rawId() | 0x1c000) ;
	 
    SumMapFG::iterator itr = theFGSumMap.find(fgid);
   
    if(itr == theFGSumMap.end()) {
      theFGSumMap.insert(std::make_pair(fgid, samples));
    } 
    else {
      // wish CaloSamples had a +=
      for(int i = 0; i < samples.size(); ++i) {
	(itr->second)[i] += samples[i];        
      }      
    }

    // Depth =2 is the second entry in map (sum). Use its original Hcal Det Id to obtain trigger tower
    if (frame.id().depth()==2)
      {      
      for(unsigned int n = 0; n < ids.size(); n++)
      	  {
          theTowerMapFG.insert(TowerMapFG::value_type(ids[n],itr->second));
	  }
      }
    
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
      theSumMap.insert(std::make_pair(id, make_pair(samples, samples)));
    else {

      IntegerCaloSamples emptySamples(id, samples.size());
      emptySamples.setPresamples (samples.presamples());

      theSumMap.insert(std::make_pair(id, make_pair(samples, emptySamples)));

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

void HcalSLHCTriggerPrimitiveAlgo::analyze(IntegerCaloSamplesPair & samples, 
					   HcalSLHCTriggerPrimitiveDigi & result){

  //------------------------------------------------------
  // First find output samples length and new presamples
  //------------------------------------------------------
  
  int shrink       = weights_.size()-1; 
  int outlength    = samples.first.size() - shrink;
  int newprelength = ((samples.first.presamples()+1)-weights_.size())+latency_;
  
  std::vector<bool>  finegrain(TPSize_,false);
  
  //------------------------------------------------------
  // Sum over all samples and make a summed sample
  // 
  // doSampleSum returns a bool that says if the input
  //   IntegerCaloSample is above the threshold value
  //------------------------------------------------------
  
  IntegerCaloSamples allSum(samples.first.id(), outlength);
  IntegerCaloSamples isoSum(samples.first.id(), outlength);
  
  bool allSOI_pegged = doSampleSum (samples.first , allSum, outlength);
  bool isoSOI_pegged = doSampleSum (samples.second, isoSum, outlength);
  
  //------------------------------------------------------
  // Collapse the sample
  //------------------------------------------------------
  
  IntegerCaloSamples allCollapsed (samples.first.id(), TPSize_  );  
  IntegerCaloSamples isoCollapsed (samples.first.id(), TPSize_  );
  
  doSampleCollapse (samples.first , allSum, allCollapsed, newprelength, allSOI_pegged);
  doSampleCollapse (samples.second, isoSum, isoCollapsed, newprelength, isoSOI_pegged);
  
  //------------------------------------------------------
  // Compress into an HcalSLHCTriggerPrimitiveDigi
  // FineGrain is left blank intentionally -- this should 
  // be filled in later.
  //------------------------------------------------------
  
  std::vector<int> nullFineGrain ( allCollapsed.size(), 0 );

  compress ( allCollapsed, isoCollapsed, nullFineGrain, result );
  
}


void HcalSLHCTriggerPrimitiveAlgo::analyzeHF(IntegerCaloSamplesPair & samples, 
					     HcalSLHCTriggerPrimitiveDigi & result){
  
  /*
  std::vector<bool> finegrain(TPSize_,false);
  
  HcalTrigTowerDetId detId_(samples.id()); 
   
  // get information from Tower map
  for(TowerMapFG::iterator mapItr = theTowerMapFG.begin(); mapItr != theTowerMapFG.end(); ++mapItr){

    HcalTrigTowerDetId detId(mapItr->first);
    if (detId == detId_) {
   
      for (int i=firstTPSample_; i < firstTPSample_+TPSize_; ++i) {
	bool set_fg = false;
	mapItr->second[i] >= FG_threshold_ ? set_fg = true : false;
	finegrain[i - firstTPSample_] = (finegrain[i - firstTPSample_] || set_fg);
      }
    }
  }  

  */

  IntegerCaloSamples output(samples.first.id(),TPSize_);

  output.setPresamples(samples.first.presamples() - firstTPSample_);

  for(int ibin2 = 0; ibin2 < TPSize_; ++ibin2) {    
    output[ibin2]=samples.first[ibin2+firstTPSample_]/4;
    if (output[ibin2] > 0x3FF) output[ibin2] = 0x3FF;  //Compression is 1 to 1 with saturation at 8 bits
    
  }

  // Get null isolation sample and null fine grain
  IntegerCaloSamples nullIsoSamples( samples.first.id(), output.size() );
  nullIsoSamples.setPresamples(samples.first.presamples());

  std::vector<int>   nullFineGrain ( samples.first.size(), 0 );
  
  compress ( output, nullIsoSamples, nullFineGrain, result );  

}


void HcalSLHCTriggerPrimitiveAlgo::doSampleCollapse (const IntegerCaloSamples& originalSamples,
						     const IntegerCaloSamples& summedSamples,
						     IntegerCaloSamples& collapsedSamples,
						     int newprelength, bool SOI_pegged ){
  
  if(peakfind_) {    
    
    collapsedSamples.setPresamples(newprelength - firstTPSample_);
    
    for(int ibin2 = 0; ibin2 < TPSize_; ++ibin2) {
      
      int idx = firstTPSample_ + ibin2;
      
      //if peak found
      if ( originalSamples[idx] >  originalSamples[idx-1] && 
	   originalSamples[idx] >= originalSamples[idx+1] && 
	   originalSamples[idx] >  theThreshold) 	 
	collapsedSamples[ibin2]=summedSamples[idx];
      
      //if no peak
      else collapsedSamples[ibin2]=0; 
    }
    
    if(SOI_pegged == true){
      collapsedSamples[collapsedSamples.presamples()] = 0x3FF;
    }        
  }
  
  //No peak finding
  else {
    
    collapsedSamples.setPresamples(newprelength - firstTPSample_ +1);
    
    for(int ibin2 = 0; ibin2 < TPSize_; ++ibin2) 
      collapsedSamples[ibin2]=summedSamples[ibin2+firstTPSample_];//just pass value  
    
  }   
  
}

bool HcalSLHCTriggerPrimitiveAlgo::doSampleSum (const IntegerCaloSamples& inputSamples, 
						IntegerCaloSamples& summedSamples,
						int outlength ){
  
  //Test is SOI input is pegged (saturated) before summing
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

void HcalSLHCTriggerPrimitiveAlgo::compress (const IntegerCaloSamples& etSamples,
					     const IntegerCaloSamples& isoSamples,
					     const vector<int> & fineGrainSamples,
					     HcalSLHCTriggerPrimitiveDigi & digi){
  
  // We have to go through HcalTriggerPrimitiveSample to do this
  //   (CaloTPGTranscoder::hcalCompress returns an HcalTriggerPrimitiveSample)
  // I don't know of any other way to do the compression without editing the LUT storage class

  for (int i = 0; i < etSamples.size(); ++i){
    
    int compressedEt  = outcoder_ -> compress ( etSamples.id() , etSamples [i], false ).compressedEt();
    int compressedIso = outcoder_ -> compress ( isoSamples.id(), isoSamples[i], false ).compressedEt();
    int fineGrain     = fineGrainSamples[i];

    digi.setSample( i, HcalSLHCTriggerPrimitiveSample ( compressedEt, compressedIso, fineGrain, 0, 0 ));

  }

  
  digi.setPresamples ( etSamples.presamples() );

  
}
