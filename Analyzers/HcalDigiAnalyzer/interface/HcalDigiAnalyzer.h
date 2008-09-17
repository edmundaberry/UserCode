#ifndef Analyzers_HcalDigiAnalyzer_H
#define Analyzers_HcalDigiAnalyzer_H

// system include files
#include <vector>
#include <utility>
#include <ostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/HcalObjects/interface/HcalQIECoder.h"

#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HFDataFrame.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"

#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"
#include "CalibCalorimetry/HcalTPGAlgos/interface/HcaluLUTTPGCoder.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseContainmentCorrection.h"

#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/HcalTowerAlgo/src/HcalHardcodeGeometryData.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"

#include "Analyzers/HcalDigiAnalyzer/interface/DigiTree.h"
#include "Analyzers/HcalDigiAnalyzer/interface/FillDigiTree.h"

#include "DataFormats/Candidate/interface/const_iterator.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TROOT.h"

#include <iostream>
#include <string>

//-----------------------------------------------------
// This is a purposely hardcoded value.
// Do not change it.

static double MaximumFractionalError_ =  0.0005; 

//-----------------------------------------------------

template< typename T>
class HcalDigiAnalyzer : public edm::EDAnalyzer {
 public:
  explicit HcalDigiAnalyzer(const edm::ParameterSet&);
  ~HcalDigiAnalyzer();
    
 private:

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  void printTrigPrimError(HcalTrigTowerDetId trigTowerDetId, HcalDetId cell,
			  int nTrigTowerDetIds, int ntrigprim);
  

  //-----------------------------------------------------
  // Root tree objects
  //-----------------------------------------------------

  DigiTree     m_digiTree;
  FillDigiTree m_fillDigi;
    
  //-----------------------------------------------------
  // Out file info
  //-----------------------------------------------------

  std::string outPath_;
  std::string outSuffix_;
  std::string subdetName_;
  std::string rootFile_;

  //-----------------------------------------------------
  // edm::InputTags 
  //-----------------------------------------------------

  edm::InputTag hcalTrigPrimTag_;
  edm::InputTag hcalDigiTag_;

  //-----------------------------------------------------
  // Output flags
  //-----------------------------------------------------

  bool doRHPulseCorrect_;
  bool scanForSpikes_;
  bool verbose_;  
  int subdet_;

  //-----------------------------------------------------
  // Pulse corrections 
  // (for calculating RH energy on the fly)
  //-----------------------------------------------------
  
  std::auto_ptr<HcalPulseContainmentCorrection> pulseCorrectionPtr;
  
  int firstRHSample_;
  int rhSamplesToAdd_;
  float phaseNS_;
  
};

template < typename T >
HcalDigiAnalyzer<T>::HcalDigiAnalyzer(const edm::ParameterSet& iConfig)
{

  //-----------------------------------------------------
  // Declare namespace
  //-----------------------------------------------------
  
  using namespace std;

  //-----------------------------------------------------
  // Get input from .cfg file
  //-----------------------------------------------------

  verbose_          = iConfig.getUntrackedParameter<bool>         ("verbose",false);
  doRHPulseCorrect_ = iConfig.getUntrackedParameter<bool>         ("doRecHitPulseCorrection",true);
  scanForSpikes_    = iConfig.getUntrackedParameter<bool>         ("scanForSpikes",false);
  hcalTrigPrimTag_  = iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigPrimTag");
  hcalDigiTag_      = iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag"); 
  outPath_          = iConfig.getUntrackedParameter<std::string>  ("outPath","/uscms/home/eberry/CMSSW_2_0_8/src/Analyzers/HcalDigiAnalyzer");
  outSuffix_        = iConfig.getUntrackedParameter<std::string>  ("outSuffix","");
  subdetName_       = iConfig.getUntrackedParameter<std::string>  ("subdetName","NO");

  //-----------------------------------------------------
  // Determine the subdetector from the .cfg file entry
  //-----------------------------------------------------
  
  if (subdetName_ == "HB") { subdet_ = 1; }
  if (subdetName_ == "HE") { subdet_ = 2; }
  if (subdetName_ == "HO") { subdet_ = 3; }
  if (subdetName_ == "HF") { subdet_ = 4; }
  if (subdetName_ == "NO") { 
    printf("Please pick a subdetector\n");
     exit(0);
  }

  //-----------------------------------------------------
  // Set up pulse correction
  //-----------------------------------------------------

  firstRHSample_          = 4;      // Configured for HBHE for now, may have user input later 
  rhSamplesToAdd_         = 4;      // Configured for HBHE for now, may have user input later
  phaseNS_                = 13.0;   // Configured for HBHE for now, may have user input later
  
  pulseCorrectionPtr = std::auto_ptr<HcalPulseContainmentCorrection>(new HcalPulseContainmentCorrection(rhSamplesToAdd_,phaseNS_,MaximumFractionalError_));
  
  //-----------------------------------------------------
  // Determine the path for the output
  //-----------------------------------------------------
  
  rootFile_ = outPath_ + "HcalFinalOutput_" + subdetName_+ outSuffix_ + ".root";
  
  if (verbose_) cout << "my root file is: " << rootFile_ << endl;
  
  //-----------------------------------------------------
  // Initialize the root tree you're going to fill
  //-----------------------------------------------------

  m_fillDigi.init(rootFile_, &m_digiTree);
  

}

template < typename T >
HcalDigiAnalyzer<T>::~HcalDigiAnalyzer(){}

template < typename T >
void HcalDigiAnalyzer<T>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //-----------------------------------------------------
  // Declare namespaces
  //-----------------------------------------------------

  using namespace edm; 
  using namespace std; 
  
  //-----------------------------------------------------
  // Initialize the Tree values
  //-----------------------------------------------------
  
  m_digiTree.init();

  //-----------------------------------------------------
  // Fill the run and event number information
  //-----------------------------------------------------

  int run   = iEvent.id().run();
  int event = iEvent.id().event();  

  m_digiTree.run   = run;
  m_digiTree.event = event;

  if (scanForSpikes_ && subdet_ == 1)
    if (verbose_) printf("Run = %d, Event = %d\n",run,event);

  //-----------------------------------------------------
  // Get setup information
  //-----------------------------------------------------
  
  ESHandle<HcalTPGCoder> inputCoder;
  ESHandle<CaloGeometry> geometry ;
  ESHandle<HcalDbService> conditions;
  
  iSetup.get<HcalTPGRecord>().get(inputCoder);
  iSetup.get<HcalDbRecord>().get(conditions);
  iSetup.get<IdealGeometryRecord>().get(geometry);

  const HcalPulseContainmentCorrection* pulseCorrection = pulseCorrectionPtr.get();
  const HcalQIEShape* shape = conditions->getHcalShape();
  CaloSamples tool;

  //-----------------------------------------------------
  // Declare important values
  //-----------------------------------------------------
  
  char subdetName[5];

  if (subdet_ == 1) sprintf(subdetName,"HB");
  if (subdet_ == 2) sprintf(subdetName,"HE");
  if (subdet_ == 3) sprintf(subdetName,"HO");
  if (subdet_ == 4) sprintf(subdetName,"HF");

  HcalTrigTowerGeometry theTrigTowerGeometry;

  vector<HcalTrigTowerDetId> trigTowerDetIds;
  HcalTrigTowerDetId trigTowerDetId;

  vector<int> detIdMap;

  int ndigis= 0;
  int count = 0;
  int ntptsample;
  int ntrigprim;
  int id;

  float fC_amp, GeV_amp, temp_amp;
  float correction;

  int ieta, iphi, depth;

  int nTrigTowerDetIds;
  int trigPrimDigiSize;
  int compressedEt;

  bool foundSpike;
  bool digiExists;				
  bool foundDigiBit;
  
  typename T::const_iterator ihcal;

  //-----------------------------------------------------
  // Get handles of digis and trig prims
  //-----------------------------------------------------

  Handle < T > hcalDigiCol;  
  bool hcalDigiTagExists = iEvent.getByLabel(hcalDigiTag_,hcalDigiCol);
  if (!hcalDigiTagExists){
    LogWarning("HcalDigiAnalyzer") << "Could not extract HCAL digis!";
    return;
  }
  
  Handle<HcalTrigPrimDigiCollection> HCALTrigPrimDigis;
  bool hcalTrigPrimDigiTagExists = iEvent.getByLabel(hcalTrigPrimTag_,HCALTrigPrimDigis);
  if (!hcalTrigPrimDigiTagExists){
    LogWarning("HcalDigiAnalyzer") << "Could not extract trigger primitives!";
    return;
  }

  //-----------------------------------------------------
  // Loop over the digi collections
  //-----------------------------------------------------

  for (ihcal  = hcalDigiCol->begin();
       ihcal != hcalDigiCol->end();
       ihcal ++){

    //-----------------------------------------------------
    // Count the number of digis in this collection
    //-----------------------------------------------------

    count++;

    //-----------------------------------------------------
    // Get the detector ID.  We'll need this to get cell
    // location information.
    //-----------------------------------------------------
    
    HcalDetId cell(ihcal->id()); 

    //-----------------------------------------------------
    // Make sure this is the right subdetector
    //-----------------------------------------------------

    if(cell.subdet() == subdet_) {

      ndigis++;
            
      //-----------------------------------------------------
      // Get cell location information (do not sort)
      //-----------------------------------------------------
      
      id = (int)(*ihcal).id();

      ieta  = cell.ieta();
      iphi  = cell.iphi();
      depth = cell.depth();

      m_digiTree.h_depth[count-1] = cell.depth();
      m_digiTree.h_ieta [count-1] = cell.ieta();
      m_digiTree.h_iphi [count-1] = cell.iphi();
      m_digiTree.h_id   [count-1] = (int) id;

      //-----------------------------------------------------
      // Get digi size information (number of [pre]samples)
      //-----------------------------------------------------
      
      m_digiTree.h_psam [ieta+41][iphi][depth] = (*ihcal).presamples();
      m_digiTree.h_size [ieta+41][iphi][depth] = (*ihcal).size();
      
      //-----------------------------------------------------
      // Now comes the coding...
      //-----------------------------------------------------
      
      const HcalCalibrations& calibrations = conditions->getHcalCalibrations(cell);
      const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
      HcalCoderDb coder (*channelCoder, *shape); 
      coder.adc2fC(*ihcal,tool);
      
      //-----------------------------------------------------
      // Loop over the time samples from this hit
      //-----------------------------------------------------
      
      foundSpike = false;

      fC_amp = 0.0;
      GeV_amp = 0.0;

      for( int ii=0; ii<tool.size(); ii++ ) { 

	int capid    = (*ihcal)[ii].capid();
	float ped    = calibrations.pedestal(capid);	
	float gain   = calibrations.rawgain(capid);  
	float rcgain = calibrations.respcorrgain(capid);
		
	//-----------------------------------------------------
	// Get rec hit amplitude in GeV 
	// (should agree with rec hit energy)
	//-----------------------------------------------------

	if ( ii >= firstRHSample_ &&                 
	     ii <  firstRHSample_ + rhSamplesToAdd_ ){  
	  
	  // Subtract fC pedestal from digi
	  temp_amp = tool[ii] - ped;
	  
	  // Store the fC amplitude
	  fC_amp  += temp_amp;

	  // Convert the fC amplitude to GeV
	  temp_amp *= rcgain;
	  
	  // Store the GeV amplitude
	  GeV_amp += temp_amp;
	  
	}

	m_digiTree.h_adc   [ieta+41][iphi][depth][ii] = (*ihcal)[ii].adc(); 
	m_digiTree.h_fC    [ieta+41][iphi][depth][ii] = tool[ii];	  
	m_digiTree.h_ped   [ieta+41][iphi][depth][ii] = ped;
	m_digiTree.h_gain  [ieta+41][iphi][depth][ii] = gain;
	m_digiTree.h_rcgain[ieta+41][iphi][depth][ii] = rcgain;
	m_digiTree.h_capid [ieta+41][iphi][depth][ii] = capid;
	m_digiTree.h_fiber [ieta+41][iphi][depth][ii] = (*ihcal)[ii].fiber();
	m_digiTree.h_fchan [ieta+41][iphi][depth][ii] = (*ihcal)[ii].fiberChan();	

	//-----------------------------------------------------
	// Check for spikes and print out warnings
	// (must set verbose to 'true' in .cfg file
	//-----------------------------------------------------
	
	if (scanForSpikes_){	
	  
	  if ( (tool[ii] > 25.0 && subdet_ == 4) || 
	       (tool[ii] > 10.0 && subdet_ == 1) ){
	    
	    foundSpike = true;
	    
	    if (verbose_){
	      
	      printf("%s SPIKE?: %f\n",subdetName,tool[ii]);
	      printf("  run = %d, event = %d\n",run,event);
	      printf("  ieta = %d, iphi = %d, depth = %d, ts = %d\n",ieta,iphi,depth,ii);
	      
	      if (tool[ii] > 100.0){		
		
		printf("%s SPIKE!!!: %f\n",subdetName,tool[ii]);
		printf("  run = %d, event = %d\n",run,event);
		printf("  ieta = %d, iphi = %d, depth = %d\n",ieta,iphi,depth);
	      }	
	    }
	  }	    
	}
      }
      
      //-----------------------------------------------------
      // Correct for RH amplitude for phase
      // Store both fC and GeV amplitudes
      //-----------------------------------------------------

      if (pulseCorrection != 0){
	correction = pulseCorrection -> getCorrection(fC_amp);
	GeV_amp *= correction;
      }
      else correction = 0.0;
      
      m_digiTree.h_correction [ieta+41][iphi][depth] = correction;
      m_digiTree.h_rh_GeV_amp [ieta+41][iphi][depth] = GeV_amp;
      m_digiTree.h_rh_fC_amp  [ieta+41][iphi][depth] = fC_amp;
      
      //-----------------------------------------------------
      // Done looping over the time samples...
      //
      // If you DID find a spike, get the trigger primitives
      //-----------------------------------------------------
      
      if (foundSpike){

	//-----------------------------------------------------
	// According to Jeremy Mans, you have to use the
	// trigger tower geometry to get trig tower det ids.
	//
	// There should be no more than 2 trig tower det ids
	// per HcalDetId (usually only 1, never less than 1)
	//-----------------------------------------------------

	trigTowerDetIds  = theTrigTowerGeometry.towerIds(cell);
	nTrigTowerDetIds = trigTowerDetIds.size();

	assert(nTrigTowerDetIds == 1 || nTrigTowerDetIds == 2 );

	//-----------------------------------------------------
	// Loop over the trigger primitive det id(s) we found
	//-----------------------------------------------------

	ntrigprim = 0;

	for (int iTrigTowerDetId = 0; iTrigTowerDetId < nTrigTowerDetIds; iTrigTowerDetId++){

	  //-----------------------------------------------------
	  // Get the trigger primitive id
	  //-----------------------------------------------------

	  trigTowerDetId = trigTowerDetIds[iTrigTowerDetId];

	  HcalTrigPrimDigiCollection::const_iterator iTrigPrimDigi = (*HCALTrigPrimDigis).find(trigTowerDetId);

	  //-----------------------------------------------------
	  // Did we find the trig prim digi?
	  // SortedCollections return end() if they can't find
	  // the item you're looking for...
	  //-----------------------------------------------------

	  digiExists = true;
	  if (iTrigPrimDigi == (*HCALTrigPrimDigis).end()) digiExists = false;

	  //-----------------------------------------------------
	  // If the trig prim digi exists (and it should), 
	  // loop over the time samples
	  //-----------------------------------------------------

	  ntptsample = 0;
	  
	  if (digiExists){
	    
	    foundDigiBit = 1;
	    
	    trigPrimDigiSize = (int) (*iTrigPrimDigi).size();
	    
	    m_digiTree.t_size  [ieta+41][iphi][depth][ntrigprim] = trigPrimDigiSize;  
	    m_digiTree.t_psam  [ieta+41][iphi][depth][ntrigprim] = (*iTrigPrimDigi).presamples();
	    m_digiTree.t_ieta  [ieta+41][iphi][depth][ntrigprim] = (int) trigTowerDetId.ieta();
	    m_digiTree.t_iphi  [ieta+41][iphi][depth][ntrigprim] = (int) trigTowerDetId.iphi();
	    m_digiTree.t_subdet[ieta+41][iphi][depth][ntrigprim] = (int) trigTowerDetId.subdet();	    
	    
	    for (int iTrigPrimSample = 0; iTrigPrimSample < trigPrimDigiSize; iTrigPrimSample++){
	      
	      HcalTriggerPrimitiveSample trigPrimSample = (*iTrigPrimDigi)[iTrigPrimSample];
	      
	      compressedEt = (int) trigPrimSample.compressedEt();
	      
	      m_digiTree.t_cET[ieta+41][iphi][depth][ntrigprim][ntptsample] = (int) compressedEt;
	      
	      ntptsample++;
	      
	    }
	  }

	  if (!digiExists){
	    foundDigiBit = 0;
	    printTrigPrimError(trigTowerDetId,cell,nTrigTowerDetIds,ntrigprim);
	  }
	  
	  m_digiTree.t_found [ieta+41][iphi][depth][ntrigprim] = foundDigiBit;		
	  m_digiTree.t_ntpts [ieta+41][iphi][depth][ntrigprim] = ntptsample;
	  ntrigprim++;

	}
	
	m_digiTree.t_spike[ieta+41][iphi][depth] = 1;
	m_digiTree.t_ntp  [ieta+41][iphi][depth] = ntrigprim;
	
      }
      
      else {
	m_digiTree.t_spike[ieta+41][iphi][depth] = 0;
	m_digiTree.t_ntp  [ieta+41][iphi][depth] = 0;
      }
    }
  }
    
  m_digiTree.nchn = count;
  m_digiTree.nhit = ndigis;
  
  m_fillDigi.fill();

}

template < typename T >
void HcalDigiAnalyzer<T>::printTrigPrimError
(HcalTrigTowerDetId trigTowerDetId,HcalDetId cell, int nTrigTowerDetIds, int ntrigprim){

  printf("***NO TRIG PRIM DIGI FOUND!!!***\n");
  printf("  We expected %d digis.\n",nTrigTowerDetIds);
  printf("  We did not find digi #%d\n",ntrigprim+1);
  printf("-----------------------------------------\n");
  printf("  Trig prim id had:\n");
  printf("    subdet = %d\n",trigTowerDetId.subdet());
  printf("    zside  = %d\n",trigTowerDetId.zside());
  printf("    ieta   = %d\n",trigTowerDetId.ieta());
  printf("    |ieta| = %d\n",trigTowerDetId.ietaAbs());
  printf("    iphi   = %d\n",trigTowerDetId.iphi());
  printf("-----------------------------------------\n");
  printf("  Hcal det id had:\n");
  printf("    subdet = %d\n",cell.subdet());
  printf("    zside  = %d\n",cell.zside());
  printf("    ieta   = %d\n",cell.ieta());
  printf("    |ieta| = %d\n",cell.ietaAbs());
  printf("    iphi   = %d\n",cell.iphi());
  printf("-----------------------------------------\n");
  

}

template < typename T >
void HcalDigiAnalyzer<T>::beginJob(const edm::EventSetup&){}

template < typename T >
void HcalDigiAnalyzer<T>::endJob() 
{
   using namespace edm;
   using namespace std;

   m_fillDigi.finalize();
}

#endif
