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
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

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
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/HcalTowerAlgo/src/HcalHardcodeGeometryData.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
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
#include <vector>

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
  
  template <class Digi>
  void analyzeThisDigi(Digi& digi, 
		       const HcalCoder& coder, const HcalCalibrations& calibrations,
		       const HcalPulseContainmentCorrection *corr, int count);
  
  float getThreshold(HcalDetId &id);
  void  readEventList();

  //-----------------------------------------------------
  // Root tree objects
  //-----------------------------------------------------

  DigiTree     m_digiTree;
  FillDigiTree m_fillDigi;

  //-----------------------------------------------------
  // Event list
  //-----------------------------------------------------
    
  std::string      eventListFile_;
  std::vector<int> eventList_;

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
  edm::InputTag caloTowerTag_;

  //-----------------------------------------------------
  // edm::ESHandles
  //-----------------------------------------------------

  edm::ESHandle<HcalTopology> hcalTopology_;
  edm::ESHandle<HcalDbService> conditions_;

  //-----------------------------------------------------
  // edm::Handles
  //-----------------------------------------------------

  edm::Handle <T> hcalDigiCol_;  
  edm::Handle <CaloTowerCollection> caloTowerCol_;

  //-----------------------------------------------------
  // Output flags
  //-----------------------------------------------------

  bool doRHPulseCorrect_;
  bool scanForSpikes_;
  bool verbose_;  
  bool useEventList_;
  int subdet_;

  //-----------------------------------------------------
  // Scheme B CaloTower creation thresholds
  //-----------------------------------------------------

  float theHOthreshold0     ;
  float theHOthresholdPlus1 ;
  float theHOthresholdMinus1;
  float theHOthresholdPlus2 ;
  float theHOthresholdMinus2;
  float theHBthreshold      ;
  float theHF1threshold     ;
  float theHF2threshold     ;
  float theHESthreshold     ;
  float theHEDthreshold     ;

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
  useEventList_     = iConfig.getUntrackedParameter<bool>         ("useEventList",false);
  hcalTrigPrimTag_  = iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigPrimTag");
  hcalDigiTag_      = iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag"); 
  caloTowerTag_     = iConfig.getUntrackedParameter<edm::InputTag>("caloTowerTag");
  outPath_          = iConfig.getUntrackedParameter<std::string>  ("outPath","/uscms/home/eberry/3DayLifetime/");
  outSuffix_        = iConfig.getUntrackedParameter<std::string>  ("outSuffix","");
  subdetName_       = iConfig.getUntrackedParameter<std::string>  ("subdetName","NO");
  eventListFile_    = iConfig.getUntrackedParameter<std::string>  ("eventList","/uscms/home/eberry/HcalNoiseData/CRAFT_on_QCD_NewHcal_v4a_evlist.txt");

  //-----------------------------------------------------
  // Determine the subdetector from the .cfg file entry
  //-----------------------------------------------------
  
  if (subdetName_ == "HB") { 
    subdet_         = 1; 
    phaseNS_        = 13.0;
    firstRHSample_  = 4;
    rhSamplesToAdd_ = 4;
  }
  
  if (subdetName_ == "HE") { 
    subdet_         = 2; 
    phaseNS_        = 13.0;
    firstRHSample_  = 4;
    rhSamplesToAdd_ = 4;
  }
  
  if (subdetName_ == "HO"){
    subdet_         = 3;
    phaseNS_        = 13.0;
    firstRHSample_  = 4;
    rhSamplesToAdd_ = 4;
  }

  if (subdetName_ == "HF") { 
    subdet_         = 4; 
    phaseNS_        = 13.0;
    firstRHSample_  = 3;
    rhSamplesToAdd_ = 1;
  }
  
  if (subdetName_ == "NO") { 
    printf("Please pick a subdetector\n");
     exit(0);
  }

  //-----------------------------------------------------
  // CaloTower Scheme B creation thresholds
  //-----------------------------------------------------

  theHOthreshold0      = 1.1;
  theHOthresholdPlus1  = 1.1;
  theHOthresholdMinus1 = 1.1;
  theHOthresholdPlus2  = 1.1;
  theHOthresholdMinus2 = 1.1;
  theHBthreshold       = 0.9;
  theHF1threshold      = 1.2;
  theHF2threshold      = 1.8;
  theHESthreshold      = 1.4;
  theHEDthreshold      = 1.4;
  
  //-----------------------------------------------------
  // Set up pulse correction
  //-----------------------------------------------------
  
  if (doRHPulseCorrect_)
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
  // Fill the run and event number information
  //-----------------------------------------------------

  int run   = iEvent.id().run();
  int event = iEvent.id().event();  

  //-----------------------------------------------------
  // If using an event list, only proceed if this event
  // is on the list
  //-----------------------------------------------------
  
  if (useEventList_){
    vector<int>::iterator eventItr;
    eventItr = find(eventList_.begin(),eventList_.end(),event);
    
    cout << "Run = " << run << endl;

    if (eventItr == eventList_.end()) {
      cout << "We DID NOT find " << event << " in our event list" << endl;
      return;
    }

    else {
      cout << "We found " << event << " in our event list!" << endl;
    }
  }

  //-----------------------------------------------------
  // Initialize the Tree values
  //-----------------------------------------------------
  
  m_digiTree.init();
  m_digiTree.run   = run;
  m_digiTree.event = event;

  if (scanForSpikes_ && subdet_ == 1)
    if (verbose_) printf("Run = %d, Event = %d\n",run,event);

  //-----------------------------------------------------
  // Get setup information
  //-----------------------------------------------------
  
  iSetup.get<IdealGeometryRecord>().get(hcalTopology_);
  iSetup.get<HcalDbRecord>().get(conditions_);

  const HcalPulseContainmentCorrection* pulseCorrection = pulseCorrectionPtr.get();
  const HcalQIEShape* shape = conditions_->getHcalShape();

  //-----------------------------------------------------
  // Loop over calo tower collection
  //-----------------------------------------------------

  int ieta, iphi;
  int nct = 0;

  bool caloTowerTagExists = iEvent.getByLabel(caloTowerTag_,caloTowerCol_);
  if (!caloTowerTagExists){
    LogWarning("HcalDigiAnalyzer") << "Could not extract CaloTowers!";
  } else {
    
    for (CaloTowerCollection::const_iterator ictow  = caloTowerCol_ -> begin();
	 ictow != caloTowerCol_ -> end();
	 ictow++){
      
      ieta = (int) (*ictow).id().ieta();
      iphi = (int) (*ictow).id().iphi();
      
      m_digiTree.ct_ieta[nct] = ieta;
      m_digiTree.ct_iphi[nct] = iphi;

      nct++;
      
    }
    
  }

  m_digiTree.nct = nct;
  
  //-----------------------------------------------------
  // Loop over the digi collections
  //-----------------------------------------------------

  int ndigis= 0;
  int count = 0;

  CaloSamples tool;  
  typename T::const_iterator ihcal;
 
  bool hcalDigiTagExists = iEvent.getByLabel(hcalDigiTag_,hcalDigiCol_);
  if (!hcalDigiTagExists){
    LogWarning("HcalDigiAnalyzer") << "Could not extract HCAL digis!";
  } else {

    for (ihcal  = hcalDigiCol_->begin();
	 ihcal != hcalDigiCol_->end();
	 ihcal ++){
            
      HcalDetId cell(ihcal->id()); 
            
      if(cell.subdet() == subdet_) {
	
	const HcalCalibrations& calibrations = conditions_->getHcalCalibrations(cell);
	const HcalQIECoder* channelCoder = conditions_->getHcalCoder(cell);
	HcalCoderDb coder (*channelCoder, *shape); 
		
	analyzeThisDigi(*ihcal, coder, calibrations, pulseCorrection,ndigis);
	
	ndigis++;
	
      }
      
      count++;
      
    }
  }
  
  m_digiTree.nchn = count;
  m_digiTree.nhit = ndigis;
  
  m_fillDigi.fill();

}

template <typename T>
template <class Digi>
void HcalDigiAnalyzer<T>::analyzeThisDigi(Digi& digi, 
						 const HcalCoder& coder, const HcalCalibrations& calibrations,
						 const HcalPulseContainmentCorrection *corr, int count
						 ){  

  //-----------------------------------------------------
  // Declare useful values
  //-----------------------------------------------------
  
  float signal;
  float amplitudeInFC  = 0.0;
  float amplitudeInGeV = 0.0;
  float correction;
  
  //-----------------------------------------------------
  // Get cell location information (do not sort)
  //-----------------------------------------------------
  
  HcalDetId cell(digi.id());
  
  int id    = (int) digi.id();    
  int ieta  = (int) cell.ieta();
  int iphi  = (int) cell.iphi();
  int depth = (int) cell.depth();

  m_digiTree.h_depth[count] = depth;
  m_digiTree.h_ieta [count] = ieta;
  m_digiTree.h_iphi [count] = iphi;
  m_digiTree.h_id   [count] = id;
  
  bool validID = HcalDetId::validDetId((HcalSubdetector)(subdet_),ieta,iphi,depth);
  if (!validID) {
    std::cout << "Invalid HCAL detector ID: subdetector = " << subdet_ << ", ieta = " << ieta << ", iphi = " << iphi << ", depth = " << depth << std::endl;
  }
  
  //-----------------------------------------------------
  // Declare a CaloSamples object (fC time samples),
  // and fill it
  //-----------------------------------------------------
  
  CaloSamples tool;
  coder.adc2fC(digi,tool);

  //-----------------------------------------------------
  // Get digi size information (number of [pre]samples)
  //-----------------------------------------------------
      
  m_digiTree.h_psam [ieta+41][iphi][depth] = digi.presamples();
  m_digiTree.h_size [ieta+41][iphi][depth] = digi.size();
  
  //-----------------------------------------------------
  // Loop over the time samples from this CaloSample
  //-----------------------------------------------------
  
  for( int ii=0; ii<tool.size(); ii++ ) { 
    
    int   capid  = (int)   digi[ii].capid();
    float ped    = (float) calibrations.pedestal(capid);	
    float gain   = (float) calibrations.rawgain(capid);  
    float rcgain = (float) calibrations.respcorrgain(capid);
    float fC     = (float) tool[ii];
    int   adc    = (int)   digi[ii].adc();

    m_digiTree.h_adc   [ieta+41][iphi][depth][ii] = adc;
    m_digiTree.h_fC    [ieta+41][iphi][depth][ii] = fC;
    m_digiTree.h_ped   [ieta+41][iphi][depth][ii] = ped;
    m_digiTree.h_gain  [ieta+41][iphi][depth][ii] = gain;
    m_digiTree.h_rcgain[ieta+41][iphi][depth][ii] = rcgain;
    m_digiTree.h_capid [ieta+41][iphi][depth][ii] = capid;
    m_digiTree.h_fiber [ieta+41][iphi][depth][ii] = digi[ii].fiber();
    m_digiTree.h_fchan [ieta+41][iphi][depth][ii] = digi[ii].fiberChan();	
    
    //-----------------------------------------------------
    // Calculate the amplitude for this digi
    //-----------------------------------------------------

    if ( ii >= firstRHSample_ && ii < firstRHSample_ + rhSamplesToAdd_){
      signal          = fC - ped; // Pedestal subtraction
      amplitudeInFC  += signal;   // Calculate fC amplitude
      signal         *= rcgain;   // Convert from fC to GeV
      amplitudeInGeV += signal;   // Calculate GeV amplitude
    }            
  } // end loop over time samples 
  
  if (corr!=0 && doRHPulseCorrect_) {
    correction = corr->getCorrection(amplitudeInFC);
    amplitudeInGeV *= correction;
  }
  else correction = 0.0;

  float threshold = getThreshold(cell); 
  if (threshold == -999){
    std::cout << "WARNING: Could not extract a threshold for this tower" << std::endl;
  }

  m_digiTree.h_correction [ieta+41][iphi][depth] = correction;
  m_digiTree.h_rh_GeV_amp [ieta+41][iphi][depth] = amplitudeInGeV;
  m_digiTree.h_rh_fC_amp  [ieta+41][iphi][depth] = amplitudeInFC;
  m_digiTree.h_threshold  [ieta+41][iphi][depth] = threshold;
  
}

template < typename T >
float HcalDigiAnalyzer<T>::getThreshold(HcalDetId & hcalDetId){
    
  float threshold = -999;

  HcalSubdetector subdet = hcalDetId.subdet();
  
  if(subdet == HcalBarrel) {
    threshold = theHBthreshold;
  }
  
  else if(subdet == HcalEndcap) {
    // check if it's single or double tower
    if(hcalDetId.ietaAbs() < hcalTopology_->firstHEDoublePhiRing()) {
      threshold = theHESthreshold;
    }
    else {
      threshold = theHEDthreshold;
    }
  }
  
  else if(subdet == HcalOuter) {
    //check if it's ring 0 or +1 or +2 or -1 or -2
    if(hcalDetId.ietaAbs() <= 4) threshold = theHOthreshold0;
    else if(hcalDetId.ieta() < 0) {
      // set threshold for ring -1 or -2
      threshold = (hcalDetId.ietaAbs() <= 10) ?  theHOthresholdMinus1 : theHOthresholdMinus2;
    } else {
      // set threshold for ring +1 or +2
      threshold = (hcalDetId.ietaAbs() <= 10) ?  theHOthresholdPlus1 : theHOthresholdPlus2;
    }
  } 
  
  else if(subdet == HcalForward) {
    if(hcalDetId.depth() == 1) {
      threshold = theHF1threshold;
    } else {
      threshold = theHF2threshold;
    }
  }
  
  return threshold;
  
}

template < typename T >
void HcalDigiAnalyzer<T>::readEventList(){

  FILE *file;

  const char *eventListFile = eventListFile_.c_str();

  file = fopen(eventListFile,"r");
  
  if (file == NULL){
    std::cout << "There is not file with the name: " << std::endl;
    std::cout << eventListFile_ << std::endl;
    exit(0);
  }
  
  int  rowNumber, run, event;
  char runString  [100];
  char eventString[100];
  
  while (true){

    int output = fscanf(file,"%d %s %d %s %d",&rowNumber,runString,&run, eventString, &event);    
    if (output == EOF) break;   
    
    eventList_.push_back(event);

  }
  
  fclose (file);

}

template < typename T >
void HcalDigiAnalyzer<T>::beginJob(const edm::EventSetup&){
  
  if (useEventList_) readEventList();

}

template < typename T >
void HcalDigiAnalyzer<T>::endJob() 
{
   using namespace edm;
   using namespace std;

   m_fillDigi.finalize();
}

#endif
