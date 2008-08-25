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

#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/HcalTowerAlgo/src/HcalHardcodeGeometryData.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

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
  
  // ----------member data ---------------------------
  DigiTree     m_digiTree;
  FillDigiTree m_fillDigi;
  
  bool adc2fcDone_;
  //std::vector<double> perpIeta_;
  
  std::string subdetName_;
  std::string outPath_;
  std::string outSuffix_;
  std::string rootFile_;

  edm::InputTag hcalTrigPrimTag_;
  edm::InputTag hcalUnsuppressedDigiTag_;
  edm::InputTag hcalDigiTag_;

  bool scanForSpikes_;

  double noiseArray_ [4][3000][10];
  
  int subdet_;
  
};

template < typename T >
HcalDigiAnalyzer<T>::HcalDigiAnalyzer(const edm::ParameterSet& iConfig):
  adc2fcDone_(false),
  outPath_(iConfig.getUntrackedParameter<std::string>("outPath","/uscms/home/eberry/CMSSW_2_0_8/src/Analyzers/HcalDigiAnalyzer")),
  outSuffix_(iConfig.getUntrackedParameter<std::string>("outSuffix",""))    
{

  using namespace std;
  
  std::vector<std::string> defaultInFiles;

  scanForSpikes_ = iConfig.getUntrackedParameter<bool>("scanForSpikes",false);
  
  subdetName_ = iConfig.getUntrackedParameter<std::string>("subdetName","NO");
  
  if (subdetName_ == "HB") { subdet_ = 1; }
  if (subdetName_ == "HE") { subdet_ = 2; }
  if (subdetName_ == "HO") { subdet_ = 3; }
  if (subdetName_ == "HF") { subdet_ = 4; }
  if (subdetName_ == "NO") { 
    printf("Please pick a subdetector\n");
     exit(0);
  }
  
  rootFile_ = outPath_ + "HcalFinalOutput_" + subdetName_+ outSuffix_ + ".root";
  
  cout << "my root file is: " << rootFile_ << endl;
  
  m_fillDigi.init(rootFile_, &m_digiTree);

  const edm::InputTag dHcalTrigPrimTag("hcalTriggerPrimitiveDigis");
  hcalTrigPrimTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalTrigPrimTag",dHcalTrigPrimTag);

  //-----------------------------------------------------
  // IMPORTANT:
  // 
  // In SIMULATION (as of CMSSW_2_0_8):
  // 
  // -- ZERO-SUPPRESSED digis are stored with the tag 
  //    'hcalDigis'
  // -- UNSUPPRESSED digis are stored with the tag
  //    'hcalUnsuppressedDigis'
  //
  // In DATA (as of CRUZET3):
  // 
  // -- Only UNSUPPRESSED digis are stored
  // -- They must be unpacked from the FED
  // -- When they are unpacked, they are stored with
  //    the tag 'hcalDigis', the SAME tag as the 
  //    SUPPRESSED digis in simulation.
  //
  // The trick, then, is to see whether we have a digi
  // object with the tag 'hcalUnsuppressedDigis'
  //
  // - If the tag EXISTS, this is SIMULATION, and we
  //   have the digis we want. Done.
  // - If the tag DOES NOT EXIST, this is DATA, and
  //   we need to look for the tag 'hcalDigis'
  //-----------------------------------------------------

  const edm::InputTag dHcalDigiTag("hcalDigis");
  hcalDigiTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalDigiTag",dHcalDigiTag);

  const edm::InputTag dHcalUnsuppressedDigiTag("hcalUnsuppressedDigis");
  hcalUnsuppressedDigiTag_ = iConfig.getUntrackedParameter<edm::InputTag>("hcalUnsuppressedDigiTag",dHcalUnsuppressedDigiTag);

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
  m_digiTree.nadc  = 0;

  if (scanForSpikes_ && subdet_ == 1)
    printf("Run = %d, Event = %d\n",run,event);

  //-----------------------------------------------------
  // Get setup information
  //-----------------------------------------------------
  
  ESHandle<HcalTPGCoder> inputCoder;
  ESHandle<CaloGeometry> geometry ;
  ESHandle<HcalDbService> conditions;
  
  iSetup.get<HcalTPGRecord>().get(inputCoder);
  iSetup.get<HcalDbRecord>().get(conditions);
  iSetup.get<IdealGeometryRecord>().get(geometry);

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

  int nTrigTowerDetIds;
  int trigPrimDigiSize;
  int compressedEt;

  bool foundSpike;
  bool digiExists;				
  bool foundDigiBit;
  
  //-----------------------------------------------------
  // Get Handle of digis
  //-----------------------------------------------------

  typename T::const_iterator ihcal;

  bool hcalUnsuppressedDigiTagExists = false;
  bool hcalDigiTagExists = false;

  Handle < T > hcalDigiCol;
  hcalUnsuppressedDigiTagExists = iEvent.getByLabel(hcalUnsuppressedDigiTag_,hcalDigiCol);
  
  if (!hcalUnsuppressedDigiTagExists)
    hcalDigiTagExists = iEvent.getByLabel(hcalDigiTag_,hcalDigiCol);

  if (!hcalUnsuppressedDigiTagExists &&
      !hcalDigiTagExists){
    LogWarning("HcalDigiAnalyzer") << "Could not extract HCAL digis with EITHER tag!";
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
      // Get the Cell Geometry and make sure it's sane
      //-----------------------------------------------------
      
      const CaloCellGeometry* cellGeometry =
	geometry->getSubdetectorGeometry (cell)->getGeometry (cell);
      if( cellGeometry == 0 ) LogInfo("") << "No cell geometry " << cell.rawId();
      
      //-----------------------------------------------------
      // Get cell location information
      //-----------------------------------------------------

      double fEta   = cellGeometry->getPosition().eta();
      double fPhi   = cellGeometry->getPosition().phi();
      
      id = (int)(*ihcal).id();

      m_digiTree.h_psam [count-1] = (*ihcal).presamples();
      m_digiTree.h_size [count-1] = (*ihcal).size();
      
      m_digiTree.h_depth[count-1] = cell.depth();
      m_digiTree.h_eta  [count-1] = fEta;
      m_digiTree.h_phi  [count-1] = fPhi; 
      m_digiTree.h_ieta [count-1] = cell.ieta();
      m_digiTree.h_iphi [count-1] = cell.iphi();
      m_digiTree.h_id   [count-1] = (int) id;
      
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

      for( int ii=0; ii<tool.size(); ii++ ) { 

	int capid  = (*ihcal)[ii].capid();
	float ped  = calibrations.pedestal(capid);
	float gain = calibrations.rawgain(capid);  
	
	m_digiTree.h_adc     [count-1][ii] = (*ihcal)[ii].adc();
	m_digiTree.h_fC      [count-1][ii] = tool[ii];	  
	m_digiTree.h_ped     [count-1][ii] = ped;
	m_digiTree.h_gain    [count-1][ii] = gain;
	m_digiTree.h_capid   [count-1][ii] = capid;
	m_digiTree.h_fiber   [count-1][ii] = (*ihcal)[ii].fiber();
	m_digiTree.h_fchan   [count-1][ii] = (*ihcal)[ii].fiberChan();
	
	//-----------------------------------------------------
	// Check for spikes and print out warnings
	//-----------------------------------------------------
	
	if (scanForSpikes_){	
	  
	  if ( (tool[ii] > 25.0 && subdet_ == 4) || 
	       (tool[ii] > 10.0 && subdet_ == 1) ){
	    
	    foundSpike = true;
	    
	    printf("%s SPIKE?: %f\n",subdetName,tool[ii]);
	    printf("  run = %d, event = %d\n",run,event);
	    printf("  ieta = %d, iphi = %d, depth = %d, ts = %d\n",cell.ieta(),cell.iphi(),cell.depth(),ii);
	    
	    if (tool[ii] > 100.0){		
	      printf("%s SPIKE!!!: %f\n",subdetName,tool[ii]);
	      printf("  run = %d, event = %d\n",run,event);
	      printf("  ieta = %d, iphi = %d, depth = %d\n",cell.ieta(),cell.iphi(),cell.depth());
	    }
	  }	    
	}
      }
      
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
	    
	    m_digiTree.t_size [count-1][ntrigprim] = trigPrimDigiSize;  
	    m_digiTree.t_psam [count-1][ntrigprim] = (*iTrigPrimDigi).presamples();
	    
	    for (int iTrigPrimSample = 0; iTrigPrimSample < trigPrimDigiSize; iTrigPrimSample++){
	      
	      HcalTriggerPrimitiveSample trigPrimSample = (*iTrigPrimDigi)[iTrigPrimSample];
	      
	      compressedEt = (int) trigPrimSample.compressedEt();
	      
	      m_digiTree.t_cET[count-1][ntrigprim][ntptsample] = (int) compressedEt;
	      
	      ntptsample++;
	      
	    }
	  }

	  if (!digiExists){
	    foundDigiBit = 0;
	    printTrigPrimError(trigTowerDetId,cell,nTrigTowerDetIds,ntrigprim);
	  }
	  
	  m_digiTree.t_found [count-1][ntrigprim] = foundDigiBit;		
	  m_digiTree.t_ntpts [count-1][ntrigprim] = ntptsample;
	  ntrigprim++;

	}
	
	m_digiTree.t_spike[count-1] = 1;
	m_digiTree.t_ntp  [count-1] = ntrigprim;
	
      }
      
      else {
	m_digiTree.t_spike[count-1] = 0;
	m_digiTree.t_ntp  [count-1] = 0;
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
