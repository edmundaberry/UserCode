#include "SimCalorimetry/HcalSimAlgos/interface/HcalAmplifier.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameters.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CondFormats/HcalObjects/interface/HcalPedestal.h"
#include "CondFormats/HcalObjects/interface/HcalGain.h"
#include "CondFormats/HcalObjects/interface/HcalPedestalWidth.h"
#include "CondFormats/HcalObjects/interface/HcalGainWidth.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TROOT.h"
 
#include<iostream>

HcalAmplifier::HcalAmplifier(const CaloVSimParameterMap * parameters, bool addNoise):
  theDbService(0), 
  theRandGaussQ(0),
  theParameterMap(parameters),
  theStartingCapId(0), 
  addNoise_(addNoise),
  amplifierEvent_(1)
{}

void HcalAmplifier::setAmplifierEvent(int iEvent){
  
  //------------------------------------------------------
  // Sets the event that the amplifier should read in from
  // the ROOT file
  //
  // Called once per event
  //------------------------------------------------------
  
  amplifierEvent_ = iEvent;
  
}

void HcalAmplifier::setFileNames(std::string hbFile, 
				std::string heFile, 
				std::string hoFile, 
				std::string hfFile){

  //------------------------------------------------------
  // Tells the amplifier where to look for files
  //
  // Called once per job
  //------------------------------------------------------

  hbFile_ = hbFile;
  heFile_ = heFile;
  hoFile_ = hoFile;
  hfFile_ = hfFile;

}

void HcalAmplifier::readTree(){

  //------------------------------------------------------
  //  Fills two arrays
  //    1. (floats) zero-bias info for working cells 
  //    2. (bools)  says whether or not given cell was 
  //                working during CRUZET
  //
  //  Called once per event
  //------------------------------------------------------

  const int maxNumChannels   = 2592; // Number of channels   in HB/HE (more than HF)
  const int maxNumTimeSlices = 10;   // Number of timeslices in HB/HE (more than HF)

  for (int iSubdet = 1; iSubdet < 5; iSubdet++){
    for (int tempIeta = 0; tempIeta < 83; tempIeta++){
      for (int tempIphi = 0; tempIphi < 73; tempIphi++){
	for (int tempDepth = 0; tempDepth < 4; tempDepth++){
	  haveDataForThisChannel_[iSubdet][tempIeta][tempIphi][tempDepth] = false;
	}}}}
  
  float h_fC   [maxNumChannels][maxNumTimeSlices];
  int   h_ieta [maxNumChannels];
  int   h_iphi [maxNumChannels];
  int   h_depth[maxNumChannels];  
  int   nhit;
  
  float tempNoise;
  int tempIeta, tempIphi, tempDepth,  nEvent;

  //------------------------------------------------------
  // Get trees
  //------------------------------------------------------

  TChain chainArray[5];
  char fileNameArray[5][200];
  
  TString hbTString = TString(hbFile_); sprintf(fileNameArray[1],"%s/dqmtree",hbTString.Data());
  TString heTString = TString(heFile_);	sprintf(fileNameArray[2],"%s/dqmtree",heTString.Data());
  TString hoTString = TString(hoFile_);	sprintf(fileNameArray[3],"%s/dqmtree",hoTString.Data());
  TString hfTString = TString(hfFile_);	sprintf(fileNameArray[4],"%s/dqmtree",hfTString.Data());
			
  for (int iSubdet = 1; iSubdet < 5; iSubdet++)
    chainArray[iSubdet].Add(fileNameArray[iSubdet]);
  
  //------------------------------------------------------
  // Loop over the tree for a given event
  //------------------------------------------------------

  nEvent = amplifierEvent_;

  setAmplifierEvent(nEvent+1); 

  for (int iSubdet = 1; iSubdet < 5; iSubdet++){

    //------------------------------------------------------
    // Skip the HE, it has too many missing channels
    //------------------------------------------------------

    if (iSubdet == 2) continue; 

    //------------------------------------------------------
    // Get relevant info from the root trees
    //------------------------------------------------------

    chainArray[iSubdet].SetBranchAddress("h_fC"   ,  h_fC    );
    chainArray[iSubdet].SetBranchAddress("h_iphi" ,  h_iphi  );
    chainArray[iSubdet].SetBranchAddress("h_ieta" ,  h_ieta  );
    chainArray[iSubdet].SetBranchAddress("h_depth",  h_depth );
    chainArray[iSubdet].SetBranchAddress("nhit"   , &nhit    );

    //------------------------------------------------------
    // Get this event
    //------------------------------------------------------
    
    chainArray[iSubdet].GetEntry(nEvent);
    
    //------------------------------------------------------
    // Finally, loop over all channels and get fC for all
    // time slices
    //------------------------------------------------------

    for (int iCHN = 0; iCHN < nhit; iCHN++){
	
      tempIeta  =  h_ieta [iCHN];
      tempIeta  += 41;
      tempIphi  =  h_iphi [iCHN];
      tempDepth =  h_depth[iCHN];

      haveDataForThisChannel_[iSubdet][tempIeta][tempIphi][tempDepth] = true;
      
      for (int iTS = 0; iTS < 10; iTS++){
	
	tempNoise =  h_fC[iCHN][iTS];
	noiseArray_[iSubdet][tempIeta][tempIphi][tempDepth][iTS] = tempNoise;
	
      } 
    }
  } 
} 

void HcalAmplifier::amplifyCRUZET(CaloSamples & frame) const {

  //------------------------------------------------------
  // Get the channel info from the frame id
  //------------------------------------------------------

  HcalDetId tempHcalDetId = HcalDetId(frame.id());

  int tempIeta  = (int) tempHcalDetId.ieta() + 41;
  int tempIphi  = (int) tempHcalDetId.iphi();
  int subdetId  = (int) tempHcalDetId.subdet();
  int tempDepth = (int) tempHcalDetId.depth();

  //------------------------------------------------------
  // IF we don't have info on this channel
  //   (this combo of ieta, iphi, and depth),
  // THEN pass this channel on to the normal 
  //    HcalAmplifier::amplify function
  //------------------------------------------------------

  if (!haveDataForThisChannel_[subdetId][tempIeta][tempIphi][tempDepth]){
    amplify(frame);
    LogDebug("HcalAmplifier") << frame;
    return;
  }

  //------------------------------------------------------
  // Fill the time slice
  //------------------------------------------------------

  for (int tbin = 0 ; tbin < frame.size(); ++tbin)
    frame[tbin] += noiseArray_[subdetId][tempIeta][tempIphi][tempDepth][tbin];
  
  LogDebug("HcalAmplifier") << frame;
}

//------------------------------------------------------
// This is the normal HcalAmplifier::amplify function
//------------------------------------------------------

void HcalAmplifier::setRandomEngine(CLHEP::HepRandomEngine & engine)
{
  theRandGaussQ = new CLHEP::RandGaussQ(engine);
}

void HcalAmplifier::amplify(CaloSamples & frame) const {
  const CaloSimParameters & parameters = theParameterMap->simParameters(frame.id());
  assert(theDbService != 0);
  HcalGenericDetId hcalGenDetId(frame.id());
  const HcalPedestal* peds = theDbService->getPedestal(hcalGenDetId);
  const HcalPedestalWidth* pwidths = theDbService->getPedestalWidth(hcalGenDetId);
  if (!peds || !pwidths )
  {
    edm::LogError("HcalAmplifier") << "Could not fetch HCAL conditions for channel " << hcalGenDetId;
  }

  double gauss [32]; //big enough
  double noise [32]; //big enough
  double fCperPE = parameters.photoelectronsToAnalog(frame.id());

  for (int i = 0; i < frame.size(); i++) gauss[i] = theRandGaussQ->fire(0., 1.);
  pwidths->makeNoise (frame.size(), gauss, noise);
  for(int tbin = 0; tbin < frame.size(); ++tbin) {
    int capId = (theStartingCapId + tbin)%4;
    double pedestal = peds->getValue (capId);
    if(addNoise_) {
      pedestal += noise [tbin];
    }
    frame[tbin] *= fCperPE;
    frame[tbin] += pedestal;
  }
  LogDebug("HcalAmplifier") << frame;
}



