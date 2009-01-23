
const int nJobsPerThreshold = 10;

const int nEventsPerJob = 10000;

const int eventCutOff = -1;
const int nPtBins = 13;

const int jetThreshold = 10;

const int minPtBin = 1;
const int maxPtBin = 13;

const int maxHltJets = 100;

const int nHltTriggers  = 3;
const int nL1Triggers   = 6;

const int MAXNGCTHT = 12;

const int     nHtBins      = 300;
const double  minHtBinEdge = 0.0;
const double  maxHtBinEdge = 300.0;
const double  htBinRange   = maxHtBinEdge - minHtBinEdge;

double htBins[nHtBins+1];

int hltTriggerValues[nHltTriggers+1] = {-1,50,100,150};
int l1TriggerValues [nL1Triggers+1]  = {-1,10,15,20,50,100,150};

TH1F* h_gctHT[nPtBins+1];

int getHtBinMax (int iHtBin){

  if (iHtBin < 0){
    cout << "ERROR: You are asking for the maximum of nonsensical bin: " << iHtBin << endl;
    return -999;
  }

  if (iHtBin > nHtBins) {
    cout << "ERROR: You are asking for the maximum of bin " << iHtBin << ", but there are only " << nHtBins << " HLT bins " << endl;
    return -999;
  }

  return htBins[iHtBin];

}

int getHtBinMin (int iHtBin){

  double tempBinMin, tempBinMax;

  if (iHtBin < 0){
    cout << "ERROR: You are asking for the maximum of nonsensical bin: " << iHtBin << endl;
    return -999;
  }

  if (iHtBin > nHtBins) {
    cout << "ERROR: You are asking for the maximum of bin " << iHtBin << ", but there are only " << nHtBins << " HLT bins " << endl;
    return -999;
  }
  
  return htBins[iHtBin - 1];

}

int getHtBinNum(double ht){

  double tempBinMin, tempBinMax;
  double binNumber;

  for (int iHtBin = 1; iHtBin <= nHtBins; iHtBin++){    
    
    tempBinMin = htBins[iHtBin - 1];
    tempBinMax = htBins[iHtBin];
    binNumber  = iHtBin;

    if (ht >= tempBinMin && 
	ht <  tempBinMax   ) return binNumber;

  }

  return -999;

}

int setupHtBins(){
  for (int iHtBin = 0; iHtBin <= nHtBins; iHtBin++)
    htBins[iHtBin] = minHtBinEdge + iHtBin * (htBinRange / ((double) nHtBins));
}

void measureEfficiencies_withBinning_MC(){

  setupHtBins();
    
  //--------------------------------------------------
  // Sanity check
  //--------------------------------------------------

  if (eventCutOff > nEventsPerJob*nJobsPerThreshold) {

    cout << "eventCutOff = " << eventCutOff 
	 << ", but three are only " << nEventsPerJob*nJobsPerThreshold << " events per file" << endl;
    return;
  }

  //--------------------------------------------------
  // Get file and tree
  //--------------------------------------------------
  
  TChain tree[maxPtBin+1];

  char tempFileName[300];
  char histName[300];

  for (int iPtBin = minPtBin; iPtBin <= maxPtBin; iPtBin++){
    
    sprintf(histName,"gctHT_ptBin%d",iPtBin);

    h_gctHT[iPtBin] = new TH1F(histName,"",1000,0,1000);

    tree[iPtBin] = new TChain("l1SkimTree","");

    for (int iJob = 1; iJob <= nJobsPerThreshold; iJob++){

      sprintf(tempFileName,
	      //"/uscms/home/eberry/data/L1AndHLTOnMC/L1SkimAnalyzerOutput_L1EmulatorOnMC_PtBin%d_%dGeV_job%d_%devents_withHLT.root",
	      "/uscms/home/eberry/data/L1AndHLTOnMC_AllJetThresholds/L1SkimAnalyzerOutput_L1EmulatorOnMC_PtBin%d_%dGeV_job%d_%devents_withHLT.root",
	      iPtBin,
	      jetThreshold,
	      iJob,
	      nEventsPerJob);

      
      tree[iPtBin].Add(tempFileName);
      
    }
  }

  //--------------------------------------------------
  // Declare useful values and branch addresses
  //--------------------------------------------------

  int nEventsInPtBin [maxPtBin+1];
  int passHLT        [maxPtBin+1][nHtBins+1][nHltTriggers+1];
  int passL1AndHLT   [maxPtBin+1][nHtBins+1][nL1Triggers+1][nHltTriggers+1];
  float l1_gctHT     [maxPtBin+1][MAXNGCTHT];
  int event          [maxPtBin+1];  

  int nHltJets       [maxPtBin+1];
  float hltJet_pt    [maxPtBin+1][maxHltJets];
  int  nEventsInPtHtBinCounter[maxPtBin+1][nHtBins+1];
  FILE *rateFile     [maxPtBin+1][nHtBins+1];
  char rateFileName  [500];
  float hltJetPt, gctHT, hltHT;
  bool pass;
  int thisHtBinNum;

  //--------------------------------------------------
  // Loop over pt bins
  //--------------------------------------------------
  
  for (int iPtBin = minPtBin; iPtBin <= maxPtBin; iPtBin++){
    
    //--------------------------------------------------
    // Set branch addresses
    //--------------------------------------------------
    
    tree[iPtBin].SetBranchAddress("gctHT"    , l1_gctHT[iPtBin]);
    tree[iPtBin].SetBranchAddress("nHLTJetCands" , &nHltJets [iPtBin]);
    tree[iPtBin].SetBranchAddress("hltJet_pt",  hltJet_pt[iPtBin]);
    
    //--------------------------------------------------
    // Zero-out arrays
    //--------------------------------------------------
    
    for (int iHtBin = 1; iHtBin <= nHtBins; iHtBin++){
      nEventsInPtHtBinCounter[iPtBin][iHtBin] = 0;
      
      for (int iHltTrigger = 1; iHltTrigger <= nHltTriggers; iHltTrigger++){
	passHLT[iPtBin][iHtBin][iHltTrigger] = 0;
	for (int iL1Trigger = 1; iL1Trigger <= nL1Triggers; iL1Trigger++){
	  passL1AndHLT[iPtBin][iHtBin][iL1Trigger][iHltTrigger] = 0;
	}
      }
    }
     
    //--------------------------------------------------
    // Loop over events
    //--------------------------------------------------
      
    if (eventCutOff < 0)
      nEventsInPtBin[iPtBin] = tree[iPtBin].GetEntries();
    if (eventCutOff > 0)
      nEventsInPtBin[iPtBin] = eventCutOff;
    
    cout << "PtBin = " << iPtBin << ", nEvents = " << nEventsInPtBin[iPtBin] << endl;
    
    for (int iEvent = 1; iEvent <= nEventsInPtBin[iPtBin]; iEvent ++){      
      
      tree[iPtBin].GetEntry(iEvent);

      if (iEvent%10000==0) cout << "Pt Bin = " << iPtBin << ", Processing Event " << iEvent << endl;
      
      //--------------------------------------------------
      // Calculate GCT HT
      //--------------------------------------------------

      gctHT = l1_gctHT[iPtBin][1];
      
      //--------------------------------------------------
      // Calculate HLT HT
      //--------------------------------------------------

      hltHT = 0.0;
      
      for (int iHltJet = 0; iHltJet < nHltJets[iPtBin]; iHltJet++){
	hltJetPt = hltJet_pt[iPtBin][iHltJet];
	if (hltJetPt > jetThreshold) hltHT += hltJetPt;
      }
      
      if (hltHT > 20.0){
	h_gctHT[iPtBin] -> Fill(gctHT);
      }

      //--------------------------------------------------
      // Get the HT bin and make sure it makes sense
      //--------------------------------------------------
      
      thisHtBinNum = getHtBinNum(hltHT);
      
      if (thisHtBinNum == -999) continue;

      if ( hltHT < getHtBinMin(thisHtBinNum) || hltHT > getHtBinMax(thisHtBinNum)) {
	cout << "ERROR: " << endl;
	cout << "    ht            = " << hltHT << endl;
	cout << "    thisHtBinNum  = " << thisHtBinNum << endl;
	cout << "    bin min       = " << getHtBinMin(thisHtBinNum) << endl;
	cout << "    bin max       = " << getHtBinMax(thisHtBinNum) << endl;
      }
      
      nEventsInPtHtBinCounter[iPtBin][thisHtBinNum]++;
      
      //--------------------------------------------------
      // Loop over HLT triggers
      //--------------------------------------------------
      
      for (int iHltTrigger = 1; iHltTrigger <= nHltTriggers; iHltTrigger++){
	
	if (hltHT > hltTriggerValues[iHltTrigger]){
	  passHLT[iPtBin][thisHtBinNum][iHltTrigger]++;

	  for (int iL1Trigger = 1; iL1Trigger <= nL1Triggers; iL1Trigger++){	   
	    if (gctHT >= l1TriggerValues[iL1Trigger]) {
	      passL1AndHLT[iPtBin][thisHtBinNum][iL1Trigger][iHltTrigger]++;
	    }	  
	  }
	}
      } // end loop over HLT triggers
    }// end loop over events
    
    //--------------------------------------------------
    // Loop over ht bins for writing out files
    //--------------------------------------------------
    
    for (int iHtBin = 1; iHtBin <= nHtBins; iHtBin++){
      
      //--------------------------------------------------
      // Open rate file for writing
      //--------------------------------------------------
      
      sprintf(rateFileName,
	      "/uscms/home/eberry/CMSSW_2_1_11/src/Analyzers/L1SkimAnalyzer/test/txt/effiFile_MC_PtBin%d_HtBin%d_threshold%dGeV.txt",
	      iPtBin,iHtBin,jetThreshold);
      
      rateFile[iPtBin][iHtBin] = fopen(rateFileName,"w");
      
      //--------------------------------------------------
      // Print out the results
      //--------------------------------------------------
      
      fprintf(rateFile[iPtBin][iHtBin],"nEventsInPtBin %d\n",nEventsInPtBin[iPtBin]);
      fprintf(rateFile[iPtBin][iHtBin],"nEventsInPtHtBin %d\n",nEventsInPtHtBinCounter[iPtBin][iHtBin]);
      fprintf(rateFile[iPtBin][iHtBin],"htBinRange %1.1f %1.1f\n",getHtBinMin(iHtBin),getHtBinMax(iHtBin));
      
      fprintf(rateFile[iPtBin][iHtBin],"jetThreshold %d\n",jetThreshold);
      fprintf(rateFile[iPtBin][iHtBin],"L1Triggers ");
      
      for (int iL1Trigger = 1; iL1Trigger <= nL1Triggers; iL1Trigger++){
	fprintf(rateFile[iPtBin][iHtBin],"%d ",l1TriggerValues[iL1Trigger]);
      }
      fprintf(rateFile[iPtBin][iHtBin],"\n");
      
      fprintf(rateFile[iPtBin][iHtBin],"hltTriggersVsPassBothAndNPassHLT\n");
      
      for (int iHltTrigger = 1; iHltTrigger <= nHltTriggers; iHltTrigger++){
	
	if (passL1AndHLT[iPtBin][iHtBin][5][iHltTrigger] != 
	    passL1AndHLT[iPtBin][iHtBin][6][iHltTrigger] 	  ){
	  cout << "DIFF: PtBin = " << iPtBin << ", HtBin = " << iHtBin << ", iHltTrigger = " << iHltTrigger << endl;
	  cout << "L1 = 5: " << passL1AndHLT[iPtBin][iHtBin][5][iHltTrigger] << endl;
	  cout << "L1 = 6: " << passL1AndHLT[iPtBin][iHtBin][6][iHltTrigger] << endl;
	}

	
	fprintf(rateFile[iPtBin][iHtBin],"%d ", hltTriggerValues[iHltTrigger]);
	
	for (int iL1Trigger = 1; iL1Trigger <= nL1Triggers; iL1Trigger++){
	  fprintf(rateFile[iPtBin][iHtBin],"%d ",passL1AndHLT[iPtBin][iHtBin][iL1Trigger][iHltTrigger]);
	}
	
	fprintf(rateFile[iPtBin][iHtBin],"%d\n", passHLT[iPtBin][iHtBin][iHltTrigger]);
	
      }      
      
      //--------------------------------------------------
      // Write to a text file for each pt bin
      //--------------------------------------------------
      
      fclose(rateFile[iPtBin][iHtBin]);
    } // End loop over ht bins
  } // End loop over pt bins

  TFile *htFile = new TFile("htDistros.root","RECREATE");
  htFile -> cd();
  for (int iPtBin = minPtBin; iPtBin <= maxPtBin; iPtBin++){

    h_gctHT[iPtBin] -> Write(h_gctHT[iPtBin] -> GetName());

  }

} // End program
