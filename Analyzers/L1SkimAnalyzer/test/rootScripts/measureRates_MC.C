const int MAXNGCTHT = 12;
const int nJobsPerThreshold = 10;
const int nEventsPerJob = 10000;
const int nPtBins = 13;

const int nL1HtThresholds = 11;
const int l1HtThresholdValues[nL1HtThresholds+1] = {10,15,20,30,40,50, 75, 100, 200, 300, 400, 500};

void measureRates_MC(int tempMinPtBin, int tempMaxPtBin, 
		     int tempMinL1JetThreshold,
		     int tempMaxL1JetThreshold ){

  const int minPtBin = tempMinPtBin;
  const int maxPtBin = tempMaxPtBin;
  
  const int minL1JetThreshold = tempMinL1JetThreshold;
  const int maxL1JetThreshold = tempMaxL1JetThreshold;    

  //--------------------------------------------------
  // Get file and tree
  //--------------------------------------------------
  
  TChain *tree[maxPtBin+1][maxL1JetThreshold+1];

  char tempFileName[500];

  for (int iPtBin = minPtBin; iPtBin <= maxPtBin; iPtBin++){
    for (int iL1JetThreshold = minL1JetThreshold; iL1JetThreshold <= maxL1JetThreshold; iL1JetThreshold++){
      
      tree[iPtBin][iL1JetThreshold] = new TChain("l1SkimTree","");
      
      for (int iJob = 1; iJob <= nJobsPerThreshold; iJob++){

	sprintf(tempFileName,
		//"/uscms/home/eberry/3DayLifetime/tmp/Analyzer/output/L1SkimAnalyzerOutput_L1EmulatorOnMC_PtBin%d_%dGeV_job%d_%devents_withCorJets.root",
		"/uscms/home/eberry/3DayLifetime/tmp/Analyzer/output/
		iPtBin,
		iL1JetThreshold,
		iJob,
		nEventsPerJob);

	cout << tempFileName << endl;
	
	tree[iPtBin][iL1JetThreshold] -> Add(tempFileName);
	
      }
    }
  }

  cout << "Done adding files" << endl;

  //--------------------------------------------------
  // Declare useful values and branch addresses
  //--------------------------------------------------

  int nEvents;
  int triggerCount[maxPtBin+1][maxL1JetThreshold+1][nL1HtThresholds+1];
  int event       [maxPtBin+1][maxL1JetThreshold+1];  
  int l1_HTT100   [maxPtBin+1][maxL1JetThreshold+1];
  int l1_HTT200   [maxPtBin+1][maxL1JetThreshold+1];
  int l1_HTT300   [maxPtBin+1][maxL1JetThreshold+1];
  int l1_HTT400   [maxPtBin+1][maxL1JetThreshold+1];
  int l1_HTT500   [maxPtBin+1][maxL1JetThreshold+1];
  float l1_gctHT  [maxPtBin+1][maxL1JetThreshold+1][MAXNGCTHT];
  
  float gctHT;

  FILE *rateFile[maxPtBin+1];
  char rateFileName[100];

  //--------------------------------------------------
  // Loop over pt bins
  //--------------------------------------------------

  for (int iPtBin = minPtBin; iPtBin <= maxPtBin; iPtBin++){

    //--------------------------------------------------
    // Open rate file for writing
    //--------------------------------------------------
    
    sprintf(rateFileName,"txt/rateFile_MC_PtBin%d_new.txt",iPtBin);
    
    rateFile[iPtBin] = fopen(rateFileName,"w");
    
    fprintf(rateFile[iPtBin],"L1HtThresholds ");
    for (int iL1HtThreshold = 0; iL1HtThreshold < nL1HtThresholds; iL1HtThreshold++)
      fprintf(rateFile[iPtBin],"%d ",l1HtThresholdValues[iL1HtThreshold]);
    fprintf(rateFile[iPtBin],"\n");
    

    //--------------------------------------------------
    // Loop over thresholds
    //--------------------------------------------------

    for (int iL1JetThreshold = minL1JetThreshold; iL1JetThreshold <= maxL1JetThreshold; iL1JetThreshold++){
 
      event[iPtBin][iL1JetThreshold] = 0;

      for (int i = 0; i < nL1HtThresholds; i++){
	triggerCount[iPtBin][iL1JetThreshold][i] = 0;
      }
      
      //--------------------------------------------------
      // Set branch addresses
      //--------------------------------------------------
      
      tree[iPtBin][iL1JetThreshold] -> SetBranchAddress("l1_HTT100", &l1_HTT100[iPtBin][iL1JetThreshold]);
      tree[iPtBin][iL1JetThreshold] -> SetBranchAddress("l1_HTT200", &l1_HTT200[iPtBin][iL1JetThreshold]);
      tree[iPtBin][iL1JetThreshold] -> SetBranchAddress("l1_HTT300", &l1_HTT300[iPtBin][iL1JetThreshold]);
      tree[iPtBin][iL1JetThreshold] -> SetBranchAddress("l1_HTT400", &l1_HTT400[iPtBin][iL1JetThreshold]);
      tree[iPtBin][iL1JetThreshold] -> SetBranchAddress("l1_HTT500", &l1_HTT500[iPtBin][iL1JetThreshold]);
      tree[iPtBin][iL1JetThreshold] -> SetBranchAddress("gctHT"    , &l1_gctHT  [iPtBin][iL1JetThreshold]);

      //--------------------------------------------------
      // Loop over events
      //--------------------------------------------------
      
      nEvents = tree[iPtBin][iL1JetThreshold] -> GetEntries();
      cout << nEvents << endl;
      
      for (int iEvent = 1; iEvent <= nEvents; iEvent ++){      
	
	event[iPtBin][iL1JetThreshold]++;
      
	tree[iPtBin][iL1JetThreshold] -> GetEntry(iEvent);
	
	gctHT = l1_gctHT[iPtBin][iL1JetThreshold][1] ;//* (1.0 / gctLSB);

	if (l1_HTT500[iPtBin][iL1JetThreshold] == 1 && gctHT < 100 ||
	    l1_HTT500[iPtBin][iL1JetThreshold] == 0 && gctHT > 100 ) cout << "ERROR" << endl;
				

	if (iEvent%10000==0) cout << "Pt Bin = " << iPtBin << ", Threshold = " << iL1JetThreshold << ", Processing Event " << iEvent << endl;

	//--------------------------------------------------
	// Make counts
	//--------------------------------------------------
	
  	for (int iL1HtThreshold = 0; iL1HtThreshold < nL1HtThresholds; iL1HtThreshold++)
	  if (gctHT > l1HtThresholdValues[iL1HtThreshold]) triggerCount[iPtBin][iL1JetThreshold][iL1HtThreshold]++;
	
      }
      
      //--------------------------------------------------
      // Print out the results
      //--------------------------------------------------
      
      printf("L1 jet threshold = %d, We looked at %d events\n",iL1JetThreshold,event[iPtBin][iL1JetThreshold]);
      for (int iL1HtThreshold = 0; iL1HtThreshold < nL1HtThresholds; iL1HtThreshold++){
	printf("     %1.5f percent (%d out of %d events) passed l1_HTT%d\n", 
	       100.0 * (double) triggerCount[iPtBin][iL1JetThreshold][iL1HtThreshold] / (double) event[iPtBin][iL1JetThreshold], 
	       triggerCount[iPtBin][iL1JetThreshold][iL1HtThreshold], 
	       event[iPtBin][iL1JetThreshold],l1HtThresholdValues[iL1HtThreshold]);
      }
      
      //--------------------------------------------------
      // Print to the file
      //--------------------------------------------------
      
      fprintf(rateFile[iPtBin],"%d ",iL1JetThreshold);

      for (int iL1HtThreshold = 0; iL1HtThreshold < nL1HtThresholds; iL1HtThreshold++)
	fprintf(rateFile[iPtBin],"%d ",triggerCount[iPtBin][iL1JetThreshold][iL1HtThreshold]);
      
      fprintf(rateFile[iPtBin],"%d\n",event[iPtBin][iL1JetThreshold]);

    } // End of the loop over l1 jet thresholds 
    
    fclose(rateFile[iPtBin]);

  } // End loop over Pt Bins

}
