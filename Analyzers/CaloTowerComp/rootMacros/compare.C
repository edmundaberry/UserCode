#include <iostream>

#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH2F.h>
#include <TH1F.h>

//-------------------------------------------------
// User set information
//-------------------------------------------------

const int MAX_N_EVENTS = -1;

const bool SAVE_PLOTS = true;

const char* OUT_FILE_NAME = "plots/plots.root";

const char* CREATED_FILE_NAME = 
  "/uscms/home/eberry/TowerAlgo/CMSSW_2_2_7/src/Producers/CaloTowersFromTrigPrimsCreator/data/CaloTowerInfo_Created.root";

const char* DEFAULT_FILE_NAME = 
  "/uscms/home/eberry/TowerAlgo/CMSSW_2_2_7/src/Producers/CaloTowersFromTrigPrimsCreator/data/CaloTowerInfo_Default.root";

const char* TREE_NAME = "caloTowerCompTree";

//-------------------------------------------------
// Analyzer settings
//------------------------------------------------- 

const int MAXNCT = 6000;

//-------------------------------------------------
// Global variables
//-------------------------------------------------

TFile *m_createdFile;
TFile *m_defaultFile;

TTree *m_createdTree;
TTree *m_defaultTree;

TH2F  *h_occupancy_default;
TH2F  *h_occupancy_created;

TH1F  *h_eEm_default;  TH1F  *h_eEm_created; 
TH1F  *h_eHad_default; TH1F  *h_eHad_created;
TH1F  *h_eOut_default; TH1F  *h_eOut_created;

//-------------------------------------------------
// Set plots
//-------------------------------------------------

bool setPlots(){

  h_occupancy_default = new TH2F("occupancy_default","",83,-41.5,41.5,73,-0.5,72.5);
  h_occupancy_created = new TH2F("occupancy_created","",83,-41.5,41.5,73,-0.5,72.5);
                      
  h_eEm_default       = new TH1F("eEm_default" ,"",100,0,200);
  h_eHad_default      = new TH1F("eHad_default","",100,0,200);
  h_eOut_default      = new TH1F("eOut_default","",100,0,200);
  		      			       
  h_eEm_created       = new TH1F("eEm_created" ,"",100,0,200);
  h_eHad_created      = new TH1F("eHad_created","",100,0,200);
  h_eOut_created      = new TH1F("eOut_created","",100,0,200);

  return true;

}

//-------------------------------------------------
// Get Trees from ROOT files
//-------------------------------------------------

bool getTrees(){

  m_createdFile = new TFile(CREATED_FILE_NAME);
  m_defaultFile = new TFile(DEFAULT_FILE_NAME);
  
  m_createdTree = (TTree*) m_createdFile -> Get(TREE_NAME);
  m_defaultTree = (TTree*) m_defaultFile -> Get(TREE_NAME);
  
  if (!m_createdTree) {
    std::cout << "No CaloTowerComp analyzer output for CREATED CaloTowers" << std::endl;
    std::cout << "  I looked here: " << CREATED_FILE_NAME << endl;
    return false;
  }
  
  else if (!m_defaultTree) {
    std::cout << "No CaloTowerComp analyzer output for DEFAULT CaloTowers" << std::endl;
    std::cout << "  I looked here: " << DEFAULT_FILE_NAME << endl;
    return false;
  }
  
  else return true;
  
}

//-------------------------------------------------
// Analyze
//-------------------------------------------------

bool analyze(){

  //-------------------------------------------------
  // Declare arrays and integers for memory
  //-------------------------------------------------

  int   created_run;		    int   default_run;		  
  int   created_event;  	    int   default_event;  	  
  int   created_nct;		    int   default_nct;		  
                          	                            	  
  int   ct_created_ieta  [MAXNCT];  int   ct_default_ieta  [MAXNCT];
  int   ct_created_iphi  [MAXNCT];  int   ct_default_iphi  [MAXNCT];
  int   ct_created_nhcon [MAXNCT];  int   ct_default_nhcon [MAXNCT];
  int   ct_created_necon [MAXNCT];  int   ct_default_necon [MAXNCT];
          			            			  
  float ct_created_eEm   [MAXNCT];  float ct_default_eEm   [MAXNCT];
  float ct_created_eHad  [MAXNCT];  float ct_default_eHad  [MAXNCT];
  float ct_created_eOut  [MAXNCT];  float ct_default_eOut  [MAXNCT];

  //-------------------------------------------------
  // Set branch addresses
  //-------------------------------------------------

  m_createdTree -> SetBranchAddress("run",      &created_run    );	
  m_createdTree -> SetBranchAddress("event",    &created_event  );
  m_createdTree -> SetBranchAddress("nct",      &created_nct    );	
					     
  m_createdTree -> SetBranchAddress("ct_ieta",  ct_created_ieta );
  m_createdTree -> SetBranchAddress("ct_iphi",  ct_created_iphi );
  m_createdTree -> SetBranchAddress("ct_nhcon", ct_created_nhcon); 
  m_createdTree -> SetBranchAddress("ct_necon", ct_created_necon); 
					       		     
  m_createdTree -> SetBranchAddress("ct_eEm",   ct_created_eEm  );
  m_createdTree -> SetBranchAddress("ct_eHad",  ct_created_eHad ); 
  m_createdTree -> SetBranchAddress("ct_eOut",  ct_created_eOut );

  m_defaultTree -> SetBranchAddress("run",      &default_run    );	
  m_defaultTree -> SetBranchAddress("event",    &default_event  );
  m_defaultTree -> SetBranchAddress("nct",      &default_nct    );	
					     
  m_defaultTree -> SetBranchAddress("ct_ieta",  ct_default_ieta );
  m_defaultTree -> SetBranchAddress("ct_iphi",  ct_default_iphi );
  m_defaultTree -> SetBranchAddress("ct_nhcon", ct_default_nhcon); 
  m_defaultTree -> SetBranchAddress("ct_necon", ct_default_necon); 
					       		     
  m_defaultTree -> SetBranchAddress("ct_eEm",   ct_default_eEm  );
  m_defaultTree -> SetBranchAddress("ct_eHad",  ct_default_eHad ); 
  m_defaultTree -> SetBranchAddress("ct_eOut",  ct_default_eOut );

  //-------------------------------------------------
  // Get maximum events
  //-------------------------------------------------

  int nEvents;
  
  int n_created_entries = m_createdTree -> GetEntries();
  int n_default_entries = m_defaultTree -> GetEntries();
  
  if (n_created_entries != n_default_entries){
    std::cout << "Differing numbers of events for these ROOT files" << std::endl;
    std::cout << "  Created ROOT file is: " << CREATED_FILE_NAME << std::endl;
    std::cout << "  And it has " << n_created_entries << " events" << std::endl;
    std::cout << "  Default ROOT file is: " << DEFAULT_FILE_NAME << std::endl;
    std::cout << "  And it has " << n_default_entries << " events" << std::endl;
    return false;
  }

  else if ( MAX_N_EVENTS < 0 ) nEvents = n_created_entries;
  
  else nEvents = TMath::Min( MAX_N_EVENTS, n_created_entries );  

  //-------------------------------------------------
  // Loop over the events
  //-------------------------------------------------

  for (int iEvent = 1; iEvent <= nEvents; iEvent++){

    //-------------------------------------------------
    // Report progress 
    //-------------------------------------------------

    if (iEvent%10 == 0) std::cout << "Processing Event " << iEvent << "/" << nEvents << std::endl;

    //-------------------------------------------------
    // Get tree entries
    //-------------------------------------------------

    m_createdTree -> GetEntry(iEvent);
    m_defaultTree -> GetEntry(iEvent);
    
    //-------------------------------------------------
    // Make sure these events are the same
    //-------------------------------------------------

    if (default_run   != created_run ||
	default_event != created_event ){
      std::cout << "Run/Event mismatch!" << std::endl;
      std::cout << "  This is event number: " << iEvent << std::endl;
      std::cout << "  The default ROOT file is :" << DEFAULT_FILE_NAME << std::endl;
      std::cout << "  The created ROOT file is :" << CREATED_FILE_NAME << std::endl;
      std::cout << "  On this event, the default run/event is: " << default_run << "/" << default_event << std::endl;
      std::cout << "  On this event, the created run/event is: " << created_run << "/" << created_event << std::endl;
      return false;
    }

    //-------------------------------------------------
    // Loop over created plots
    //-------------------------------------------------
    
    for (int ict = 0; ict < created_nct; ict++){

      h_occupancy_created -> Fill(ct_created_ieta[ict], ct_created_iphi[ict]);
      h_eEm_created       -> Fill(ct_created_eEm [ict]);
      h_eHad_created      -> Fill(ct_created_eHad[ict]);
      h_eOut_created      -> Fill(ct_created_eOut[ict]);

    }
   
    //-------------------------------------------------
    // Loop over default plots
    //-------------------------------------------------
    
    for (int ict = 0; ict < default_nct; ict++){

      h_occupancy_default -> Fill(ct_default_ieta[ict], ct_default_iphi[ict]);
      h_eEm_default       -> Fill(ct_default_eEm [ict]);
      h_eHad_default      -> Fill(ct_default_eHad[ict]);
      h_eOut_default      -> Fill(ct_default_eOut[ict]);

    }
  }

  return true;

}

void savePlots(){

  TFile *file = new TFile(OUT_FILE_NAME, "RECREATE");
  
  file -> cd();

  h_occupancy_default -> Write();
  h_occupancy_created -> Write();
                      
  h_eEm_default       -> Write();
  h_eHad_default      -> Write();
  h_eOut_default      -> Write();
  		      
  h_eEm_created       -> Write();
  h_eHad_created      -> Write();
  h_eOut_created      -> Write();

  file -> Close();

}


//-------------------------------------------------
// Main function
//-------------------------------------------------

void compare(){
  
  bool plotsSet = setPlots();   if (!plotsSet) return;
  bool gotTrees = getTrees();   if (!gotTrees) return;
  bool analyzed = analyze();    if (!analyzed) return;

  if (SAVE_PLOTS) savePlots();

}
