#include <iostream>

#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH2F.h>
#include <TH1F.h>

//-------------------------------------------------
// Shouldn't change
//-------------------------------------------------

const char* TREE_NAME = "CaloTowersFromTrigPrimsAnalyzerTree";

//-------------------------------------------------
// Script running information
//-------------------------------------------------

const int MAX_N_EVENTS = -1;
const bool SAVE_PLOTS = true;  

//-------------------------------------------------
// Analyzer settings
//------------------------------------------------- 

const int MAXNCT  = 11000;
const int MAXNTPG = 10000;

//-------------------------------------------------
// Global variables
//-------------------------------------------------

TFile *m_inputFile;
TTree *m_inputTree;

TH2F  *h_cct_occupancy;
TH2F  *h_dct_occupancy;

TH1F  *h_cct_totalE;      
TH1F  *h_cct_totalEEm;    
TH1F  *h_cct_totalEHad;   

TH1F  *h_dct_totalE;      
TH1F  *h_dct_totalEEm;    
TH1F  *h_dct_totalEHad;   

TH1F  *h_tpg_totalE;     
TH1F  *h_tpg_totalEEm;   
TH1F  *h_tpg_totalEHad;  

TH1F  *h_ratio_totalE;   
TH1F  *h_ratio_totalEEm; 
TH1F  *h_ratio_totalEHad;


char INPUT_FILE_NAME  [300];
char OUTPUT_FILE_NAME [300];

//-------------------------------------------------
// Set file names
//-------------------------------------------------

bool setFileNames(){

  sprintf(INPUT_FILE_NAME ,"CaloTowersFromTrigPrimsAnalyzer.root");
  sprintf(OUTPUT_FILE_NAME,"Plots.root");
  
  return true;
}

//-------------------------------------------------
// Set plots
//-------------------------------------------------

bool setPlots(){
  
  h_cct_occupancy   = new TH2F("cct_occupancy","",83,-41.5,41.5,73,-0.5,72.5);
  h_dct_occupancy   = new TH2F("dct_occupancy","",83,-41.5,41.5,73,-0.5,72.5);

  h_cct_totalE      = new TH1F("cct_totalE"   ,"",100,0,400);
  h_cct_totalEEm    = new TH1F("cct_totalEEm" ,"",100,0,400);
  h_cct_totalEHad   = new TH1F("cct_totalEHad","",100,0,400);

  h_dct_totalE      = new TH1F("dct_totalE"   ,"",100,0,400);
  h_dct_totalEEm    = new TH1F("dct_totalEEm" ,"",100,0,400);
  h_dct_totalEHad   = new TH1F("dct_totalEHad","",100,0,400);
  		    
  h_tpg_totalE      = new TH1F("tpg_totalE"   ,"",100,0,400);
  h_tpg_totalEEm    = new TH1F("tpg_totalEEm" ,"",100,0,400);
  h_tpg_totalEHad   = new TH1F("tpg_totalEHad","",100,0,400);
  		   
  h_ratio_totalE    = new TH1F ("ratio_totalE"   ,"",100,0.0,2.0);
  h_ratio_totalEEm  = new TH1F ("ratio_totalEEm" ,"",100,0.0,2.0);
  h_ratio_totalEHad = new TH1F ("ratio_totalEHad","",100,0.0,2.0);


  return true;

}

//-------------------------------------------------
// Get Trees from ROOT files
//-------------------------------------------------

bool getTrees(){

  m_inputFile = new TFile(INPUT_FILE_NAME);
  
  m_inputTree = (TTree*) m_inputFile -> Get(TREE_NAME);
  
  if (!m_inputTree) {
    std::cout << "No tree found!" << std::endl;
    std::cout << "  I looked here: " << INPUT_FILE_NAME << std::endl;
    std::cout << "  And I used tree name: " << TREE_NAME << std::endl;
    return false;
  }
  
  else return true;
  
}

//-------------------------------------------------
// Analyze
//-------------------------------------------------

bool analysis(){

  int run, event, ntpg, nct;
  
  int   ct_ieta     [MAXNCT];
  int   ct_iphi     [MAXNCT];
  int   ct_nhcon    [MAXNCT];
  int   ct_necon    [MAXNCT];
  int   ct_isMine   [MAXNCT];
  		    
  float ct_eEm      [MAXNCT];
  float ct_eHad     [MAXNCT];
  float ct_etEm     [MAXNCT];
  float ct_etHad    [MAXNCT];
  
  int   tpg_ieta    [MAXNTPG];
  int   tpg_iphi    [MAXNTPG];

  int   tpg_isHcal  [MAXNTPG];
  int   tpg_isEcal  [MAXNTPG];
  int   tpg_isHF    [MAXNTPG];
  
  float tpg_eta     [MAXNTPG];
  float tpg_e       [MAXNTPG];
  float tpg_et      [MAXNTPG];

  //-------------------------------------------------
  // Set branch addresses
  //-------------------------------------------------

  m_inputTree -> SetBranchAddress("run"        , &run         );	
  m_inputTree -> SetBranchAddress("event"      , &event       );
  m_inputTree -> SetBranchAddress("nct"        , &nct         );	
  m_inputTree -> SetBranchAddress("ntpg"       , &ntpg        );	
  
  m_inputTree -> SetBranchAddress("ct_ieta"    ,  ct_ieta     );
  m_inputTree -> SetBranchAddress("ct_iphi"    ,  ct_iphi     );
  m_inputTree -> SetBranchAddress("ct_nhcon"   ,  ct_nhcon    ); 
  m_inputTree -> SetBranchAddress("ct_necon"   ,  ct_necon    ); 
  
  m_inputTree -> SetBranchAddress("ct_eEm"     ,  ct_eEm      );
  m_inputTree -> SetBranchAddress("ct_eHad"    ,  ct_eHad     );       
  m_inputTree -> SetBranchAddress("ct_etEm"    ,  ct_etEm     );
  m_inputTree -> SetBranchAddress("ct_etHad"   ,  ct_etHad    ); 
  m_inputTree -> SetBranchAddress("ct_isMine"  ,  ct_isMine   );
  

  m_inputTree -> SetBranchAddress("tpg_ieta"   ,  tpg_ieta    );
  m_inputTree -> SetBranchAddress("tpg_iphi"   ,  tpg_iphi    );
  
  m_inputTree -> SetBranchAddress("tpg_isHcal" ,  tpg_isHcal  );
  m_inputTree -> SetBranchAddress("tpg_isEcal" ,  tpg_isEcal  );
  m_inputTree -> SetBranchAddress("tpg_isHF"   ,  tpg_isHF    );
  
  m_inputTree -> SetBranchAddress("tpg_meanEta",  tpg_eta     );
  m_inputTree -> SetBranchAddress("tpg_energy" ,  tpg_e       );
  m_inputTree -> SetBranchAddress("tpg_et"     ,  tpg_et      );

  float eEm, eHad, etEm, etHad;
  float e, et;
  int ieta, iphi;
  
  float event_cct_eEm   ;
  float event_cct_eHad  ;
  float event_dct_eEm   ;
  float event_dct_eHad  ;
  float event_tpg_eEm  ;
  float event_tpg_eHad ;

  float event_cct_setSumEm;
  float event_cct_setSumHad;
  float event_dct_setSumEm;
  float event_dct_setSumHad;
  float event_tpg_setSumEm;
  float event_tpg_setSumHad;

  bool isHcal, isEcal, isHF;
  bool isMine;

  //-------------------------------------------------
  // Get maximum events
  //-------------------------------------------------

  int nEvents;

  int n_entries = m_inputTree -> GetEntries();

  if ( MAX_N_EVENTS < 0 ) nEvents = n_entries;
  
  else nEvents = TMath::Min( MAX_N_EVENTS, n_entries );  

  //-------------------------------------------------
  // Loop over the events
  //-------------------------------------------------

  for (int iEvent = 0; iEvent <= nEvents; iEvent++){

    //-------------------------------------------------
    // Report progress 
    //-------------------------------------------------

    if (iEvent%100 == 0 && iEvent != 0 && iEvent != nEvents)
      std::cout << "Processing Event " << iEvent << "/" << nEvents << std::endl;

    //-------------------------------------------------
    // Get tree entries
    //-------------------------------------------------

    m_inputTree -> GetEntry(iEvent);

    //-------------------------------------------------
    // Zero-out total-event values
    //-------------------------------------------------

    event_cct_eEm   = 0.0;
    event_dct_eEm   = 0.0;
    event_tpg_eEm   = 0.0;

    event_cct_eHad  = 0.0;
    event_dct_eHad  = 0.0;
    event_tpg_eHad  = 0.0;

    //-------------------------------------------------
    // Loop over CaloTowers and set single-tower values
    //-------------------------------------------------

    for (int iCT = 0; iCT < nct; iCT++){

      eEm   = ct_eEm   [iCT];
      eHad  = ct_eHad  [iCT];
      etEm  = ct_etEm  [iCT];

      etHad = ct_etHad [iCT];
      ieta  = ct_ieta  [iCT];
      iphi  = ct_iphi  [iCT];

      isMine = (ct_isMine[iCT] == 1);      

      //-------------------------------------------------
      // Is this tower mine?
      //-------------------------------------------------
      
      if (isMine){
	event_cct_eEm  += eEm;
	event_cct_eHad += eHad;      
	event_cct_setSumEm   += etEm;
	event_cct_setSumHad  += etHad;
	
	h_cct_occupancy -> Fill (ieta, iphi);
      }
      
      if (!isMine){
	event_dct_eEm  += eEm;
	event_dct_eHad += eHad;      
	event_dct_setSumEm   += etEm;
	event_dct_setSumHad  += etHad;

	h_dct_occupancy -> Fill (ieta, iphi);
      }
      
      
      
    }
    
    //-------------------------------------------------
    // Loop over TPG's
    //-------------------------------------------------

    for (int iTPG = 0; iTPG < ntpg; iTPG++){
      
      isHcal = ( tpg_isHcal[iTPG] == 1 );
      isEcal = ( tpg_isEcal[iTPG] == 1 );
      isHF   = ( tpg_isHF  [iTPG] == 1 );
      
      e  = tpg_e  [iTPG];
      et = tpg_et [iTPG];
      
      ieta  = tpg_ieta  [iTPG];
      iphi  = tpg_iphi  [iTPG];
      
      if (e != e ) continue;
      
      if (isHcal) {
	event_tpg_eHad += e;
	event_tpg_setSumHad += et;
      }
      
      if (isEcal) {
	event_tpg_eEm  += e;
	event_tpg_setSumEm  += et;
      }
    }
    
    float event_cct_e  = event_cct_eEm  + event_cct_eHad ;
    float event_dct_e  = event_dct_eEm  + event_dct_eHad ;
    float event_tpg_e  = event_tpg_eEm + event_tpg_eHad;
    
    h_cct_totalE       -> Fill (event_cct_e);
    h_cct_totalEEm     -> Fill (event_cct_eEm);		    
    h_cct_totalEHad    -> Fill (event_cct_eHad);		    
    
    h_dct_totalE       -> Fill (event_dct_e);
    h_dct_totalEEm     -> Fill (event_dct_eEm);		    
    h_dct_totalEHad    -> Fill (event_dct_eHad);		    
    
    h_tpg_totalE      -> Fill (event_tpg_e);
    h_tpg_totalEEm    -> Fill (event_tpg_eEm);		    
    h_tpg_totalEHad   -> Fill (event_tpg_eHad);		    
    		     
    h_ratio_totalE    -> Fill ( event_cct_e    / event_dct_e    );
    h_ratio_totalEEm  -> Fill ( event_cct_eEm  / event_dct_eEm  );
    h_ratio_totalEHad -> Fill ( event_cct_eHad / event_dct_eHad );    
    
    if (iEvent == nEvents) 
      std::cout << "Processed Event " << iEvent << "/" << nEvents << " ... DONE!" << std::endl;
    
  }
  
  return true;
  
}

void savePlots(){

  TFile *file = new TFile(OUTPUT_FILE_NAME, "RECREATE");
  
  file -> cd();

  h_dct_totalE      -> Write();
  h_dct_totalEEm    -> Write();
  h_dct_totalEHad   -> Write();

  h_cct_totalE      -> Write();
  h_cct_totalEEm    -> Write();
  h_cct_totalEHad   -> Write(); 
  		   
  h_tpg_totalE      -> Write();
  h_tpg_totalEEm    -> Write();
  h_tpg_totalEHad   -> Write();
  		   
  h_ratio_totalE    -> Write();
  h_ratio_totalEEm  -> Write();
  h_ratio_totalEHad -> Write();

  h_cct_occupancy   -> Write();
  h_dct_occupancy   -> Write();

  file -> Close();

  std::cout << OUTPUT_FILE_NAME << std::endl;
  
}


//-------------------------------------------------
// Main function
//-------------------------------------------------

void analyze(){

  setFileNames();

  bool plotsSet = setPlots();   if (!plotsSet) return;
  bool gotTrees = getTrees();   if (!gotTrees) return;
  bool analyzed = analysis();    if (!analyzed) return;
  
  if (SAVE_PLOTS) savePlots();
  
}
