#ifndef Producers_CaloTowersFromTrigPrimsCreator_CaloTowersFromTrigPrimsAnalyzerTree_h
#define Producers_CaloTowersFromTrigPrimsCreator_CaloTowersFromTrigPrimsAnalyzerTree_h

class CaloTowersFromTrigPrimsAnalyzerTree {
  
 public:
   
  enum { MAXNTPG  = 10000 };
  enum { MAXNCT   = 11000 };
  enum { MAXNRJET = 500   };
  enum { MAXNGJET = 500   };
  enum { MAXNJCT  = 50    };
  
  CaloTowersFromTrigPrimsAnalyzerTree(); 
  virtual ~CaloTowersFromTrigPrimsAnalyzerTree();
  
  int   run;		  
  int   event;  	  
  int   nct;		  
  int   ntpg;

  int   nrjet;	   
  int   ngjet;
  int   nCaloJet;  
  int   nTPGJet;   
  int   nCaloCorJet;
  int   nTPGCorJet;
  
  int   tpg_isHcal   [MAXNTPG];
  int   tpg_isEcal   [MAXNTPG];
  int   tpg_isHF     [MAXNTPG];
  int   tpg_isEB     [MAXNTPG];
  int   tpg_isEE     [MAXNTPG];
  
  int   tpg_ieta     [MAXNTPG];
  int   tpg_iphi     [MAXNTPG];
  
  float tpg_energy   [MAXNTPG];
  float tpg_et       [MAXNTPG];
  float tpg_meanEta  [MAXNTPG];
  
  int   ct_ieta  [MAXNCT];
  int   ct_iphi  [MAXNCT];
  int   ct_nhcon [MAXNCT];
  int   ct_necon [MAXNCT];

  int   ct_isMine[MAXNCT];
                          
  float ct_eEm   [MAXNCT];
  float ct_eHad  [MAXNCT];
  float ct_eOut  [MAXNCT];

  float ct_etEm  [MAXNCT];
  float ct_etHad [MAXNCT];
  float ct_etOut [MAXNCT];

  int   rjet_nct    [MAXNRJET];
  int   rjet_isCor  [MAXNRJET];
  int   rjet_isMine [MAXNRJET];

  int   rjet_ct_ieta[MAXNRJET][MAXNJCT];
  int   rjet_ct_iphi[MAXNRJET][MAXNJCT];

  float rjet_pt     [MAXNRJET];
  float rjet_et     [MAXNRJET];
  float rjet_eta    [MAXNRJET];
  float rjet_phi    [MAXNRJET];

  float gjet_pt     [MAXNGJET];
  float gjet_et     [MAXNGJET];
  float gjet_eta    [MAXNGJET];
  float gjet_phi    [MAXNGJET];

  void init() {

    run         = -999;
    event       = -999;
    ntpg        = 0;
    nct         = 0;
    
    ngjet       = 0;
    nrjet       = 0;	   
    nCaloJet    = 0;  
    nTPGJet     = 0;   
    nCaloCorJet = 0;
    nTPGCorJet  = 0;

    for (int iTPG = 0; iTPG < MAXNTPG; iTPG++){
 
      tpg_isHcal   [iTPG] = -999;
      tpg_isEcal   [iTPG] = -999;
      tpg_isHF     [iTPG] = -999;   
      tpg_isEB     [iTPG] = -999;   
      tpg_isEE     [iTPG] = -999;   
                     
      tpg_ieta     [iTPG] = -999;
      tpg_iphi     [iTPG] = -999;
                        
      tpg_energy   [iTPG] = -999.0;
      tpg_et       [iTPG] = -999.0;
      tpg_meanEta  [iTPG] = -999.0;

    }

    for (int iCT = 0; iCT < MAXNCT; iCT++){
 
      ct_ieta  [iCT] = -999;
      ct_iphi  [iCT] = -999;
      ct_nhcon [iCT] = -999;
      ct_necon [iCT] = -999;
      
      ct_isMine[iCT] = -999;

      ct_eEm   [iCT] = -999.0;
      ct_eHad  [iCT] = -999.0;
      ct_eOut  [iCT] = -999.0;

      ct_etEm  [iCT] = -999.0;
      ct_etHad [iCT] = -999.0;
      ct_etOut [iCT] = -999.0;

    }

    for (int iRJET = 0; iRJET < MAXNRJET; iRJET++){
      rjet_isCor [iRJET] = -999;
      rjet_isMine[iRJET] = -999;
      
      rjet_pt    [iRJET] = -999.0;
      rjet_et    [iRJET] = -999.0;
      rjet_eta   [iRJET] = -999.0;
      rjet_phi   [iRJET] = -999.0;

      rjet_nct   [iRJET] = -999;

      for (int iJCT = 0; iJCT < MAXNJCT; iJCT++){
	rjet_ct_ieta [iRJET][iJCT] = -999;
	rjet_ct_iphi [iRJET][iJCT] = -999;
      }
    }

    for (int iGJET = 0; iGJET < MAXNGJET; iGJET++){

      gjet_pt  [iGJET] = -999.0;
      gjet_et  [iGJET] = -999.0;
      gjet_eta [iGJET] = -999.0;
      gjet_phi [iGJET] = -999.0;

    }
    
  }
  
  
 private:
  
  CaloTowersFromTrigPrimsAnalyzerTree(const CaloTowersFromTrigPrimsAnalyzerTree&); // stop default
  const CaloTowersFromTrigPrimsAnalyzerTree& operator=(const CaloTowersFromTrigPrimsAnalyzerTree&); // stop default

};


#endif
