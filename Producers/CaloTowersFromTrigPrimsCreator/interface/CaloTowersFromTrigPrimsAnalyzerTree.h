#ifndef Producers_CaloTowersFromTrigPrimsCreator_CaloTowersFromTrigPrimsAnalyzerTree_h
#define Producers_CaloTowersFromTrigPrimsCreator_CaloTowersFromTrigPrimsAnalyzerTree_h

class CaloTowersFromTrigPrimsAnalyzerTree {
  
 public:
   
  enum { MAXNTPG = 10000 };
  enum { MAXNCT  = 11000 };
  
  CaloTowersFromTrigPrimsAnalyzerTree(); 
  virtual ~CaloTowersFromTrigPrimsAnalyzerTree();
  
  int   run;		  
  int   event;  	  
  int   nct;		  
  int   ntpg;		  
  
  int   tpg_isHcal   [MAXNTPG];
  int   tpg_isEcal   [MAXNTPG];
  int   tpg_isHF     [MAXNTPG];
  
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

  void init() {

    run     = -999;
    event   = -999;
    ntpg    = -999;
    nct     = -999;

    for (int iTPG = 0; iTPG < MAXNTPG; iTPG++){
 
      tpg_isHcal   [iTPG] = -999;
      tpg_isEcal   [iTPG] = -999;
      tpg_isHF     [iTPG] = -999;   
                     
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

    
  }
  
  
 private:

  CaloTowersFromTrigPrimsAnalyzerTree(const CaloTowersFromTrigPrimsAnalyzerTree&); // stop default
  const CaloTowersFromTrigPrimsAnalyzerTree& operator=(const CaloTowersFromTrigPrimsAnalyzerTree&); // stop default

};


#endif
