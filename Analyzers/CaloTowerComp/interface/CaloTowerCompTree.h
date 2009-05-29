#ifndef Analyzers_CaloTowerComp_CaloTowerCompTree_h
#define Analyzers_CaloTowerComp_CaloTowerCompTree_h

class CaloTowerCompTree {
  
 public:
   
  enum { MAXNCT = 6000 };
  
  CaloTowerCompTree(); 
  virtual ~CaloTowerCompTree();
  
  int   run;		  
  int   event;  	  
  int   nct;		  
                          
  int   ct_ieta  [MAXNCT];
  int   ct_iphi  [MAXNCT];
  int   ct_nhcon [MAXNCT];
  int   ct_necon [MAXNCT];
                          
  float ct_eEm   [MAXNCT];
  float ct_eHad  [MAXNCT];
  float ct_eOut  [MAXNCT];

  void init() {

    run     = -999;
    event   = -999;
    nct     = -999;

    for (int iCT = 0; iCT < MAXNCT; iCT++){
 
      ct_ieta  [iCT] = -999;
      ct_iphi  [iCT] = -999;
      ct_nhcon [iCT] = -999;
      ct_necon [iCT] = -999;
      
      ct_eEm   [iCT] = -999.0;
      ct_eHad  [iCT] = -999.0;
      ct_eOut  [iCT] = -999.0;

    }
    
  }
  
  
 private:

  CaloTowerCompTree(const CaloTowerCompTree&); // stop default
  const CaloTowerCompTree& operator=(const CaloTowerCompTree&); // stop default

};


#endif
