#ifndef Analyzers_CaloTowerAna_CaloTowerTree_h
#define Analyzers_CaloTowerAna_CaloTowerTree_h

class CaloTowerTree {
  
 public:
   
  enum {MAXNJETS      = 100  };
  enum {MAXNJTOWERS   = 100  };
  enum {MAXNGTOWERS   = 1000 };
  enum {MAXNDIGIS     = 10   };
  enum {MAXNRECHITS   = 20   };
  enum {MAXTIMESLICES = 10   };
  enum {MAXNSAMPLES   = 10   };
  enum {MAXNGTOWERCONS= 5    };
  
  CaloTowerTree(); 
  virtual ~CaloTowerTree();
  
  int run;
  int event;  
  int nrjet;
  		            				
  int   rjet_nct            [MAXNJETS];
  int   rjet_ct_ieta        [MAXNJETS][MAXNJTOWERS];
  int   rjet_ct_iphi        [MAXNJETS][MAXNJTOWERS];
  int   rjet_ct_ndigi       [MAXNJETS][MAXNJTOWERS];
			    
  int   rjet_ct_digi_ieta   [MAXNJETS][MAXNJTOWERS][MAXNDIGIS];
  int   rjet_ct_digi_iphi   [MAXNJETS][MAXNJTOWERS][MAXNDIGIS];
  int   rjet_ct_digi_depth  [MAXNJETS][MAXNJTOWERS][MAXNDIGIS];
  int   rjet_ct_digi_size   [MAXNJETS][MAXNJTOWERS][MAXNDIGIS];
  int   rjet_ct_digi_ps     [MAXNJETS][MAXNJTOWERS][MAXNDIGIS];

  int   rjet_ct_digi_adc    [MAXNJETS][MAXNJTOWERS][MAXNDIGIS][MAXNSAMPLES];
  int   rjet_ct_digi_capid  [MAXNJETS][MAXNJTOWERS][MAXNDIGIS][MAXNSAMPLES];

  float rjet_p              [MAXNJETS];
  float rjet_pt             [MAXNJETS];
  float rjet_e              [MAXNJETS];
  float rjet_et             [MAXNJETS];  	             
  float rjet_eta            [MAXNJETS];
  float rjet_phi            [MAXNJETS];

  float rjet_ct_emE         [MAXNJETS][MAXNJTOWERS];
  float rjet_ct_hadE        [MAXNJETS][MAXNJTOWERS];
  float rjet_ct_outE        [MAXNJETS][MAXNJTOWERS];

  float rjet_ct_digi_fC     [MAXNJETS][MAXNJTOWERS][MAXNDIGIS][MAXNSAMPLES];
  float rjet_ct_digi_ped    [MAXNJETS][MAXNJTOWERS][MAXNDIGIS][MAXNSAMPLES];
  float rjet_ct_digi_gain   [MAXNJETS][MAXNJTOWERS][MAXNDIGIS][MAXNSAMPLES];
  float rjet_ct_digi_rcgain [MAXNJETS][MAXNJTOWERS][MAXNDIGIS][MAXNSAMPLES];

  void init() {

    run     = -999;
    event   = -999;
    nrjet   = -999;
    
    for (int i = 0; i<MAXNJETS; i++){
      
      rjet_p  [i] = -999.0;
      rjet_pt [i] = -999.0;
      rjet_e  [i] = -999.0;
      rjet_et [i] = -999.0;      
      rjet_eta[i] = -999.0;
      rjet_phi[i] = -999.0;

      rjet_nct[i] = -999;

      for (int j = 0; j < MAXNJTOWERS; j++){
	
	rjet_ct_ieta [i][j] = -999;
	rjet_ct_iphi [i][j] = -999;
	rjet_ct_ndigi[i][j] = -999;

	rjet_ct_emE  [i][j] = -999.0;
	rjet_ct_hadE [i][j] = -999.0;
	rjet_ct_outE [i][j] = -999.0;
	
	for (int k = 0; k < MAXNDIGIS; k++){

	  rjet_ct_digi_ieta  [i][j][k] = -999;
	  rjet_ct_digi_iphi  [i][j][k] = -999;
	  rjet_ct_digi_depth [i][j][k] = -999;
	  rjet_ct_digi_size  [i][j][k] = -999;
	  rjet_ct_digi_ps    [i][j][k] = -999;

	  for (int l = 0; l < MAXNSAMPLES; l++){

	    rjet_ct_digi_adc    [i][j][k][l] = -999;
	    rjet_ct_digi_capid  [i][j][k][l] = -999;
	    rjet_ct_digi_fC     [i][j][k][l] = -999.0;
	    rjet_ct_digi_ped    [i][j][k][l] = -999.0;
	    rjet_ct_digi_gain   [i][j][k][l] = -999.0;
	    rjet_ct_digi_rcgain [i][j][k][l] = -999.0;

	  }	  
	}	
      }      
    }
  }
  
  
 private:

  CaloTowerTree(const CaloTowerTree&); // stop default
  const CaloTowerTree& operator=(const CaloTowerTree&); // stop default

};


#endif
