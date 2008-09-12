#ifndef Analyzers_DijetAnalyzer_MyTree_h
#define Analyzers_DijetAnalyzer_MyTree_h

class MyTree {
  
 public:
   
  enum {MAXNJETS      = 100  };
  enum {MAXNJTOWERS   = 100  };
  enum {MAXNGTOWERS   = 1000 };
  enum {MAXNRECHITS   = 20   };
  enum {MAXTIMESLICES = 10   };
  enum {MAXNSAMPLES   = 10   };
  enum {MAXNGTOWERCONS= 5    };
  
  MyTree(); 
  virtual ~MyTree();
  
  int run;
  int event;
  
  int nrjet;
  
  float jet_p        [MAXNJETS];
  float jet_pt       [MAXNJETS];
  float jet_px       [MAXNJETS];
  float jet_py       [MAXNJETS];
  float jet_pz       [MAXNJETS];	             
  float jet_e        [MAXNJETS];
  float jet_et       [MAXNJETS];  	             
  float jet_eta      [MAXNJETS];
  float jet_phi      [MAXNJETS];
  int   jet_njtowers [MAXNJETS];
  
  int ngtower;
  
  int   genTower_ieta      [MAXNGTOWERS];
  int   genTower_ietaAbs   [MAXNGTOWERS];
  int   genTower_iphi      [MAXNGTOWERS];            
  int   genTower_zside     [MAXNGTOWERS];
  int   genTower_ncon      [MAXNGTOWERS];
  
  float genTower_phi       [MAXNGTOWERS];
  float genTower_eta       [MAXNGTOWERS];
  float genTower_hadEnergy [MAXNGTOWERS];
  float genTower_emEnergy  [MAXNGTOWERS];
  float genTower_hadEt     [MAXNGTOWERS];
  float genTower_emEt      [MAXNGTOWERS];

  int   genTower_con_subdet[MAXNGTOWERS][MAXNGTOWERCONS];
  int   genTower_con_ieta  [MAXNGTOWERS][MAXNGTOWERCONS];
  int   genTower_con_iphi  [MAXNGTOWERS][MAXNGTOWERCONS];
  int   genTower_con_depth [MAXNGTOWERS][MAXNGTOWERCONS];
  int   genTower_con_zside [MAXNGTOWERS][MAXNGTOWERCONS];

  int   jetTower_ieta      [MAXNJETS][MAXNJTOWERS];
  int   jetTower_ietaAbs   [MAXNJETS][MAXNJTOWERS];
  int   jetTower_iphi      [MAXNJETS][MAXNJTOWERS];            
  int   jetTower_zside     [MAXNJETS][MAXNJTOWERS];
  int   jetTower_nrhits    [MAXNJETS][MAXNJTOWERS];
					 
  float jetTower_phi       [MAXNJETS][MAXNJTOWERS];
  float jetTower_eta       [MAXNJETS][MAXNJTOWERS];
  float jetTower_hadEnergy [MAXNJETS][MAXNJTOWERS];
  float jetTower_hadEt     [MAXNJETS][MAXNJTOWERS];
  float jetTower_emEnergy  [MAXNJETS][MAXNJTOWERS];
  float jetTower_emEt      [MAXNJETS][MAXNJTOWERS];

  int   jetTower_depth     [MAXNJETS][MAXNJTOWERS][MAXNRECHITS];
  int   jetTower_rhSubdet  [MAXNJETS][MAXNJTOWERS][MAXNRECHITS];
  int   jetTower_ntslice   [MAXNJETS][MAXNJTOWERS][MAXNRECHITS];
  float jetTower_rhEnergy  [MAXNJETS][MAXNJTOWERS][MAXNRECHITS];

  int   jetTower_digi_ADC  [MAXNJETS][MAXNJTOWERS][MAXNRECHITS][MAXTIMESLICES];
  float jetTower_digi_fC   [MAXNJETS][MAXNJTOWERS][MAXNRECHITS][MAXTIMESLICES];
  
  void init() {

    run     = -999;
    event   = -999;
    nrjet   = -999;
    
    for (int i = 0; i<MAXNJETS; i++){
      
      jet_p  [i] = -999.0;
      jet_p  [i] = -999.0;
      jet_pt [i] = -999.0;
      jet_px [i] = -999.0;
      jet_py [i] = -999.0;
      jet_pz [i] = -999.0;
      
      jet_e  [i] = -999.0;
      jet_et [i] = -999.0;
      
      jet_eta[i] = -999.0;
      jet_phi[i] = -999.0;
      
      jet_njtowers[i] = -999;
      
      for (int j = 0; j<MAXNJTOWERS; j++){
	
	jetTower_ieta     [i][j] = -999;
	jetTower_ietaAbs  [i][j] = -999;
	jetTower_iphi     [i][j] = -999;
	jetTower_zside    [i][j] = -999;
	
	jetTower_phi      [i][j] = -999.0;
	jetTower_eta      [i][j] = -999.0;
	jetTower_hadEnergy[i][j] = -999.0;
	jetTower_hadEt    [i][j] = -999.0;
	jetTower_emEnergy [i][j] = -999.0;
	jetTower_emEt     [i][j] = -999.0;

	jetTower_nrhits   [i][j] = -999;

	for (int k = 0; k<MAXNRECHITS; k++){
	  
	  jetTower_rhSubdet[i][j][k] = -999;
	  jetTower_rhEnergy[i][j][k] = -999.0;
	  jetTower_depth   [i][j][k] = -999;
	  jetTower_ntslice [i][j][k] = -999;

	  for (int l = 0; l < MAXTIMESLICES; l++){
	    jetTower_digi_fC [i][j][k][l] = -999.0;
	    jetTower_digi_ADC[i][j][k][l] = -999;
	  }
	}
      }
    }

    ngtower = -999;

    for (int i = 0 ; i<MAXNGTOWERS;i++){
      genTower_ieta     [i] = -999;
      genTower_ietaAbs  [i] = -999;
      genTower_iphi     [i] = -999;
      genTower_zside    [i] = -999;
      genTower_ncon     [i] = -999;
                     
      genTower_phi      [i] = -999.0;
      genTower_eta      [i] = -999.0;
      genTower_hadEnergy[i] = -999.0;
      genTower_hadEt    [i] = -999.0;
      genTower_emEnergy [i] = -999.0;
      genTower_emEt     [i] = -999.0;
      
      for (int j = 0; j < MAXNGTOWERCONS; j++){

	genTower_con_subdet[i][j] = -999;
	genTower_con_ieta  [i][j] = -999;
	genTower_con_iphi  [i][j] = -999;
	genTower_con_depth [i][j] = -999;
	genTower_con_zside [i][j] = -999;

      }
    }
  }

  
 private:

  MyTree(const MyTree&); // stop default
  const MyTree& operator=(const MyTree&); // stop default

};


#endif
