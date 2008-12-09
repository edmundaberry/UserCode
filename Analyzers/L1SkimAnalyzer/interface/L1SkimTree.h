

#ifndef Analyzers_L1SkimAnalyzer_L1SkimTree_h
#define Analyzers_L1SkimAnalyzer_L1SkimTree_h

class L1SkimTree
{

 public:
  
  enum {MAXNCALOEM    = 200};
  enum {MAXNHLTJETS   = 100};
  enum {MAXNHLTJTOWERS= 100};
  enum {MAXNGCTHT     = 12 };
  enum {MAXNGCTCENJET = 12 };
  enum {MAXNGCTFORJET = 12 };
  enum {MAXNGCTTAUJET = 12 };
  enum {MAXNL1CENJET  = 12 };
  enum {MAXNL1FORJET  = 12 };
  enum {MAXNL1TAUJET  = 12 };

  L1SkimTree();
  virtual ~L1SkimTree();

  int   run;
  int   event;

  int   nCaloEm;
  int   nHLTJetCands;					 
  int   nGctCenJet;
  int   nGctForJet;
  int   nGctTauJet;
  int   nL1GctEtHads;
  int   nL1CenJet;
  int   nL1TauJet;
  int   nL1ForJet;
	  
  int   l1_HTT100;
  int   l1_HTT200;
  int   l1_HTT300;
  int   l1_HTT400;
  int   l1_HTT500;

  int   l1_SingleJet15;
  int   l1_SingleJet20;
  int   l1_SingleJet30;
  int   l1_SingleJet50;

  int   caloEm_ieta        [MAXNCALOEM]; 
  int   caloEm_iphi        [MAXNCALOEM]; 
  int   caloEm_rank        [MAXNCALOEM]; 
  int   caloEm_raw         [MAXNCALOEM]; 
  int   caloEm_rctCrate    [MAXNCALOEM]; 
  int   caloEm_rctRegion   [MAXNCALOEM]; 
  int   caloEm_rctCard     [MAXNCALOEM]; 
  int   caloEm_cindex      [MAXNCALOEM]; 
  int   caloEm_bx          [MAXNCALOEM]; 
  int   caloEm_iso         [MAXNCALOEM]; 

  int   gctCenJet_iphiMin  [MAXNGCTCENJET];
  int   gctCenJet_iphiMax  [MAXNGCTCENJET];
  int   gctCenJet_ietaMin  [MAXNGCTCENJET];
  int   gctCenJet_ietaMax  [MAXNGCTCENJET];
  int   gctCenJet_ietaSign [MAXNGCTCENJET];
  int   gctCenJet_etaIndex [MAXNGCTCENJET];
  int   gctCenJet_phiIndex [MAXNGCTCENJET];
  int   gctCenJet_capBlock [MAXNGCTCENJET];
  int   gctCenJet_capIndex [MAXNGCTCENJET];
  int   gctCenJet_bx       [MAXNGCTCENJET];
  int   gctCenJet_rank     [MAXNGCTCENJET];
 
  int   gctForJet_iphiMin  [MAXNGCTFORJET];
  int   gctForJet_iphiMax  [MAXNGCTFORJET];
  int   gctForJet_ietaMin  [MAXNGCTCENJET];
  int   gctForJet_ietaMax  [MAXNGCTCENJET];
  int   gctForJet_ietaSign [MAXNGCTFORJET];
  int   gctForJet_etaIndex [MAXNGCTFORJET];
  int   gctForJet_phiIndex [MAXNGCTFORJET];
  int   gctForJet_capBlock [MAXNGCTFORJET];
  int   gctForJet_capIndex [MAXNGCTFORJET];
  int   gctForJet_bx       [MAXNGCTFORJET];
  int   gctForJet_rank     [MAXNGCTFORJET];
  
  int   gctTauJet_iphiMin  [MAXNGCTTAUJET];
  int   gctTauJet_iphiMax  [MAXNGCTTAUJET];
  int   gctTauJet_ietaMin  [MAXNGCTCENJET];
  int   gctTauJet_ietaMax  [MAXNGCTCENJET];
  int   gctTauJet_ietaSign [MAXNGCTTAUJET];
  int   gctTauJet_etaIndex [MAXNGCTTAUJET];
  int   gctTauJet_phiIndex [MAXNGCTTAUJET];
  int   gctTauJet_capBlock [MAXNGCTTAUJET];
  int   gctTauJet_capIndex [MAXNGCTTAUJET];
  int   gctTauJet_bx       [MAXNGCTTAUJET];
  int   gctTauJet_rank     [MAXNGCTTAUJET];

  int   nHLTJetTowers      [MAXNHLTJETS];		 
  int   hltJetTower_ieta   [MAXNHLTJETS][MAXNHLTJTOWERS];
  int   hltJetTower_iphi   [MAXNHLTJETS][MAXNHLTJTOWERS];

  float l1CenJet_px        [MAXNL1CENJET];
  float l1CenJet_py        [MAXNL1CENJET];
  float l1CenJet_pz        [MAXNL1CENJET];
  float l1CenJet_pt        [MAXNL1CENJET];
  float l1CenJet_et        [MAXNL1CENJET];
  float l1CenJet_eta       [MAXNL1CENJET];
  float l1CenJet_phi       [MAXNL1CENJET];

  float l1TauJet_px        [MAXNL1TAUJET];
  float l1TauJet_py        [MAXNL1TAUJET];
  float l1TauJet_pz        [MAXNL1TAUJET];
  float l1TauJet_pt        [MAXNL1TAUJET];
  float l1TauJet_et        [MAXNL1TAUJET];
  float l1TauJet_eta       [MAXNL1TAUJET];
  float l1TauJet_phi       [MAXNL1TAUJET];

  float l1ForJet_px        [MAXNL1FORJET];
  float l1ForJet_py        [MAXNL1FORJET];
  float l1ForJet_pz        [MAXNL1FORJET];
  float l1ForJet_pt        [MAXNL1FORJET];
  float l1ForJet_et        [MAXNL1FORJET];
  float l1ForJet_eta       [MAXNL1FORJET];
  float l1ForJet_phi       [MAXNL1FORJET];

  float gctCenJet_et       [MAXNGCTCENJET];
  float gctForJet_et       [MAXNGCTFORJET];
  float gctTauJet_et       [MAXNGCTTAUJET];
  
  float gctHT              [MAXNGCTHT];
  float gctHT_UnCorr       [MAXNGCTHT];
  
  float hltJet_pt          [MAXNHLTJETS];		 
  float hltJet_et          [MAXNHLTJETS];		 
  float hltJet_phi         [MAXNHLTJETS];		 
  float hltJet_eta         [MAXNHLTJETS];		 	                                               

  void init (){
    
    run            = -999;
    event          = -999;
		   
    nL1GctEtHads   = -999;
    nCaloEm        = -999;
    nGctCenJet     = -999;
    nGctTauJet     = -999;
    nHLTJetCands   = -999;
    nL1CenJet      = -999;
    nL1TauJet      = -999;
    nL1ForJet      = -999;

    l1_SingleJet15 = -999;
    l1_SingleJet20 = -999;
    l1_SingleJet30 = -999;
    l1_SingleJet50 = -999;

    l1_HTT100      = -999;
    l1_HTT200      = -999;
    l1_HTT300      = -999;
    l1_HTT400      = -999;
    l1_HTT500      = -999;

    for (int i = 0; i < MAXNL1CENJET; i++){

      l1CenJet_px        [i] = -999.0;
      l1CenJet_py        [i] = -999.0;
      l1CenJet_pz        [i] = -999.0;
      l1CenJet_pt        [i] = -999.0;
      l1CenJet_et        [i] = -999.0;
      l1CenJet_eta       [i] = -999.0;
      l1CenJet_phi       [i] = -999.0;
      
    }
 
    for (int i = 0; i < MAXNL1TAUJET; i++){

      l1TauJet_px        [i] = -999.0;
      l1TauJet_py        [i] = -999.0;
      l1TauJet_pz        [i] = -999.0;
      l1TauJet_pt        [i] = -999.0;
      l1TauJet_et        [i] = -999.0;
      l1TauJet_eta       [i] = -999.0;
      l1TauJet_phi       [i] = -999.0;
      
    }

    for (int i = 0; i < MAXNL1FORJET; i++){

      l1ForJet_px        [i] = -999.0;
      l1ForJet_py        [i] = -999.0;
      l1ForJet_pz        [i] = -999.0;
      l1ForJet_pt        [i] = -999.0;
      l1ForJet_et        [i] = -999.0;
      l1ForJet_eta       [i] = -999.0;
      l1ForJet_phi       [i] = -999.0;
      
    }
    

    for (int i = 0; i < MAXNGCTHT; i++){

      gctHT              [i] = -999;
      gctHT_UnCorr       [i] = -999;

    }

    for (int i = 0; i < MAXNGCTCENJET; i++){

      gctCenJet_iphiMin  [i] = -999;
      gctCenJet_iphiMax  [i] = -999;
      gctCenJet_ietaMin  [i] = -999;
      gctCenJet_ietaMax  [i] = -999;
      gctCenJet_ietaSign [i] = -999;
      gctCenJet_etaIndex [i] = -999;
      gctCenJet_phiIndex [i] = -999;
      gctCenJet_capBlock [i] = -999;
      gctCenJet_capIndex [i] = -999;
      gctCenJet_bx       [i] = -999;
      gctCenJet_rank     [i] = -999;
      gctCenJet_et       [i] = -999.0;

    }
    
    for (int i = 0; i < MAXNGCTFORJET; i++){

      gctForJet_iphiMin  [i] = -999;
      gctForJet_iphiMax  [i] = -999;
      gctForJet_ietaMin  [i] = -999;
      gctForJet_ietaMax  [i] = -999;
      gctForJet_ietaSign [i] = -999;
      gctForJet_etaIndex [i] = -999;
      gctForJet_phiIndex [i] = -999;
      gctForJet_capBlock [i] = -999;
      gctForJet_capIndex [i] = -999;
      gctForJet_bx       [i] = -999;
      gctForJet_rank     [i] = -999;
      gctForJet_et       [i] = -999.0;

    }

    for (int i = 0; i < MAXNGCTTAUJET; i++){

      gctTauJet_iphiMin  [i] = -999;
      gctTauJet_iphiMax  [i] = -999;
      gctTauJet_ietaMin  [i] = -999;
      gctTauJet_ietaMax  [i] = -999;
      gctTauJet_ietaSign [i] = -999;
      gctTauJet_etaIndex [i] = -999;
      gctTauJet_phiIndex [i] = -999;
      gctTauJet_capBlock [i] = -999;
      gctTauJet_capIndex [i] = -999;
      gctTauJet_bx       [i] = -999;
      gctTauJet_rank     [i] = -999;
      gctTauJet_et       [i] = -999.0;

    }
    
    for (int i = 0; i < MAXNCALOEM; i++){

      caloEm_ieta        [i] = -999;
      caloEm_iphi        [i] = -999;
      caloEm_rank        [i] = -999;
      caloEm_raw         [i] = -999;
      caloEm_rctCrate    [i] = -999;
      caloEm_rctRegion   [i] = -999;
      caloEm_rctCard     [i] = -999;
      caloEm_cindex      [i] = -999;
      caloEm_bx          [i] = -999;
      caloEm_iso         [i] = -999;

    }

    for (int i = 0; i < MAXNHLTJETS; i++){

      hltJet_pt      [i] = -999.0;
      hltJet_et      [i] = -999.0;
      hltJet_phi     [i] = -999.0;
      hltJet_eta     [i] = -999.0;
      nHLTJetTowers  [i] = -999;

      for (int j = 0; j < MAXNHLTJTOWERS; j++){
	
	hltJetTower_ieta [i][j] = -999;  
	hltJetTower_iphi [i][j] = -999;
 
      }
    }
  }

 private:
  
  L1SkimTree (const L1SkimTree&);
  const L1SkimTree& operator=(const L1SkimTree&);

};

#endif
