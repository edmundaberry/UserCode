#ifndef Analyzers_L1SkimAnalyzer_L1SkimTree_h
#define Analyzers_L1SkimAnalyzer_L1SkimTree_h

class L1SkimTree
{
      
   public:
      
      enum {MAXNCALOEM    = 200};
      enum {MAXNGENJETS   = 100};
      enum {MAXNHLTJETS   = 100};
      enum {MAXNHLTJTOWERS= 100};
      enum {MAXNGCTHT     = 12 };
      enum {MAXNL1CENJET  = 12 };
      enum {MAXNL1FORJET  = 12 }; 
      enum {MAXNL1TAUJET  = 12 };
      
      L1SkimTree();
      virtual ~L1SkimTree();
      
      int   run;
      int   event;
      
      int   nGenJets;
      int   nHLTJetCands;					 
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
      
      int   l1CenJet_ietaSign  [MAXNL1CENJET];
      int   l1CenJet_etaIndex  [MAXNL1CENJET];
      int   l1CenJet_phiIndex  [MAXNL1CENJET];
      int   l1CenJet_capBlock  [MAXNL1CENJET];
      int   l1CenJet_capIndex  [MAXNL1CENJET];
      int   l1CenJet_bx        [MAXNL1CENJET];
      int   l1CenJet_rank      [MAXNL1CENJET];        
      
      int   l1TauJet_ietaSign  [MAXNL1TAUJET];
      int   l1TauJet_etaIndex  [MAXNL1TAUJET];
      int   l1TauJet_phiIndex  [MAXNL1TAUJET];
      int   l1TauJet_capBlock  [MAXNL1TAUJET];
      int   l1TauJet_capIndex  [MAXNL1TAUJET];
      int   l1TauJet_bx        [MAXNL1TAUJET];
      int   l1TauJet_rank      [MAXNL1TAUJET];
      
      int   l1ForJet_ietaSign  [MAXNL1FORJET];
      int   l1ForJet_etaIndex  [MAXNL1FORJET];
      int   l1ForJet_phiIndex  [MAXNL1FORJET];
      int   l1ForJet_capBlock  [MAXNL1FORJET];
      int   l1ForJet_capIndex  [MAXNL1FORJET];
      int   l1ForJet_bx        [MAXNL1FORJET];
      int   l1ForJet_rank      [MAXNL1FORJET];
      
      int nHLTJetTowers        [MAXNHLTJETS];		 
      int hltJetTower_ieta     [MAXNHLTJETS][MAXNHLTJTOWERS];
      int hltJetTower_iphi     [MAXNHLTJETS][MAXNHLTJTOWERS];

      float genJet_p           [MAXNGENJETS];
      float genJet_px 	       [MAXNGENJETS];
      float genJet_py 	       [MAXNGENJETS];
      float genJet_pz 	       [MAXNGENJETS];
      float genJet_pt 	       [MAXNGENJETS];
      float genJet_et 	       [MAXNGENJETS];
      float genJet_eta	       [MAXNGENJETS];
      float genJet_phi         [MAXNGENJETS];

      float l1CenJet_etaMin    [MAXNL1CENJET];
      float l1CenJet_etaMax    [MAXNL1CENJET];
      float l1CenJet_phiMin    [MAXNL1CENJET];
      float l1CenJet_phiMax    [MAXNL1CENJET];
      float l1CenJet_px        [MAXNL1CENJET];
      float l1CenJet_py        [MAXNL1CENJET];
      float l1CenJet_pz        [MAXNL1CENJET];
      float l1CenJet_pt        [MAXNL1CENJET];
      float l1CenJet_et        [MAXNL1CENJET];
      float l1CenJet_eta       [MAXNL1CENJET];
      float l1CenJet_phi       [MAXNL1CENJET];
      
      float l1TauJet_etaMin    [MAXNL1TAUJET];
      float l1TauJet_etaMax    [MAXNL1TAUJET];
      float l1TauJet_phiMin    [MAXNL1TAUJET];
      float l1TauJet_phiMax    [MAXNL1TAUJET];
      float l1TauJet_px        [MAXNL1TAUJET];
      float l1TauJet_py        [MAXNL1TAUJET];
      float l1TauJet_pz        [MAXNL1TAUJET];
      float l1TauJet_pt        [MAXNL1TAUJET];
      float l1TauJet_et        [MAXNL1TAUJET];
      float l1TauJet_eta       [MAXNL1TAUJET];
      float l1TauJet_phi       [MAXNL1TAUJET];

      float l1ForJet_etaMin    [MAXNL1FORJET];
      float l1ForJet_etaMax    [MAXNL1FORJET];
      float l1ForJet_phiMin    [MAXNL1FORJET];
      float l1ForJet_phiMax    [MAXNL1FORJET];
      float l1ForJet_px        [MAXNL1FORJET];
      float l1ForJet_py        [MAXNL1FORJET];
      float l1ForJet_pz        [MAXNL1FORJET];
      float l1ForJet_pt        [MAXNL1FORJET];
      float l1ForJet_et        [MAXNL1FORJET];
      float l1ForJet_eta       [MAXNL1FORJET];
      float l1ForJet_phi       [MAXNL1FORJET];
      
      float gctHT              [MAXNGCTHT];
      float gctHT_UnCorr       [MAXNGCTHT];

      float hltHT;
      
      float hltJet_pt          [MAXNHLTJETS];		 
      float hltJet_et          [MAXNHLTJETS];		 
      float hltJet_phi         [MAXNHLTJETS];		 
      float hltJet_eta         [MAXNHLTJETS];		 	                                               
      
      void init (){
	 
	 run            = -999;
	 event          = -999;
	 
	 nGenJets       = -999;
	 nL1GctEtHads   = -999;
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

	 hltHT          = -999.0;

	 for (int i = 0; i < MAXNGENJETS; i++){

	   genJet_p            [i] = -999.0;
	   genJet_px           [i] = -999.0;
	   genJet_py           [i] = -999.0;
	   genJet_pz           [i] = -999.0;
	   genJet_pt           [i] = -999.0;
	   genJet_et           [i] = -999.0;
	   genJet_eta          [i] = -999.0;
	   genJet_phi          [i] = -999.0;

	 }
	 
	 for (int i = 0; i < MAXNL1CENJET; i++){
	    
	    l1CenJet_ietaSign  [i] = -999;
	    l1CenJet_etaIndex  [i] = -999;
	    l1CenJet_phiIndex  [i] = -999;
	    l1CenJet_capBlock  [i] = -999;
	    l1CenJet_capIndex  [i] = -999;
	    l1CenJet_bx        [i] = -999;
	    l1CenJet_rank      [i] = -999;
	  
	    l1CenJet_etaMin    [i] = -999.0;    
	    l1CenJet_etaMax    [i] = -999.0;
	    l1CenJet_phiMin    [i] = -999.0;
	    l1CenJet_phiMax    [i] = -999.0;
	    l1CenJet_px        [i] = -999.0;
	    l1CenJet_py        [i] = -999.0;
	    l1CenJet_pz        [i] = -999.0;
	    l1CenJet_pt        [i] = -999.0;
	    l1CenJet_et        [i] = -999.0;
	    l1CenJet_eta       [i] = -999.0;
	    l1CenJet_phi       [i] = -999.0;
	    
	 }
 
	 for (int i = 0; i < MAXNL1TAUJET; i++){
	    
	    l1TauJet_ietaSign  [i] = -999;
	    l1TauJet_etaIndex  [i] = -999;
	    l1TauJet_phiIndex  [i] = -999;
	    l1TauJet_capBlock  [i] = -999;
	    l1TauJet_capIndex  [i] = -999;
	    l1TauJet_bx        [i] = -999;
	    l1TauJet_rank      [i] = -999;

	    l1TauJet_etaMin    [i] = -999.0;    
	    l1TauJet_etaMax    [i] = -999.0;
	    l1TauJet_phiMin    [i] = -999.0;
	    l1TauJet_phiMax    [i] = -999.0;	    
	    l1TauJet_px        [i] = -999.0;
	    l1TauJet_py        [i] = -999.0;
	    l1TauJet_pz        [i] = -999.0;
	    l1TauJet_pt        [i] = -999.0;
	    l1TauJet_et        [i] = -999.0;
	    l1TauJet_eta       [i] = -999.0;
	    l1TauJet_phi       [i] = -999.0;
      
	 }
	 
	 for (int i = 0; i < MAXNL1FORJET; i++){
	    
	    l1ForJet_ietaSign  [i] = -999;
	    l1ForJet_etaIndex  [i] = -999;
	    l1ForJet_phiIndex  [i] = -999;
	    l1ForJet_capBlock  [i] = -999;
	    l1ForJet_capIndex  [i] = -999;
	    l1ForJet_bx        [i] = -999;
	    l1ForJet_rank      [i] = -999;
	    
	    l1ForJet_etaMin    [i] = -999.0;    
	    l1ForJet_etaMax    [i] = -999.0;
	    l1ForJet_phiMin    [i] = -999.0;
	    l1ForJet_phiMax    [i] = -999.0;
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
	 
	 for (int i = 0; i < MAXNHLTJETS; i++){
	    
	    hltJet_pt          [i] = -999.0;
	    hltJet_et          [i] = -999.0;
	    hltJet_phi         [i] = -999.0;
	    hltJet_eta         [i] = -999.0;
	    nHLTJetTowers      [i] = -999;
	    
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
