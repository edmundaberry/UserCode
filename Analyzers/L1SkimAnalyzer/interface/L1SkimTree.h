#ifndef Analyzers_L1SkimAnalyzer_L1SkimTree_h
#define Analyzers_L1SkimAnalyzer_L1SkimTree_h

class L1SkimTree
{

 public:
  
  enum {MAXNCALOEM    = 200};
  enum {MAXNGCTCENJET = 4  };
  enum {MAXNGCTFORJET = 4  };
  enum {MAXNGCTTAUJET = 4  };

  L1SkimTree();
  virtual ~L1SkimTree();

  int   nCaloEm;
        
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

  int   nGctCenJet;

  int   gctCenJet_ietaSign [MAXNGCTCENJET];
  int   gctCenJet_ieta     [MAXNGCTCENJET];
  int   gctCenJet_iphi     [MAXNGCTCENJET];
  int   gctCenJet_capBlock [MAXNGCTCENJET];
  int   gctCenJet_capIndex [MAXNGCTCENJET];
  int   gctCenJet_bx       [MAXNGCTCENJET];
  int   gctCenJet_rank     [MAXNGCTCENJET];
  float gctCenJet_et       [MAXNGCTCENJET];

  int   nGctForJet;
 
  int   gctForJet_ietaSign [MAXNGCTFORJET];
  int   gctForJet_ieta     [MAXNGCTFORJET];
  int   gctForJet_iphi     [MAXNGCTFORJET];
  int   gctForJet_capBlock [MAXNGCTFORJET];
  int   gctForJet_capIndex [MAXNGCTFORJET];
  int   gctForJet_bx       [MAXNGCTFORJET];
  int   gctForJet_rank     [MAXNGCTFORJET];
  float gctForJet_et       [MAXNGCTFORJET];

  int   nGctTauJet;
  
  int   gctTauJet_ietaSign [MAXNGCTTAUJET];
  int   gctTauJet_ieta     [MAXNGCTTAUJET];
  int   gctTauJet_iphi     [MAXNGCTTAUJET];
  int   gctTauJet_capBlock [MAXNGCTTAUJET];
  int   gctTauJet_capIndex [MAXNGCTTAUJET];
  int   gctTauJet_bx       [MAXNGCTTAUJET];
  int   gctTauJet_rank     [MAXNGCTTAUJET];
  float gctTauJet_et       [MAXNGCTTAUJET];

  void init (){

    nCaloEm    = -999;
    nGctCenJet = -999;
    nGctTauJet = -999;

    for (int i = 0; i < MAXNGCTCENJET; i++){

      gctCenJet_ietaSign [i] = -999;
      gctCenJet_ieta     [i] = -999;
      gctCenJet_iphi     [i] = -999;
      gctCenJet_capBlock [i] = -999;
      gctCenJet_capIndex [i] = -999;
      gctCenJet_bx       [i] = -999;
      gctCenJet_rank     [i] = -999;
      gctCenJet_et       [i] = -999;

    }
    
    for (int i = 0; i < MAXNGCTFORJET; i++){

      gctForJet_ietaSign [i] = -999;
      gctForJet_ieta     [i] = -999;
      gctForJet_iphi     [i] = -999;
      gctForJet_capBlock [i] = -999;
      gctForJet_capIndex [i] = -999;
      gctForJet_bx       [i] = -999;
      gctForJet_rank     [i] = -999;
      gctForJet_et       [i] = -999;

    }

    for (int i = 0; i < MAXNGCTTAUJET; i++){

      gctTauJet_ietaSign [i] = -999;
      gctTauJet_ieta     [i] = -999;
      gctTauJet_iphi     [i] = -999;
      gctTauJet_capBlock [i] = -999;
      gctTauJet_capIndex [i] = -999;
      gctTauJet_bx       [i] = -999;
      gctTauJet_rank     [i] = -999;
      gctTauJet_et       [i] = -999;

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
  }

 private:
  
  L1SkimTree (const L1SkimTree&);
  const L1SkimTree& operator=(const L1SkimTree&);

};

#endif
