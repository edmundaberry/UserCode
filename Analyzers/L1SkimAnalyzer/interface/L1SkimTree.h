#ifndef Analyzers_L1SkimAnalyzer_L1SkimTree_h
#define Analyzers_L1SkimAnalyzer_L1SkimTree_h

class L1SkimTree
{

 public:
  
  enum {MAXNCALOEM = 200};

  L1SkimTree();
  virtual ~L1SkimTree();

  int nCaloEm;

  int caloEm_ieta     [MAXNCALOEM]; 
  int caloEm_iphi     [MAXNCALOEM]; 
  int caloEm_rank     [MAXNCALOEM]; 
  int caloEm_raw      [MAXNCALOEM]; 
  int caloEm_rctCrate [MAXNCALOEM]; 
  int caloEm_rctRegion[MAXNCALOEM]; 
  int caloEm_rctCard  [MAXNCALOEM]; 
  int caloEm_cindex   [MAXNCALOEM]; 
  int caloEm_bx       [MAXNCALOEM]; 
  int caloEm_iso      [MAXNCALOEM]; 

  void init (){

    nCaloEm = -999;

    for (int i = 0; i < MAXNCALOEM; i++){

      caloEm_ieta     [i] = -999;
      caloEm_iphi     [i] = -999;
      caloEm_rank     [i] = -999;
      caloEm_raw      [i] = -999;
      caloEm_rctCrate [i] = -999;
      caloEm_rctRegion[i] = -999;
      caloEm_rctCard  [i] = -999;
      caloEm_cindex   [i] = -999;
      caloEm_bx       [i] = -999;
      caloEm_iso      [i] = -999;

    }
  }

 private:
  
  L1SkimTree (const L1SkimTree&);
  const L1SkimTree& operator=(const L1SkimTree&);

};

#endif
