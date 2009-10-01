#ifndef ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONTREE_H
#define ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONTREE_H

class CSCTFMuonTree {
  
 public:
  
  CSCTFMuonTree(); 
  virtual ~CSCTFMuonTree();

  enum { MAXNL1DETID = 10   }; 
  enum { MAXNGMU     = 1    }; 
  enum { MAXNL1DIGI  = 2    };
  enum { MAXNCOMBO   = 11   };
  enum { MAXNSIMHIT  = 100  };

  int run;
  int event;

  int nsh;
  int ngmu;
  int nl1detid;
  int ncombo;
  int filled;

  int nskip;
  int ptBin;	 
  int evtComboHitID;
  int evtComboFrbID;
  
  float gmu_pt                 [MAXNGMU];
  float gmu_eta                [MAXNGMU];
  float gmu_phi                [MAXNGMU];
  
  float l1detid_digi_eta       [MAXNL1DETID];
  float l1detid_digi_phi       [MAXNL1DETID];
  int   l1detid_digi_strip     [MAXNL1DETID];
  int   l1detid_digi_keyWG     [MAXNL1DETID];
  int   l1detid_digi_quality   [MAXNL1DETID];
  int   l1detid_digi_badphi    [MAXNL1DETID];
  int   l1detid_digi_pattern   [MAXNL1DETID];
  int   l1detid_digi_bend      [MAXNL1DETID];
    					    
  int   l1detid_digi_lclPhi    [MAXNL1DETID];
  int   l1detid_digi_lclPhiBend[MAXNL1DETID];
  int   l1detid_digi_gblPhi    [MAXNL1DETID];
  int   l1detid_digi_gblEta    [MAXNL1DETID];

  int   l1detid_frBit          [MAXNL1DETID];
  int   l1detid_cscid          [MAXNL1DETID];
  int   l1detid_stat           [MAXNL1DETID];
  int   l1detid_ring           [MAXNL1DETID];
  int   l1detid_cham           [MAXNL1DETID];
  int   l1detid_layr           [MAXNL1DETID];
  int   l1detid_endc           [MAXNL1DETID];
  int   l1detid_sect           [MAXNL1DETID];

  float sh_eta                 [MAXNSIMHIT];
  float sh_phi                 [MAXNSIMHIT];
  int   sh_pdg                 [MAXNSIMHIT];
  int   sh_frBit               [MAXNSIMHIT];
  int   sh_stat                [MAXNSIMHIT];
  int   sh_ring                [MAXNSIMHIT];
  int   sh_cham                [MAXNSIMHIT];
  int   sh_layr                [MAXNSIMHIT];
  int   sh_endc                [MAXNSIMHIT];
  int   sh_sect                [MAXNSIMHIT];
  int   sh_cscid               [MAXNSIMHIT];

  int   combo_hitId            [MAXNCOMBO];
  int   combo_frbId            [MAXNCOMBO];
  int   combo_isFirst          [MAXNCOMBO];
  int   combo_sector           [MAXNCOMBO];
  int   combo_etaBin           [MAXNCOMBO];
  int   combo_dphi             [MAXNCOMBO];
  float combo_phi              [MAXNCOMBO][2];
  float combo_eta              [MAXNCOMBO][2];

  void init() {

    run   = -999;
    event = -999;
    
    ngmu  = 0;
    nl1detid = 0;
    ncombo = 0;
    nsh = 0;
    nskip = 0;

    filled = 0;

    ptBin      = -999;
    evtComboHitID = -1;
    evtComboFrbID = -1;
    
    for (int i = 0; i < MAXNGMU; i++){
      
      gmu_pt [i] = -999.;
      gmu_eta[i] = -999.;
      gmu_phi[i] = -999.;
      
    }

    for (int j = 0; j < MAXNL1DETID; j++){
      
      l1detid_frBit [j] = -999;
      l1detid_cscid [j] = -999;
      l1detid_stat  [j] = -999;
      l1detid_ring  [j] = -999;
      l1detid_cham  [j] = -999;
      l1detid_layr  [j] = -999;
      l1detid_endc  [j] = -999;      
      l1detid_sect  [j] = -999;      
      
      l1detid_digi_eta       [j] = -999.;
      l1detid_digi_phi       [j] = -999.;
      l1detid_digi_strip     [j] = -999;
      l1detid_digi_keyWG     [j] = -999;
      l1detid_digi_quality   [j] = -999;
      l1detid_digi_badphi    [j] = -999;
      l1detid_digi_pattern   [j] = -999;
      l1detid_digi_bend      [j] = -999;
      l1detid_digi_lclPhi    [j] = -999;
      l1detid_digi_lclPhiBend[j] = -999;
      l1detid_digi_gblPhi    [j] = -999;
      l1detid_digi_gblEta    [j] = -999;
      
    }

    for (int l = 0; l < MAXNCOMBO; l++){
      combo_isFirst [l] = -999;
      combo_hitId   [l] = -999;
      combo_frbId   [l] = -999;
      combo_dphi    [l] = -999;
      combo_etaBin  [l] = -999;
      combo_sector  [l] = -999;
      
      for (int m = 0; m < 2; m++){
	combo_eta[l][m] = -999;
	combo_phi[l][m] = -999;
      }
    }

    for (int n = 0; n < MAXNSIMHIT; ++n){
      
      sh_eta   [n] = -999.;
      sh_phi   [n] = -999.;
      sh_pdg   [n] = -999;
      sh_frBit [n] = -999;
      sh_stat  [n] = -999;
      sh_ring  [n] = -999;
      sh_cham  [n] = -999;
      sh_layr  [n] = -999;
      sh_endc  [n] = -999;
      sh_sect  [n] = -999;
      sh_cscid [n] = -999;

    }
    
  }


 private:
  
  CSCTFMuonTree(const CSCTFMuonTree&); // stop default
  const CSCTFMuonTree& operator=(const CSCTFMuonTree&); // stop default

};

#endif

