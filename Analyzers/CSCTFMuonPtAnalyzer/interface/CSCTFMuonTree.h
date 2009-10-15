#ifndef ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONTREE_H
#define ANALYZERS_CSCTFMUONPTANALYZER_CSCTFMUONTREE_H

class CSCTFMuonTree {
  
 public:
  
  CSCTFMuonTree(); 
  virtual ~CSCTFMuonTree();

  enum { MAXNSIMHIT  = 100  };
  enum { MAXNGMU     = 1    }; 
  enum { MAXNSTUB    = 10   };
  enum { MAXNSECTOR  = 10   }; 
  enum { MAXNCOMBO   = 11   };

  int run;
  int event;

  int nsh;
  int ngmu;
  int nstub;
  int nsector;

  int quality;
  int ptBin;	 
  int evtComboHitID;
  int evtComboFrbID;
  
  float gmu_pt                 [MAXNGMU];
  float gmu_eta                [MAXNGMU];
  float gmu_phi                [MAXNGMU];
  
  float stub_digi_eta          [MAXNSTUB];
  float stub_digi_phi          [MAXNSTUB];
  int   stub_digi_strip        [MAXNSTUB];
  int   stub_digi_keyWG        [MAXNSTUB];
  int   stub_digi_quality      [MAXNSTUB];
  int   stub_digi_badphi       [MAXNSTUB];
  int   stub_digi_pattern      [MAXNSTUB];
  int   stub_digi_bend         [MAXNSTUB];
    			       		    
  int   stub_digi_lclPhi       [MAXNSTUB];
  int   stub_digi_lclPhiBend   [MAXNSTUB];
  int   stub_digi_gblPhi       [MAXNSTUB];
  int   stub_digi_gblEta       [MAXNSTUB];
			       
  int   stub_frBit             [MAXNSTUB];
  int   stub_cscid             [MAXNSTUB];
  int   stub_stat              [MAXNSTUB];
  int   stub_ring              [MAXNSTUB];
  int   stub_cham              [MAXNSTUB];
  int   stub_endc              [MAXNSTUB];
  int   stub_sect              [MAXNSTUB];

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

  int   sect_quality           [MAXNSECTOR];
  int   sect_ncombo            [MAXNSECTOR];
  int   sect_num               [MAXNSECTOR];
  int   sect_combo_hitId       [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_frbId       [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_isFirst     [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_etaBin      [MAXNSECTOR][MAXNCOMBO]; 
  int   sect_combo_s1ring      [MAXNSECTOR][MAXNCOMBO]; 
  int   sect_combo_dphi        [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_deta        [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_phi1        [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_eta1        [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_stat1       [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_phi2        [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_eta2        [MAXNSECTOR][MAXNCOMBO];
  int   sect_combo_stat2       [MAXNSECTOR][MAXNCOMBO];

  void init() {

    run   = -999;
    event = -999;
    
    nstub    = 0;
    ngmu     = 0;
    nsh      = 0;
    nsector  = 0;

    quality       = -999;
    ptBin         = -999;
    evtComboHitID = -1;
    evtComboFrbID = -1;
    
    for (int i = 0; i < MAXNGMU; i++){
      
      gmu_pt [i] = -999.;
      gmu_eta[i] = -999.;
      gmu_phi[i] = -999.;
      
    }

    for (int j = 0; j < MAXNSTUB; j++){
      
      stub_frBit [j] = -999;
      stub_cscid [j] = -999;
      stub_stat  [j] = -999;
      stub_ring  [j] = -999;
      stub_cham  [j] = -999;
      stub_endc  [j] = -999;      
      stub_sect  [j] = -999;      
      
      stub_digi_eta       [j] = -999.;
      stub_digi_phi       [j] = -999.;
      stub_digi_strip     [j] = -999;
      stub_digi_keyWG     [j] = -999;
      stub_digi_quality   [j] = -999;
      stub_digi_badphi    [j] = -999;
      stub_digi_pattern   [j] = -999;
      stub_digi_bend      [j] = -999;
      stub_digi_lclPhi    [j] = -999;
      stub_digi_lclPhiBend[j] = -999;
      stub_digi_gblPhi    [j] = -999;
      stub_digi_gblEta    [j] = -999;
      
    }

    for (int l = 0; l < MAXNSECTOR; l++){

      sect_ncombo [l] = 0;
      sect_num    [l] = -999;
      sect_quality[l] = -999;

      for (int l1 = 0; l1 < MAXNCOMBO; ++l1 ){
	
	sect_combo_isFirst [l][l1] = -999;
	sect_combo_hitId   [l][l1] = -999;
	sect_combo_frbId   [l][l1] = -999;
	sect_combo_dphi    [l][l1] = -999;
	sect_combo_deta    [l][l1] = -999;
	sect_combo_etaBin  [l][l1] = -999;	
	sect_combo_s1ring  [l][l1] = -999;	
	sect_combo_eta1    [l][l1] = -999;
	sect_combo_phi1    [l][l1] = -999;
	sect_combo_stat1   [l][l1] = -999;	
	sect_combo_eta2    [l][l1] = -999;
	sect_combo_phi2    [l][l1] = -999;
	sect_combo_stat2   [l][l1] = -999;	

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

