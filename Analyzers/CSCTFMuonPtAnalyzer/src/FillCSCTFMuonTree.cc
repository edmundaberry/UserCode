#include <iostream>

#include "Analyzers/CSCTFMuonPtAnalyzer/interface/FillCSCTFMuonTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

FillCSCTFMuonTree::FillCSCTFMuonTree(){}

FillCSCTFMuonTree::~FillCSCTFMuonTree(){}

void FillCSCTFMuonTree::init(std::string filename, CSCTFMuonTree* treePtr) {

  m_file = new TFile( filename.c_str(), "RECREATE" );
  
  m_tree = new TTree( "csctfMuonTree", "CSCTF Muon Info", 1 );
  
  m_treePtr = treePtr;

  m_tree->Branch("run"                  , &(m_treePtr->run                  ), "run/I"                            );
  m_tree->Branch("event"                , &(m_treePtr->event                ), "event/I"                          );

  m_tree->Branch("ngmu"                 , &(m_treePtr->ngmu                 ), "ngmu/I"                           );
  m_tree->Branch("nstub"                , &(m_treePtr->nstub                ), "nstub/I"                          );
  m_tree->Branch("nsector"              , &(m_treePtr->nsector              ), "nsector/I"                        );
  m_tree->Branch("nsh"                  , &(m_treePtr->nsh                  ), "nsh/I"                            );

  m_tree->Branch("ptBin"                , &(m_treePtr->ptBin                ), "ptBin/I"                          );
  m_tree->Branch("evtComboHitID"        , &(m_treePtr->evtComboHitID        ), "evtComboHitID/I"                  );
  m_tree->Branch("evtComboFrbID"        , &(m_treePtr->evtComboFrbID        ), "evtComboFrbID/I"                  );
  
  m_tree->Branch("gmu_pt"               , &(m_treePtr->gmu_pt               ), "gmu_pt[ngmu]/F"                   );
  m_tree->Branch("gmu_eta"              , &(m_treePtr->gmu_eta              ), "gmu_eta[ngmu]/F"                  );
  m_tree->Branch("gmu_phi"              , &(m_treePtr->gmu_phi              ), "gmu_phi[ngmu]/F"                  );
  
  m_tree->Branch("sect_ncombo"          , &(m_treePtr->sect_ncombo          ), "sect_ncombo[nsector]/I"           );
  m_tree->Branch("sect_quality"         , &(m_treePtr->sect_quality         ), "sect_quality[nsector]/I"          );
  m_tree->Branch("sect_combo_dphi"      , &(m_treePtr->sect_combo_dphi      ), "sect_combo_dphi[nsector][11]/I"   );
  m_tree->Branch("sect_combo_deta"      , &(m_treePtr->sect_combo_deta      ), "sect_combo_deta[nsector][11]/I"   );
  m_tree->Branch("sect_combo_isFirst"   , &(m_treePtr->sect_combo_isFirst   ), "sect_combo_isFirst[nsector][11]/I");
  m_tree->Branch("sect_combo_etaBin"    , &(m_treePtr->sect_combo_etaBin    ), "sect_combo_etaBin[nsector][11]/I" );
  m_tree->Branch("sect_combo_hitId"     , &(m_treePtr->sect_combo_hitId     ), "sect_combo_hitId[nsector][11]/I"  );
  m_tree->Branch("sect_combo_frbId"     , &(m_treePtr->sect_combo_frbId     ), "sect_combo_frbId[nsector][11]/I"  );
  m_tree->Branch("sect_combo_phi1"      , &(m_treePtr->sect_combo_phi1      ), "sect_combo_phi1[nsector][11]/I"   );
  m_tree->Branch("sect_combo_eta1"      , &(m_treePtr->sect_combo_eta1      ), "sect_combo_eta1[nsector][11]/I"   );
  m_tree->Branch("sect_combo_stat1"     , &(m_treePtr->sect_combo_stat1     ), "sect_combo_stat1[nsector][11]/I"  );
  m_tree->Branch("sect_combo_phi2"      , &(m_treePtr->sect_combo_phi2      ), "sect_combo_phi2[nsector][11]/I"   );
  m_tree->Branch("sect_combo_eta2"      , &(m_treePtr->sect_combo_eta2      ), "sect_combo_eta2[nsector][11]/I"   );
  m_tree->Branch("sect_combo_stat2"     , &(m_treePtr->sect_combo_stat2     ), "sect_combo_stat2[nsector][11]/I"  );
  m_tree->Branch("sect_combo_s1ring"    , &(m_treePtr->sect_combo_s1ring    ), "sect_combo_s1ring[nsector][11]/I" );
  
  m_tree->Branch("stub_stat"            , &(m_treePtr->stub_stat            ), "stub_stat[nstub]/I"               );
  m_tree->Branch("stub_ring"            , &(m_treePtr->stub_ring            ), "stub_ring[nstub]/I"               );
  m_tree->Branch("stub_cham"            , &(m_treePtr->stub_cham            ), "stub_cham[nstub]/I"               );
  m_tree->Branch("stub_endc"            , &(m_treePtr->stub_endc            ), "stub_endc[nstub]/I"               );
  m_tree->Branch("stub_sect"            , &(m_treePtr->stub_sect            ), "stub_sect[nstub]/I"               );
  m_tree->Branch("stub_frBit"           , &(m_treePtr->stub_frBit           ), "stub_frBit[nstub]/I"              );
  m_tree->Branch("stub_cscid"           , &(m_treePtr->stub_cscid           ), "stub_cscid[nstub]/I"              );

  m_tree->Branch("stub_digi_eta"        , &(m_treePtr->stub_digi_eta        ), "stub_digi_eta[nstub]/F"           );
  m_tree->Branch("stub_digi_phi"        , &(m_treePtr->stub_digi_phi        ), "stub_digi_phi[nstub]/F"           );
  m_tree->Branch("stub_digi_bend"       , &(m_treePtr->stub_digi_bend       ), "stub_digi_bend[nstub]/I"          );
  m_tree->Branch("stub_digi_strip"      , &(m_treePtr->stub_digi_strip      ), "stub_digi_strip[nstub]/I"         );
  m_tree->Branch("stub_digi_keyWG"      , &(m_treePtr->stub_digi_keyWG      ), "stub_digi_keyWG[nstub]/I"         );
  m_tree->Branch("stub_digi_pattern"    , &(m_treePtr->stub_digi_pattern    ), "stub_digi_pattern[nstub]/I"       );
  m_tree->Branch("stub_digi_quality"    , &(m_treePtr->stub_digi_quality    ), "stub_digi_quality[nstub]/I"       );
  m_tree->Branch("stub_digi_badphi"     , &(m_treePtr->stub_digi_badphi     ), "stub_digi_badphi[nstub]/I"        );
  m_tree->Branch("stub_digi_gblPhi"     , &(m_treePtr->stub_digi_gblPhi     ), "stub_digi_gblPhi[nstub]/I"        );
  m_tree->Branch("stub_digi_gblEta"     , &(m_treePtr->stub_digi_gblEta     ), "stub_digi_gblEta[nstub]/I"        );
  
  /*
  m_tree->Branch("sh_eta"                  , &(m_treePtr->sh_eta                  ), "sh_eta[nsh]/F"                          );
  m_tree->Branch("sh_phi"                  , &(m_treePtr->sh_phi                  ), "sh_phi[nsh]/F"                          );
  m_tree->Branch("sh_pdg"                  , &(m_treePtr->sh_pdg                  ), "sh_pdg[nsh]/I"                          );
  m_tree->Branch("sh_stat"                 , &(m_treePtr->sh_stat                 ), "sh_stat[nsh]/I"                         );
  m_tree->Branch("sh_ring"                 , &(m_treePtr->sh_ring                 ), "sh_ring[nsh]/I"                         );
  m_tree->Branch("sh_cham"                 , &(m_treePtr->sh_cham                 ), "sh_cham[nsh]/I"                         );
  m_tree->Branch("sh_layr"                 , &(m_treePtr->sh_layr                 ), "sh_layr[nsh]/I"                         );
  m_tree->Branch("sh_endc"                 , &(m_treePtr->sh_endc                 ), "sh_endc[nsh]/I"                         );
  m_tree->Branch("sh_sect"                 , &(m_treePtr->sh_sect                 ), "sh_sect[nsh]/I"                         );
  m_tree->Branch("sh_frBit"                , &(m_treePtr->sh_frBit                ), "sh_frBit[nsh]/I"                        );
  m_tree->Branch("sh_cscid"                , &(m_treePtr->sh_cscid                ), "sh_cscid[nsh]/I"                        );
  */

}

CSCTFMuonTree* 
FillCSCTFMuonTree::getTreePtr() {
   return m_treePtr;
}

void 
FillCSCTFMuonTree::fill() {
   m_tree->Fill();
}

void 
FillCSCTFMuonTree::finalize() {
   m_file->Write();
   m_file->Close();
}



