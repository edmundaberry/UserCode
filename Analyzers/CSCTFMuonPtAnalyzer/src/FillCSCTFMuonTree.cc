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

  m_tree->Branch("run"                     , &(m_treePtr->run                     ), "run/I"                                  );
  m_tree->Branch("event"                   , &(m_treePtr->event                   ), "event/I"                                );
  m_tree->Branch("filled"                  , &(m_treePtr->filled                  ), "filled/I"                               );
  m_tree->Branch("ngmu"                    , &(m_treePtr->ngmu                    ), "ngmu/I"                                 );
  m_tree->Branch("nl1detid"                , &(m_treePtr->nl1detid                ), "nl1detid/I"                             );
  m_tree->Branch("ncombo"                  , &(m_treePtr->ncombo                  ), "ncombo/I"                               );
  m_tree->Branch("evtComboHitID"           , &(m_treePtr->evtComboHitID           ), "evtComboHitID/I"                        );
  m_tree->Branch("evtComboFrbID"           , &(m_treePtr->evtComboFrbID           ), "evtComboFrbID/I"                        );
  m_tree->Branch("ptBin"                   , &(m_treePtr->ptBin                   ), "ptBin/I"                                );
  m_tree->Branch("nsh"                     , &(m_treePtr->nsh                     ), "nsh/I"                                  );
  m_tree->Branch("nskip"                   , &(m_treePtr->nskip                   ), "nskip/I"                                );
  
  m_tree->Branch("gmu_pt"                  , &(m_treePtr->gmu_pt                  ), "gmu_pt[ngmu]/F"                         );
  m_tree->Branch("gmu_eta"                 , &(m_treePtr->gmu_eta                 ), "gmu_eta[ngmu]/F"                        );
  m_tree->Branch("gmu_phi"                 , &(m_treePtr->gmu_phi                 ), "gmu_phi[ngmu]/F"                        );
  
  m_tree->Branch("combo_phi"               , &(m_treePtr->combo_phi               ), "combo_phi[ncombo][2]/F"                 );
  m_tree->Branch("combo_eta"               , &(m_treePtr->combo_eta               ), "combo_eta[ncombo][2]/F"                 );
  m_tree->Branch("combo_dphi"              , &(m_treePtr->combo_dphi              ), "combo_dphi[ncombo]/I"                   );
  m_tree->Branch("combo_isFirst"           , &(m_treePtr->combo_isFirst           ), "combo_isFirst[ncombo]/I"                );
  m_tree->Branch("combo_etaBin"            , &(m_treePtr->combo_etaBin            ), "combo_etaBin[ncombo]/I"                 );
  m_tree->Branch("combo_sector"            , &(m_treePtr->combo_sector            ), "combo_sector[ncombo]/I"                 );
  m_tree->Branch("combo_hitId"             , &(m_treePtr->combo_hitId             ), "combo_hitId[ncombo]/I"                  );
  m_tree->Branch("combo_frbId"             , &(m_treePtr->combo_frbId             ), "combo_frbId[ncombo]/I"                  );
  				            				     	   				      	       
  m_tree->Branch("l1detid_stat"            , &(m_treePtr->l1detid_stat            ), "l1detid_stat[nl1detid]/I"               );
  m_tree->Branch("l1detid_ring"            , &(m_treePtr->l1detid_ring            ), "l1detid_ring[nl1detid]/I"               );
  m_tree->Branch("l1detid_cham"            , &(m_treePtr->l1detid_cham            ), "l1detid_cham[nl1detid]/I"               );
  m_tree->Branch("l1detid_layr"            , &(m_treePtr->l1detid_layr            ), "l1detid_layr[nl1detid]/I"               );
  m_tree->Branch("l1detid_endc"            , &(m_treePtr->l1detid_endc            ), "l1detid_endc[nl1detid]/I"               );
  m_tree->Branch("l1detid_sect"            , &(m_treePtr->l1detid_sect            ), "l1detid_sect[nl1detid]/I"               );
  m_tree->Branch("l1detid_frBit"           , &(m_treePtr->l1detid_frBit           ), "l1detid_frBit[nl1detid]/I"              );
  m_tree->Branch("l1detid_cscid"           , &(m_treePtr->l1detid_cscid           ), "l1detid_cscid[nl1detid]/I"              );
  m_tree->Branch("l1detid_digi_eta"        , &(m_treePtr->l1detid_digi_eta        ), "l1detid_digi_eta[nl1detid]/F"           );
  m_tree->Branch("l1detid_digi_phi"        , &(m_treePtr->l1detid_digi_phi        ), "l1detid_digi_phi[nl1detid]/F"           );
  m_tree->Branch("l1detid_digi_bend"       , &(m_treePtr->l1detid_digi_bend       ), "l1detid_digi_bend[nl1detid]/I"          );
  m_tree->Branch("l1detid_digi_strip"      , &(m_treePtr->l1detid_digi_strip      ), "l1detid_digi_strip[nl1detid]/I"         );
  m_tree->Branch("l1detid_digi_keyWG"      , &(m_treePtr->l1detid_digi_keyWG      ), "l1detid_digi_keyWG[nl1detid]/I"         );
  m_tree->Branch("l1detid_digi_pattern"    , &(m_treePtr->l1detid_digi_pattern    ), "l1detid_digi_pattern[nl1detid]/I"       );
  m_tree->Branch("l1detid_digi_quality"    , &(m_treePtr->l1detid_digi_quality    ), "l1detid_digi_quality[nl1detid]/I"       );
  m_tree->Branch("l1detid_digi_badphi"     , &(m_treePtr->l1detid_digi_badphi     ), "l1detid_digi_badphi[nl1detid]/I"        );

  m_tree->Branch("l1detid_digi_lclPhi"     , &(m_treePtr->l1detid_digi_lclPhi     ), "l1detid_digi_lclPhi[nl1detid]/I"        );
  m_tree->Branch("l1detid_digi_lclPhiBend" , &(m_treePtr->l1detid_digi_lclPhiBend ), "l1detid_digi_lclPhiBend[nl1detid]/I"    );
  m_tree->Branch("l1detid_digi_gblPhi"     , &(m_treePtr->l1detid_digi_gblPhi     ), "l1detid_digi_gblPhi[nl1detid]/I"        );
  m_tree->Branch("l1detid_digi_gblEta"     , &(m_treePtr->l1detid_digi_gblEta     ), "l1detid_digi_gblEta[nl1detid]/I"        );
  
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



