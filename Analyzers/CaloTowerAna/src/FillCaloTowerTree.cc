#include "Analyzers/CaloTowerAna/interface/FillCaloTowerTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

FillCaloTowerTree::FillCaloTowerTree(){}

FillCaloTowerTree::~FillCaloTowerTree(){}

void FillCaloTowerTree::init(std::string filename, CaloTowerTree* treePtr) {

   m_file = new TFile( filename.c_str(), "RECREATE" );

   m_tree = new TTree( "caloTowerTree", "CaloTower Info", 1 );

   m_treePtr = treePtr;
   
   m_tree->Branch("run"                , &(m_treePtr->run                ), "run/I"                                    );
   m_tree->Branch("event"              , &(m_treePtr->event              ), "event/I"                                  );
   m_tree->Branch("nrjet"              , &(m_treePtr->nrjet              ), "nrjet/I"                                  );
   m_tree->Branch("rjet_nct"           , &(m_treePtr->rjet_nct           ), "rjet_nct[nrjet]/I"                        ); 
														       
   m_tree->Branch("rjet_ct_ndigi"      , &(m_treePtr->rjet_ct_ndigi      ), "rjet_ct_ndigi[nrjet][100]/I"              );
   m_tree->Branch("rjet_ct_ieta"       , &(m_treePtr->rjet_ct_ieta       ), "rjet_ct_ieta[nrjet][100]/I"               );
   m_tree->Branch("rjet_ct_iphi"       , &(m_treePtr->rjet_ct_iphi       ), "rjet_ct_iphi[nrjet][100]/I"               );
   m_tree->Branch("rjet_ct_digi_ieta"  , &(m_treePtr->rjet_ct_digi_ieta  ), "rjet_cd_digi_ieta[nrjet][100][10]/I"      );
   m_tree->Branch("rjet_ct_digi_iphi"  , &(m_treePtr->rjet_ct_digi_iphi  ), "rjet_cd_digi_iphi[nrjet][100][10]/I"      );
   m_tree->Branch("rjet_ct_digi_depth" , &(m_treePtr->rjet_ct_digi_depth ), "rjet_cd_digi_depth[nrjet][100][10]/I"     );
   m_tree->Branch("rjet_ct_digi_size"  , &(m_treePtr->rjet_ct_digi_size  ), "rjet_cd_digi_size[nrjet][100][10]/I"      );
   m_tree->Branch("rjet_ct_digi_ps"    , &(m_treePtr->rjet_ct_digi_ps    ), "rjet_cd_digi_ps[nrjet][100][10]/I"        );   
   m_tree->Branch("rjet_ct_digi_adc"   , &(m_treePtr->rjet_ct_digi_adc   ), "rjet_ct_digi_adc[nrjet][100][10][10]/I"   );
   m_tree->Branch("rjet_ct_digi_capid" , &(m_treePtr->rjet_ct_digi_capid ), "rjet_ct_digi_capid[nrjet][100][10][10]/I" );

   m_tree->Branch("rjet_pt"            , &(m_treePtr->rjet_pt            ), "rjet_pt[nrjet]/F"                         );
   m_tree->Branch("rjet_e"             , &(m_treePtr->rjet_e             ), "rjet_e[nrjet]/F"                          );
   m_tree->Branch("rjet_et"            , &(m_treePtr->rjet_et            ), "rjet_et[nrjet]/F"                         );
   m_tree->Branch("rjet_eta"           , &(m_treePtr->rjet_eta           ), "rjet_eta[nrjet]/F"                        );
   m_tree->Branch("rjet_phi"           , &(m_treePtr->rjet_phi           ), "rjet_phi[nrjet]/F"                        );
   m_tree->Branch("rjet_ct_hadE"       , &(m_treePtr->rjet_ct_hadE       ), "rjet_ct_hadE[nrjet][100]/F"               );
   m_tree->Branch("rjet_ct_outE"       , &(m_treePtr->rjet_ct_outE       ), "rjet_ct_outE[nrjet][100]/F"               );
   m_tree->Branch("rjet_ct_emE"        , &(m_treePtr->rjet_ct_emE        ), "rjet_ct_emE[nrjet][100]/F"                );
   m_tree->Branch("rjet_ct_digi_fC"    , &(m_treePtr->rjet_ct_digi_fC    ), "rjet_ct_digi_fC[nrjet][100][10][10]/F"    );
   m_tree->Branch("rjet_ct_digi_ped"   , &(m_treePtr->rjet_ct_digi_ped   ), "rjet_ct_digi_ped[nrjet][100][10][10]/F"   );
   m_tree->Branch("rjet_ct_digi_gain"  , &(m_treePtr->rjet_ct_digi_gain  ), "rjet_ct_digi_gain[nrjet][100][10][10]/F"  );
   m_tree->Branch("rjet_ct_digi_rcgain", &(m_treePtr->rjet_ct_digi_rcgain), "rjet_ct_digi_rcgain[nrjet][100][10][10]/F");

}

CaloTowerTree* 
FillCaloTowerTree::getTreePtr() {
   return m_treePtr;
}

void 
FillCaloTowerTree::fill() {
   m_tree->Fill();
}

void 
FillCaloTowerTree::finalize() {
   m_file->Write();
   m_file->Close();
}

