
#include "Analyzers/DijetAnalyzer/interface/FillMyTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

FillMyTree::FillMyTree()
{}

FillMyTree::~FillMyTree()
{}

void FillMyTree::init(std::string filename, MyTree* treePtr) {

   m_file = new TFile( filename.c_str(), "RECREATE" );

   m_tree = new TTree( "dijettree", "Reco Dijet Info", 1 );

   m_treePtr = treePtr;

   m_tree->Branch("run"    , &(m_treePtr->run    ), "run/I"           );
   m_tree->Branch("event"  , &(m_treePtr->event  ), "event/I"         );

   m_tree->Branch("nrjet"  ,    &(m_treePtr->nrjet  ),    "nrjet/I"      );
   m_tree->Branch("ngtower",    &(m_treePtr->ngtower),    "ngtower/I"    );

   m_tree->Branch("jet_njtowers"        ,&(m_treePtr->jet_njtowers        ),"jet_njtowers[nrjet]/I");
   m_tree->Branch("genTower_ieta"       ,&(m_treePtr->genTower_ieta       ),"genTower_ieta[ngtower]/I");
   m_tree->Branch("genTower_ietaAbs"    ,&(m_treePtr->genTower_ietaAbs    ),"genTower_ietaAbs[ngtower]/I"); 
   m_tree->Branch("genTower_iphi"       ,&(m_treePtr->genTower_iphi       ),"genTower_iphi[ngtower]/I"); 
   m_tree->Branch("genTower_zside"      ,&(m_treePtr->genTower_zside      ),"genTower_zside[ngtower]/I"); 
   m_tree->Branch("genTower_ncon"       ,&(m_treePtr->genTower_ncon       ),"genTower_ncon[ngtower]/I");
   
   m_tree->Branch("genTower_con_subdet" ,&(m_treePtr->genTower_con_subdet ),"genTower_con_subdet[ngtower][5]");
   m_tree->Branch("genTower_con_ieta"   ,&(m_treePtr->genTower_con_ieta   ),"genTower_con_ieta[ngtower][5]");
   m_tree->Branch("genTower_con_iphi"   ,&(m_treePtr->genTower_con_iphi   ),"genTower_con_iphi[ngtower][5]");
   m_tree->Branch("genTower_con_zside"  ,&(m_treePtr->genTower_con_zside  ),"genTower_con_zside[ngtower][5]");
   m_tree->Branch("genTower_con_depth"  ,&(m_treePtr->genTower_con_depth  ),"genTower_con_depth[ngtower][5]");
   
   m_tree->Branch("jetTower_ieta"       ,&(m_treePtr->jetTower_ieta       ),"jetTower_ieta[nrjet][100]/I");
   m_tree->Branch("jetTower_iphi"       ,&(m_treePtr->jetTower_iphi       ),"jetTower_iphi[nrjet][100]/I");
   m_tree->Branch("jetTower_zside"      ,&(m_treePtr->jetTower_zside      ),"jetTower_zside[nrjet][100]/I");
   m_tree->Branch("jetTower_ietaAbs"    ,&(m_treePtr->jetTower_ietaAbs    ),"jetTower_ietaAbs[nrjet][100]/I");
   m_tree->Branch("jetTower_nrhits"     ,&(m_treePtr->jetTower_nrhits     ),"jetTower_nrhits[nrjet][100]/I");
   m_tree->Branch("jetTower_rhSubdet"   ,&(m_treePtr->jetTower_rhSubdet   ),"jetTower_rhSubdet[nrjet][100][20]/I");
   m_tree->Branch("jetTower_ntslice"    ,&(m_treePtr->jetTower_ntslice    ),"jetTower_ntslice[nrjet][100][20]/I");
   m_tree->Branch("jetTower_depth"      ,&(m_treePtr->jetTower_depth      ),"jetTower_depth[nrjet][100][20]/I");
   m_tree->Branch("jetTower_digi_ADC"   ,&(m_treePtr->jetTower_digi_ADC   ),"jetTower_digi_ADC[nrjet][100][20][10]/I");
   
   m_tree->Branch("jet_p"             ,&(m_treePtr->jet_p              ),"jet_p[nrjet]/F"  );
   m_tree->Branch("jet_pt"            ,&(m_treePtr->jet_pt             ),"jet_pt[nrjet]/F" );
   m_tree->Branch("jet_px"            ,&(m_treePtr->jet_px             ),"jet_px[nrjet]/F" );
   m_tree->Branch("jet_py"            ,&(m_treePtr->jet_py             ),"jet_py[nrjet]/F" );
   m_tree->Branch("jet_pz"            ,&(m_treePtr->jet_pz             ),"jet_pz[nrjet]/F" );
   m_tree->Branch("jet_e"             ,&(m_treePtr->jet_e              ),"jet_e[nrjet]/F"  );
   m_tree->Branch("jet_et"            ,&(m_treePtr->jet_et             ),"jet_et[nrjet]/F" );
   m_tree->Branch("jet_eta"           ,&(m_treePtr->jet_eta            ),"jet_eta[nrjet]/F");
   m_tree->Branch("jet_phi"           ,&(m_treePtr->jet_phi            ),"jet_phi[nrjet]/F");
   
   m_tree->Branch("genTower_eta"      ,&(m_treePtr->genTower_eta       ),"genTower_eta[ngtower]/F");
   m_tree->Branch("genTower_phi"      ,&(m_treePtr->genTower_phi       ),"genTower_phi[ngtower]/F"); 
   m_tree->Branch("genTower_hadEnergy",&(m_treePtr->genTower_hadEnergy ),"genTower_hadEnergy[ngtower]/F"); 
   m_tree->Branch("genTower_hadEt"    ,&(m_treePtr->genTower_hadEt     ),"genTower_hadEt[ngtower]/F");        
   m_tree->Branch("jetTower_eta"      ,&(m_treePtr->jetTower_eta       ),"jetTower_eta[nrjet][100]/F");
   m_tree->Branch("jetTower_phi"      ,&(m_treePtr->jetTower_phi       ),"jetTower_phi[nrjet][100]/F");
   m_tree->Branch("jetTower_hadEnergy",&(m_treePtr->jetTower_hadEnergy ),"jetTower_hadEnergy[nrjet][100]/F");
   m_tree->Branch("jetTower_hadEt"    ,&(m_treePtr->jetTower_hadEt     ),"jetTower_hadEt[nrjet][100]/F");
   m_tree->Branch("jetTower_emEnergy" ,&(m_treePtr->jetTower_emEnergy  ),"jetTower_emEnergy[nrjet][100]/F");
   m_tree->Branch("jetTower_emEt"     ,&(m_treePtr->jetTower_emEt      ),"jetTower_emEt[nrjet][100]/F");
   
   m_tree->Branch("jetTower_digi_fC"  ,&(m_treePtr->jetTower_digi_fC   ),"jetTower_digi_fC[nrjet][100][20][10]/F");

   m_tree->Print();
}

MyTree* 
FillMyTree::getTreePtr() {
   return m_treePtr;
}

void 
FillMyTree::fill() {
   m_tree->Fill();
}

void 
FillMyTree::finalize() {
   m_file->Write();
   m_file->Close();
}

