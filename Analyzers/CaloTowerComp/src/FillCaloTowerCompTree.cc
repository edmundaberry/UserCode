#include "Analyzers/CaloTowerComp/interface/FillCaloTowerCompTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

FillCaloTowerCompTree::FillCaloTowerCompTree(){}

FillCaloTowerCompTree::~FillCaloTowerCompTree(){}

void FillCaloTowerCompTree::init(std::string filename, CaloTowerCompTree* treePtr) {

   m_file = new TFile( filename.c_str(), "RECREATE" );
   m_tree = new TTree( "caloTowerCompTree", "CaloTower Comparison Info", 1 );

   m_treePtr = treePtr;
   
   m_tree->Branch("run"      , &(m_treePtr->run     ), "run/I"          );
   m_tree->Branch("event"    , &(m_treePtr->event   ), "event/I"        );
   m_tree->Branch("nct"      , &(m_treePtr->nct     ), "nct/I"          );

   m_tree->Branch("ct_ieta"  , &(m_treePtr->ct_ieta ), "ct_ieta[nct]/I" );
   m_tree->Branch("ct_iphi"  , &(m_treePtr->ct_iphi ), "ct_iphi[nct]/I" );
   m_tree->Branch("ct_nhcon" , &(m_treePtr->ct_nhcon), "ct_nhcon[nct]/I");
   m_tree->Branch("ct_necon" , &(m_treePtr->ct_necon), "ct_necon[nct]/I");

   m_tree->Branch("ct_eEm"   , &(m_treePtr->ct_eEm  ), "ct_eEm[nct]/F"  );
   m_tree->Branch("ct_eHad"  , &(m_treePtr->ct_eHad ), "ct_eHad[nct]/F" );
   m_tree->Branch("ct_eOut"  , &(m_treePtr->ct_eOut ), "ct_eOut[nct]/F" );
   
}

CaloTowerCompTree* 
FillCaloTowerCompTree::getTreePtr() {
   return m_treePtr;
}

void 
FillCaloTowerCompTree::fill() {
   m_tree->Fill();
}

void 
FillCaloTowerCompTree::finalize() {
   m_file->Write();
   m_file->Close();
}

