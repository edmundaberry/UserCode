
#include "Analyzers/L1SkimAnalyzer/interface/FillL1SkimTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

FillL1SkimTree::FillL1SkimTree(){}

FillL1SkimTree::~FillL1SkimTree(){}

void FillL1SkimTree:: init(std::string filename, L1SkimTree* treePtr){

  m_treePtr = treePtr;

  m_file = new TFile( filename.c_str(), "RECREATE" );

  m_tree = new TTree("l1SkimTree","L1 Skim Info", 1);

  
  m_tree -> Branch("nCaloEm", &(m_treePtr->nCaloEm), "nCaloEm/I");

  m_tree -> Branch("caloEm_ieta"      , &(m_treePtr->caloEm_ieta     ), "caloEm_ieta[nCaloEm]/I");      
  m_tree -> Branch("caloEm_iphi"      , &(m_treePtr->caloEm_iphi     ), "caloEm_iphi[nCaloEm]/I");       
  m_tree -> Branch("caloEm_rank"      , &(m_treePtr->caloEm_rank     ), "caloEm_rank[nCaloEm]/I");       
  m_tree -> Branch("caloEm_raw"       , &(m_treePtr->caloEm_raw      ), "caloEm_raw[nCaloEm]/I");        
  m_tree -> Branch("caloEm_rctCrate"  , &(m_treePtr->caloEm_rctCrate ), "caloEm_rctCrate[nCaloEm]/I");
  m_tree -> Branch("caloEm_rctRegion" , &(m_treePtr->caloEm_rctRegion), "caloEm_rctRegion[nCaloEm]/I");  
  m_tree -> Branch("caloEm_rctCard"   , &(m_treePtr->caloEm_rctCard  ), "caloEm_rctCard[nCaloEm]/I");    
  m_tree -> Branch("caloEm_cindex"    , &(m_treePtr->caloEm_cindex   ), "caloEm_cindex[nCaloEm]/I");     
  m_tree -> Branch("caloEm_bx"        , &(m_treePtr->caloEm_bx       ), "caloEm_bx[nCaloEm]/I");         
  m_tree -> Branch("caloEm_iso"       , &(m_treePtr->caloEm_iso      ), "caloEm_iso[nCaloEm]/I");        

}

L1SkimTree* FillL1SkimTree::getTreePtr(){
  return m_treePtr;
}

void FillL1SkimTree::fill(){
  m_tree->Fill();
}

void FillL1SkimTree::finalize(){
   m_file->Write();
   m_file->Close();
}
