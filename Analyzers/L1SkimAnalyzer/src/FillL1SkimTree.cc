
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

  m_tree -> Branch("caloEm_ieta"       , &(m_treePtr->caloEm_ieta       ), "caloEm_ieta[nCaloEm]/I");      
  m_tree -> Branch("caloEm_iphi"       , &(m_treePtr->caloEm_iphi       ), "caloEm_iphi[nCaloEm]/I");       
  m_tree -> Branch("caloEm_rank"       , &(m_treePtr->caloEm_rank       ), "caloEm_rank[nCaloEm]/I");       
  m_tree -> Branch("caloEm_raw"        , &(m_treePtr->caloEm_raw        ), "caloEm_raw[nCaloEm]/I");        
  m_tree -> Branch("caloEm_rctCrate"   , &(m_treePtr->caloEm_rctCrate   ), "caloEm_rctCrate[nCaloEm]/I");
  m_tree -> Branch("caloEm_rctRegion"  , &(m_treePtr->caloEm_rctRegion  ), "caloEm_rctRegion[nCaloEm]/I");  
  m_tree -> Branch("caloEm_rctCard"    , &(m_treePtr->caloEm_rctCard    ), "caloEm_rctCard[nCaloEm]/I");    
  m_tree -> Branch("caloEm_cindex"     , &(m_treePtr->caloEm_cindex     ), "caloEm_cindex[nCaloEm]/I");     
  m_tree -> Branch("caloEm_bx"         , &(m_treePtr->caloEm_bx         ), "caloEm_bx[nCaloEm]/I");         
  m_tree -> Branch("caloEm_iso"        , &(m_treePtr->caloEm_iso        ), "caloEm_iso[nCaloEm]/I");      

  m_tree -> Branch("nGctCenJet", &(m_treePtr->nGctCenJet),"nGctCenJet/I");
  
  m_tree -> Branch("gctCenJet_ietaSign", &(m_treePtr->gctCenJet_ietaSign), "gctCenJet_ietaSign[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_ieta"    , &(m_treePtr->gctCenJet_ieta    ), "gctCenJet_ieta[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_iphi"    , &(m_treePtr->gctCenJet_iphi    ), "gctCenJet_iphi[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_capBlock", &(m_treePtr->gctCenJet_capBlock), "gctCenJet_capBlock[nGctCenJet]/I");
  m_tree -> Branch("gctCenJet_capIndex", &(m_treePtr->gctCenJet_capIndex), "gctCenJet_capIndex[nGctCenJet]/I");
  m_tree -> Branch("gctCenJet_bx"      , &(m_treePtr->gctCenJet_bx     	), "gctCenJet_bx[nGctCenJet]/I");   
  m_tree -> Branch("gctCenJet_rank"    , &(m_treePtr->gctCenJet_rank   	), "gctCenJet_rank[nGctCenJet]/I");      
  m_tree -> Branch("gctCenJet_et"      , &(m_treePtr->gctCenJet_et      ), "gctCenJet_et[nGctCenJet]/F");     

  m_tree -> Branch("nGctTauJet", &(m_treePtr->nGctTauJet),"nGctTauJet/I");
  
  m_tree -> Branch("gctTauJet_ietaSign", &(m_treePtr->gctTauJet_ietaSign), "gctTauJet_ietaSign[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_ieta"    , &(m_treePtr->gctTauJet_ieta    ), "gctTauJet_ieta[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_iphi"    , &(m_treePtr->gctTauJet_iphi    ), "gctTauJet_iphi[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_capBlock", &(m_treePtr->gctTauJet_capBlock), "gctTauJet_capBlock[nGctTauJet]/I");
  m_tree -> Branch("gctTauJet_capIndex", &(m_treePtr->gctTauJet_capIndex), "gctTauJet_capIndex[nGctTauJet]/I");
  m_tree -> Branch("gctTauJet_bx"      , &(m_treePtr->gctTauJet_bx     	), "gctTauJet_bx[nGctTauJet]/I");   
  m_tree -> Branch("gctTauJet_rank"    , &(m_treePtr->gctTauJet_rank   	), "gctTauJet_rank[nGctTauJet]/I");      
  m_tree -> Branch("gctTauJet_et"      , &(m_treePtr->gctTauJet_et      ), "gctTauJet_et[nGctTauJet]/F");     

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
