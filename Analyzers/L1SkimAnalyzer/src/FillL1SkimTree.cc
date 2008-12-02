
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

  m_tree -> Branch("nCaloEm"           , &(m_treePtr->nCaloEm           ), "nCaloEm/I");
  m_tree -> Branch("nHLTJetCands"      , &(m_treePtr->nHLTJetCands      ), "nHLTJetCands/I");
  m_tree -> Branch("nGctCenJet"        , &(m_treePtr->nGctCenJet        ), "nGctCenJet/I");
  m_tree -> Branch("nGctForJet"        , &(m_treePtr->nGctForJet        ), "nGctForJet/I");
  m_tree -> Branch("nGctTauJet"        , &(m_treePtr->nGctTauJet        ), "nGctTauJet/I");
  m_tree -> Branch("nL1GctEtHads"      , &(m_treePtr->nL1GctEtHads      ), "nL1GctEtHads/I");

  m_tree -> Branch("l1_HTT100"         , &(m_treePtr->l1_HTT100         ), "l1_HTT100/I");
  m_tree -> Branch("l1_HTT200"         , &(m_treePtr->l1_HTT200         ), "l1_HTT200/I");
  m_tree -> Branch("l1_HTT300"         , &(m_treePtr->l1_HTT300         ), "l1_HTT300/I");
  m_tree -> Branch("l1_HTT400"         , &(m_treePtr->l1_HTT400         ), "l1_HTT400/I");
  m_tree -> Branch("l1_HTT500"         , &(m_treePtr->l1_HTT500         ), "l1_HTT500/I");

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

  m_tree -> Branch("gctCenJet_iphiMin" , &(m_treePtr->gctCenJet_iphiMin ), "gctCenJet_iphiMin[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_iphiMax" , &(m_treePtr->gctCenJet_iphiMax ), "gctCenJet_iphiMax[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_ietaMin" , &(m_treePtr->gctCenJet_ietaMin ), "gctCenJet_ietaMin[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_ietaMax" , &(m_treePtr->gctCenJet_ietaMax ), "gctCenJet_ietaMax[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_ietaSign", &(m_treePtr->gctCenJet_ietaSign), "gctCenJet_ietaSign[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_etaIndex", &(m_treePtr->gctCenJet_etaIndex), "gctCenJet_etaIndex[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_phiIndex", &(m_treePtr->gctCenJet_phiIndex), "gctCenJet_phiIndex[nGctCenJet]/I");    
  m_tree -> Branch("gctCenJet_capBlock", &(m_treePtr->gctCenJet_capBlock), "gctCenJet_capBlock[nGctCenJet]/I");
  m_tree -> Branch("gctCenJet_capIndex", &(m_treePtr->gctCenJet_capIndex), "gctCenJet_capIndex[nGctCenJet]/I");
  m_tree -> Branch("gctCenJet_bx"      , &(m_treePtr->gctCenJet_bx     	), "gctCenJet_bx[nGctCenJet]/I");   
  m_tree -> Branch("gctCenJet_rank"    , &(m_treePtr->gctCenJet_rank   	), "gctCenJet_rank[nGctCenJet]/I");      
  
  m_tree -> Branch("gctForJet_iphiMin" , &(m_treePtr->gctForJet_iphiMin ), "gctForJet_iphiMin[nGctForJet]/I");    
  m_tree -> Branch("gctForJet_iphiMax" , &(m_treePtr->gctForJet_iphiMax ), "gctForJet_iphiMax[nGctForJet]/I");    
  m_tree -> Branch("gctForJet_ietaMin" , &(m_treePtr->gctForJet_ietaMin ), "gctForJet_ietaMin[nGctForJet]/I");    
  m_tree -> Branch("gctForJet_ietaMax" , &(m_treePtr->gctForJet_ietaMax ), "gctForJet_ietaMax[nGctForJet]/I");    
  m_tree -> Branch("gctForJet_ietaSign", &(m_treePtr->gctForJet_ietaSign), "gctForJet_ietaSign[nGctForJet]/I");    
  m_tree -> Branch("gctForJet_etaIndex", &(m_treePtr->gctForJet_etaIndex), "gctForJet_etaIndex[nGctForJet]/I");    
  m_tree -> Branch("gctForJet_phiIndex", &(m_treePtr->gctForJet_phiIndex), "gctForJet_phiIndex[nGctForJet]/I");    
  m_tree -> Branch("gctForJet_capBlock", &(m_treePtr->gctForJet_capBlock), "gctForJet_capBlock[nGctForJet]/I");
  m_tree -> Branch("gctForJet_capIndex", &(m_treePtr->gctForJet_capIndex), "gctForJet_capIndex[nGctForJet]/I");
  m_tree -> Branch("gctForJet_bx"      , &(m_treePtr->gctForJet_bx     	), "gctForJet_bx[nGctForJet]/I");   
  m_tree -> Branch("gctForJet_rank"    , &(m_treePtr->gctForJet_rank   	), "gctForJet_rank[nGctForJet]/I");      
  
  m_tree -> Branch("gctTauJet_iphiMin" , &(m_treePtr->gctTauJet_iphiMin ), "gctTauJet_iphiMin[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_iphiMax" , &(m_treePtr->gctTauJet_iphiMax ), "gctTauJet_iphiMax[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_ietaMin" , &(m_treePtr->gctTauJet_ietaMin ), "gctTauJet_ietaMin[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_ietaMax" , &(m_treePtr->gctTauJet_ietaMax ), "gctTauJet_ietaMax[nGctTauJet]/I");     
  m_tree -> Branch("gctTauJet_ietaSign", &(m_treePtr->gctTauJet_ietaSign), "gctTauJet_ietaSign[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_etaIndex", &(m_treePtr->gctTauJet_etaIndex), "gctTauJet_etaIndex[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_phiIndex", &(m_treePtr->gctTauJet_phiIndex), "gctTauJet_phiIndex[nGctTauJet]/I");    
  m_tree -> Branch("gctTauJet_capBlock", &(m_treePtr->gctTauJet_capBlock), "gctTauJet_capBlock[nGctTauJet]/I");
  m_tree -> Branch("gctTauJet_capIndex", &(m_treePtr->gctTauJet_capIndex), "gctTauJet_capIndex[nGctTauJet]/I");
  m_tree -> Branch("gctTauJet_bx"      , &(m_treePtr->gctTauJet_bx     	), "gctTauJet_bx[nGctTauJet]/I");   
  m_tree -> Branch("gctTauJet_rank"    , &(m_treePtr->gctTauJet_rank   	), "gctTauJet_rank[nGctTauJet]/I");      

  m_tree -> Branch("nHLTJetTowers"     , &(m_treePtr->nHLTJetTowers     ), "nHLTJetTowers[nHLTJetCands]/I");		 
                                                   
  m_tree -> Branch("hltJetTower_ieta"  , &(m_treePtr->hltJetTower_ieta  ), "hltJetTower_ieta[nHLTJetCands][100]/I");
  m_tree -> Branch("hltJetTower_iphi"  , &(m_treePtr->hltJetTower_iphi  ), "hltJetTower_iphi[nHLTJetCands][100]/I");

  m_tree -> Branch("gctCenJet_et"      , &(m_treePtr->gctCenJet_et      ), "gctCenJet_et[nGctCenJet]/F");       
  m_tree -> Branch("gctForJet_et"      , &(m_treePtr->gctForJet_et      ), "gctForJet_et[nGctForJet]/F");       
  m_tree -> Branch("gctTauJet_et"      , &(m_treePtr->gctTauJet_et      ), "gctTauJet_et[nGctTauJet]/F");     
  
  m_tree -> Branch("gctHT"             , &(m_treePtr->gctHT             ), "gctHT[nL1GctEtHads]/F");
  m_tree -> Branch("gctHT_UnCorr"      , &(m_treePtr->gctHT_UnCorr      ), "gctHT_UnCorr[nL1GctEtHads]/F");

  m_tree -> Branch("hltJet_pt"         , &(m_treePtr->hltJet_pt         ), "hltJet_pt[nHLTJetCands]/F");		 
  m_tree -> Branch("hltJet_et"         , &(m_treePtr->hltJet_et         ), "hltJet_et[nHLTJetCands]/F"); 		 
  m_tree -> Branch("hltJet_phi"        , &(m_treePtr->hltJet_phi        ), "hltJet_phi[nHLTJetCands]/F");		 
  m_tree -> Branch("hltJet_eta"        , &(m_treePtr->hltJet_eta        ), "hltJet_eta[nHLTJetCands]/F");		 

}

L1SkimTree* FillL1SkimTree::getTreePtr(){
  return m_treePtr;
}

void FillL1SkimTree::fill(){
  m_tree->Fill();
}

void FillL1SkimTree::finalize(){
  
  m_file->Write(); 
  
  //m_file->Close();

}
