
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

  m_tree -> Branch("run"                  , &(m_treePtr->run                  ), "run/I");
  m_tree -> Branch("event"                , &(m_treePtr->event                ), "event/I");
  				          				      
  m_tree -> Branch("nGenJets"             , &(m_treePtr->nGenJets             ), "nGenJets/I");
  m_tree -> Branch("nL1GctEtHads"         , &(m_treePtr->nL1GctEtHads         ), "nL1GctEtHads/I");
  m_tree -> Branch("nHLTJetCands"         , &(m_treePtr->nHLTJetCands         ), "nHLTJetCands/I");
  m_tree -> Branch("nHLTCorJetCands"      , &(m_treePtr->nHLTCorJetCands      ), "nHLTCorJetCands/I");
  m_tree -> Branch("nRecoJetCands"        , &(m_treePtr->nRecoJetCands        ), "nRecoJetCands/I");
  m_tree -> Branch("nL1CenJet"            , &(m_treePtr->nL1CenJet            ), "nL1CenJet/I");
  m_tree -> Branch("nL1TauJet"            , &(m_treePtr->nL1TauJet            ), "nL1TauJet/I");
  m_tree -> Branch("nL1ForJet"            , &(m_treePtr->nL1ForJet            ), "nL1ForJet/I");
				          				      
  m_tree -> Branch("l1_HTT100"            , &(m_treePtr->l1_HTT100            ), "l1_HTT100/I");
  m_tree -> Branch("l1_HTT200"            , &(m_treePtr->l1_HTT200            ), "l1_HTT200/I");
  m_tree -> Branch("l1_HTT300"            , &(m_treePtr->l1_HTT300            ), "l1_HTT300/I");
  m_tree -> Branch("l1_HTT400"            , &(m_treePtr->l1_HTT400            ), "l1_HTT400/I");
  m_tree -> Branch("l1_HTT500"            , &(m_treePtr->l1_HTT500            ), "l1_HTT500/I");
  				          				      
  m_tree -> Branch("l1_SingleJet15"       , &(m_treePtr->l1_SingleJet15       ), "l1_SingleJet15/I");
  m_tree -> Branch("l1_SingleJet20"       , &(m_treePtr->l1_SingleJet20       ), "l1_SingleJet20/I");
  m_tree -> Branch("l1_SingleJet30"       , &(m_treePtr->l1_SingleJet30       ), "l1_SingleJet30/I");
  m_tree -> Branch("l1_SingleJet50"       , &(m_treePtr->l1_SingleJet50       ), "l1_SingleJet50/I");
  				          				      
  m_tree -> Branch("nHLTJetTowers"        , &(m_treePtr->nHLTJetTowers        ), "nHLTJetTowers[nHLTJetCands]/I");		 
  m_tree -> Branch("nHLTCorJetTowers"     , &(m_treePtr->nHLTCorJetTowers     ), "nHLTCorJetTowers[nHLTCorJetCands]/I");		 
                                                      			      
  m_tree -> Branch("l1CenJet_ietaSign"    , &(m_treePtr->l1CenJet_ietaSign    ), "l1CenJet_ietaSign[nL1CenJet]/I");
  m_tree -> Branch("l1CenJet_etaIndex"    , &(m_treePtr->l1CenJet_etaIndex    ), "l1CenJet_etaIndex[nL1CenJet]/I");
  m_tree -> Branch("l1CenJet_phiIndex"    , &(m_treePtr->l1CenJet_phiIndex    ), "l1CenJet_phiIndex[nL1CenJet]/I");
  m_tree -> Branch("l1CenJet_capBlock"    , &(m_treePtr->l1CenJet_capBlock    ), "l1CenJet_capBlock[nL1CenJet]/I");
  m_tree -> Branch("l1CenJet_capIndex"    , &(m_treePtr->l1CenJet_capIndex    ), "l1CenJet_capIndex[nL1CenJet]/I");
  m_tree -> Branch("l1CenJet_bx"          , &(m_treePtr->l1CenJet_bx          ), "l1CenJet_bx[nL1CenJet]/I");
  m_tree -> Branch("l1CenJet_rank"        , &(m_treePtr->l1CenJet_rank        ), "l1CenJet_rank[nL1CenJet]/I");
				          				      
  m_tree -> Branch("l1TauJet_ietaSign"    , &(m_treePtr->l1TauJet_ietaSign    ), "l1TauJet_ietaSign[nL1TauJet]/I");
  m_tree -> Branch("l1TauJet_etaIndex"    , &(m_treePtr->l1TauJet_etaIndex    ), "l1TauJet_etaIndex[nL1TauJet]/I");
  m_tree -> Branch("l1TauJet_phiIndex"    , &(m_treePtr->l1TauJet_phiIndex    ), "l1TauJet_phiIndex[nL1TauJet]/I");
  m_tree -> Branch("l1TauJet_capBlock"    , &(m_treePtr->l1TauJet_capBlock    ), "l1TauJet_capBlock[nL1TauJet]/I");
  m_tree -> Branch("l1TauJet_capIndex"    , &(m_treePtr->l1TauJet_capIndex    ), "l1TauJet_capIndex[nL1TauJet]/I");
  m_tree -> Branch("l1TauJet_bx"          , &(m_treePtr->l1TauJet_bx          ), "l1TauJet_bx[nL1TauJet]/I");
  m_tree -> Branch("l1TauJet_rank"        , &(m_treePtr->l1TauJet_rank        ), "l1TauJet_rank[nL1TauJet]/I");
				          				      
  m_tree -> Branch("l1ForJet_ietaSign"    , &(m_treePtr->l1ForJet_ietaSign    ), "l1ForJet_ietaSign[nL1ForJet]/I");
  m_tree -> Branch("l1ForJet_etaIndex"    , &(m_treePtr->l1ForJet_etaIndex    ), "l1ForJet_etaIndex[nL1ForJet]/I");
  m_tree -> Branch("l1ForJet_phiIndex"    , &(m_treePtr->l1ForJet_phiIndex    ), "l1ForJet_phiIndex[nL1ForJet]/I");
  m_tree -> Branch("l1ForJet_capBlock"    , &(m_treePtr->l1ForJet_capBlock    ), "l1ForJet_capBlock[nL1ForJet]/I");
  m_tree -> Branch("l1ForJet_capIndex"    , &(m_treePtr->l1ForJet_capIndex    ), "l1ForJet_capIndex[nL1ForJet]/I");
  m_tree -> Branch("l1ForJet_bx"          , &(m_treePtr->l1ForJet_bx          ), "l1ForJet_bx[nL1ForJet]/I");
  m_tree -> Branch("l1ForJet_rank"        , &(m_treePtr->l1ForJet_rank        ), "l1ForJet_rank[nL1ForJet]/I");

  m_tree -> Branch("hltJetTower_ieta"     , &(m_treePtr->hltJetTower_ieta     ), "hltJetTower_ieta[nHLTJetCands][100]/I");
  m_tree -> Branch("hltJetTower_iphi"     , &(m_treePtr->hltJetTower_iphi     ), "hltJetTower_iphi[nHLTJetCands][100]/I");
  m_tree -> Branch("hltCorJetTower_ieta"  , &(m_treePtr->hltCorJetTower_ieta  ), "hltCorJetTower_ieta[nHLTCorJetCands][100]/I");
  m_tree -> Branch("hltCorJetTower_iphi"  , &(m_treePtr->hltCorJetTower_iphi  ), "hltCorJetTower_iphi[nHLTCorJetCands][100]/I");

  m_tree -> Branch("recoJetTower_ieta"    , &(m_treePtr->recoJetTower_ieta    ), "recoJetTower_ieta[nRecoJetCands][100]/I");
  m_tree -> Branch("recoJetTower_iphi"    , &(m_treePtr->recoJetTower_iphi    ), "recoJetTower_iphi[nRecoJetCands][100]/I");
  				          				      
  m_tree -> Branch("gctHT"                , &(m_treePtr->gctHT                ), "gctHT[nL1GctEtHads]/F");
  m_tree -> Branch("gctHT_UnCorr"         , &(m_treePtr->gctHT_UnCorr         ), "gctHT_UnCorr[nL1GctEtHads]/F");
				          				      
  m_tree -> Branch("genJet_p"             , &(m_treePtr->genJet_p             ), "genJet_p[nGenJets]/F");
  m_tree -> Branch("genJet_px"            , &(m_treePtr->genJet_px            ), "genJet_px[nGenJets]/F");
  m_tree -> Branch("genJet_py"            , &(m_treePtr->genJet_py            ), "genJet_py[nGenJets]/F");
  m_tree -> Branch("genJet_pz"            , &(m_treePtr->genJet_pz            ), "genJet_pz[nGenJets]/F");
  m_tree -> Branch("genJet_pt"            , &(m_treePtr->genJet_pt            ), "genJet_pt[nGenJets]/F");
  m_tree -> Branch("genJet_et"            , &(m_treePtr->genJet_et            ), "genJet_et[nGenJets]/F");
  m_tree -> Branch("genJet_eta"           , &(m_treePtr->genJet_eta           ), "genJet_eta[nGenJets]/F");
  m_tree -> Branch("genJet_phi"           , &(m_treePtr->genJet_phi           ), "genJet_phi[nGenJets]/F");
				          				      
  m_tree -> Branch("l1CenJet_etaMin"      , &(m_treePtr->l1CenJet_etaMin      ), "l1CenJet_etaMin[nL1CenJet]/F");  
  m_tree -> Branch("l1CenJet_etaMax"      , &(m_treePtr->l1CenJet_etaMax      ), "l1CenJet_etaMax[nL1CenJet]/F");
  m_tree -> Branch("l1CenJet_phiMin"      , &(m_treePtr->l1CenJet_phiMin      ), "l1CenJet_phiMin[nL1CenJet]/F");
  m_tree -> Branch("l1CenJet_phiMax"      , &(m_treePtr->l1CenJet_phiMax      ), "l1CenJet_phiMax[nL1CenJet]/F"); 
  m_tree -> Branch("l1CenJet_px"          , &(m_treePtr->l1CenJet_px          ), "l1CenJet_px[nL1CenJet]/F");
  m_tree -> Branch("l1CenJet_py"          , &(m_treePtr->l1CenJet_py          ), "l1CenJet_py[nL1CenJet]/F");
  m_tree -> Branch("l1CenJet_pz"          , &(m_treePtr->l1CenJet_pz          ), "l1CenJet_pz[nL1CenJet]/F");
  m_tree -> Branch("l1CenJet_pt"          , &(m_treePtr->l1CenJet_pt          ), "l1CenJet_pt[nL1CenJet]/F");
  m_tree -> Branch("l1CenJet_et"          , &(m_treePtr->l1CenJet_et          ), "l1CenJet_et[nL1CenJet]/F");
  m_tree -> Branch("l1CenJet_eta"         , &(m_treePtr->l1CenJet_eta         ), "l1CenJet_eta[nL1CenJet]/F");
  m_tree -> Branch("l1CenJet_phi"         , &(m_treePtr->l1CenJet_phi         ), "l1CenJet_phi[nL1CenJet]/F");
  				          				      
  m_tree -> Branch("l1TauJet_etaMin"      , &(m_treePtr->l1TauJet_etaMin      ), "l1TauJet_etaMin[nL1TauJet]/F");  
  m_tree -> Branch("l1TauJet_etaMax"      , &(m_treePtr->l1TauJet_etaMax      ), "l1TauJet_etaMax[nL1TauJet]/F");
  m_tree -> Branch("l1TauJet_phiMin"      , &(m_treePtr->l1TauJet_phiMin      ), "l1TauJet_phiMin[nL1TauJet]/F");
  m_tree -> Branch("l1TauJet_phiMax"      , &(m_treePtr->l1TauJet_phiMax      ), "l1TauJet_phiMax[nL1TauJet]/F"); 
  m_tree -> Branch("l1TauJet_px"          , &(m_treePtr->l1TauJet_px          ), "l1TauJet_px[nL1TauJet]/F");
  m_tree -> Branch("l1TauJet_py"          , &(m_treePtr->l1TauJet_py          ), "l1TauJet_py[nL1TauJet]/F");
  m_tree -> Branch("l1TauJet_pz"          , &(m_treePtr->l1TauJet_pz          ), "l1TauJet_pz[nL1TauJet]/F");
  m_tree -> Branch("l1TauJet_pt"          , &(m_treePtr->l1TauJet_pt          ), "l1TauJet_pt[nL1TauJet]/F");
  m_tree -> Branch("l1TauJet_et"          , &(m_treePtr->l1TauJet_et          ), "l1TauJet_et[nL1TauJet]/F");
  m_tree -> Branch("l1TauJet_eta"         , &(m_treePtr->l1TauJet_eta         ), "l1TauJet_eta[nL1TauJet]/F");
  m_tree -> Branch("l1TauJet_phi"         , &(m_treePtr->l1TauJet_phi         ), "l1TauJet_phi[nL1TauJet]/F");
  				          				      
  m_tree -> Branch("l1ForJet_etaMin"      , &(m_treePtr->l1ForJet_etaMin      ), "l1ForJet_etaMin[nL1ForJet]/F");  
  m_tree -> Branch("l1ForJet_etaMax"      , &(m_treePtr->l1ForJet_etaMax      ), "l1ForJet_etaMax[nL1ForJet]/F");
  m_tree -> Branch("l1ForJet_phiMin"      , &(m_treePtr->l1ForJet_phiMin      ), "l1ForJet_phiMin[nL1ForJet]/F");
  m_tree -> Branch("l1ForJet_phiMax"      , &(m_treePtr->l1ForJet_phiMax      ), "l1ForJet_phiMax[nL1ForJet]/F"); 
  m_tree -> Branch("l1ForJet_px"          , &(m_treePtr->l1ForJet_px          ), "l1ForJet_px[nL1ForJet]/F");
  m_tree -> Branch("l1ForJet_py"          , &(m_treePtr->l1ForJet_py          ), "l1ForJet_py[nL1ForJet]/F");
  m_tree -> Branch("l1ForJet_pz"          , &(m_treePtr->l1ForJet_pz          ), "l1ForJet_pz[nL1ForJet]/F");
  m_tree -> Branch("l1ForJet_pt"          , &(m_treePtr->l1ForJet_pt          ), "l1ForJet_pt[nL1ForJet]/F");
  m_tree -> Branch("l1ForJet_et"          , &(m_treePtr->l1ForJet_et          ), "l1ForJet_et[nL1ForJet]/F");
  m_tree -> Branch("l1ForJet_eta"         , &(m_treePtr->l1ForJet_eta         ), "l1ForJet_eta[nL1ForJet]/F");
  m_tree -> Branch("l1ForJet_phi"         , &(m_treePtr->l1ForJet_phi         ), "l1ForJet_phi[nL1ForJet]/F");
				          				      
  m_tree -> Branch("hltHT"                , &(m_treePtr->hltHT                ), "hltHT/F");
  m_tree -> Branch("hltJet_pt"            , &(m_treePtr->hltJet_pt            ), "hltJet_pt[nHLTJetCands]/F");		 
  m_tree -> Branch("hltJet_et"            , &(m_treePtr->hltJet_et            ), "hltJet_et[nHLTJetCands]/F"); 		 
  m_tree -> Branch("hltJet_phi"           , &(m_treePtr->hltJet_phi           ), "hltJet_phi[nHLTJetCands]/F");		 
  m_tree -> Branch("hltJet_eta"           , &(m_treePtr->hltJet_eta           ), "hltJet_eta[nHLTJetCands]/F");		 

  m_tree -> Branch("hltCorHT"             , &(m_treePtr->hltCorHT             ), "hltCorHT/F");
  m_tree -> Branch("hltCorJet_pt"         , &(m_treePtr->hltCorJet_pt         ), "hltCorJet_pt[nHLTCorJetCands]/F");		 
  m_tree -> Branch("hltCorJet_et"         , &(m_treePtr->hltCorJet_et         ), "hltCorJet_et[nHLTCorJetCands]/F"); 		 
  m_tree -> Branch("hltCorJet_phi"        , &(m_treePtr->hltCorJet_phi        ), "hltCorJet_phi[nHLTCorJetCands]/F");		 
  m_tree -> Branch("hltCorJet_eta"        , &(m_treePtr->hltCorJet_eta        ), "hltCorJet_eta[nHLTCorJetCands]/F");		 
				          				      
  m_tree -> Branch("recoHT"               , &(m_treePtr->recoHT               ), "recoHT/F");
  m_tree -> Branch("recoJet_pt"           , &(m_treePtr->recoJet_pt           ), "recoJet_pt[nRecoJetCands]/F");		 
  m_tree -> Branch("recoJet_et"           , &(m_treePtr->recoJet_et           ), "recoJet_et[nRecoJetCands]/F"); 		 
  m_tree -> Branch("recoJet_phi"          , &(m_treePtr->recoJet_phi          ), "recoJet_phi[nRecoJetCands]/F");		 
  m_tree -> Branch("recoJet_eta"          , &(m_treePtr->recoJet_eta          ), "recoJet_eta[nRecoJetCands]/F");		 

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
