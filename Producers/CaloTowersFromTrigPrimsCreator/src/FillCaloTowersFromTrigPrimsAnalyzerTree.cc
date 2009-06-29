#include "Producers/CaloTowersFromTrigPrimsCreator/interface/FillCaloTowersFromTrigPrimsAnalyzerTree.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

FillCaloTowersFromTrigPrimsAnalyzerTree::FillCaloTowersFromTrigPrimsAnalyzerTree(){}

FillCaloTowersFromTrigPrimsAnalyzerTree::~FillCaloTowersFromTrigPrimsAnalyzerTree(){}

void FillCaloTowersFromTrigPrimsAnalyzerTree::init(std::string filename, CaloTowersFromTrigPrimsAnalyzerTree* treePtr) {

   m_file = new TFile( filename.c_str(), "RECREATE" );
   m_tree = new TTree( "CaloTowersFromTrigPrimsAnalyzerTree", "CaloTower and TPG Energy and Location Info");

   m_treePtr = treePtr;
   
   m_tree->Branch("run"         , &(m_treePtr->run         ), "run/I"               );
   m_tree->Branch("event"       , &(m_treePtr->event       ), "event/I"             );
   m_tree->Branch("ntpg"        , &(m_treePtr->ntpg        ), "ntpg/I"              );
   m_tree->Branch("nct"         , &(m_treePtr->nct         ), "nct/I"               );
   
   m_tree->Branch("tpg_ieta"    , &(m_treePtr->tpg_ieta    ), "tpg_ieta[ntpg]/I"    );
   m_tree->Branch("tpg_iphi"    , &(m_treePtr->tpg_iphi    ), "tpg_iphi[ntpg]/I"    );
   						    
   m_tree->Branch("tpg_isHcal"  , &(m_treePtr->tpg_isHcal  ), "tpg_isHcal[ntpg]/I"  );
   m_tree->Branch("tpg_isEcal"  , &(m_treePtr->tpg_isEcal  ), "tpg_isEcal[ntpg]/I"  );
   m_tree->Branch("tpg_isHF"    , &(m_treePtr->tpg_isHF    ), "tpg_isHF[ntpg]/I"    );

   m_tree->Branch("ct_ieta"     , &(m_treePtr->ct_ieta     ), "ct_ieta[nct]/I"      );
   m_tree->Branch("ct_iphi"     , &(m_treePtr->ct_iphi     ), "ct_iphi[nct]/I"      );
   m_tree->Branch("ct_nhcon"    , &(m_treePtr->ct_nhcon    ), "ct_nhcon[nct]/I"     );
   m_tree->Branch("ct_necon"    , &(m_treePtr->ct_necon    ), "ct_necon[nct]/I"     );
   
   m_tree->Branch("ct_isMine"   , &(m_treePtr->ct_isMine   ), "ct_isMine[nct]/I"    );

   m_tree->Branch("tpg_et"      , &(m_treePtr->tpg_et      ), "tpg_et[ntpg]/F"      );
   m_tree->Branch("tpg_energy"  , &(m_treePtr->tpg_energy  ), "tpg_energy[ntpg]/F"  );
   m_tree->Branch("tpg_meanEta" , &(m_treePtr->tpg_meanEta ), "tpg_meanEta[ntpg]/F" );
			        		           		            
   m_tree->Branch("ct_eEm"      , &(m_treePtr->ct_eEm      ), "ct_eEm[nct]/F"       );
   m_tree->Branch("ct_eHad"     , &(m_treePtr->ct_eHad     ), "ct_eHad[nct]/F"      );
   m_tree->Branch("ct_eOut"     , &(m_treePtr->ct_eOut     ), "ct_eOut[nct]/F"      );
			        		           		            
   m_tree->Branch("ct_etEm"     , &(m_treePtr->ct_etEm     ), "ct_etEm[nct]/F"      );
   m_tree->Branch("ct_etHad"    , &(m_treePtr->ct_etHad    ), "ct_etHad[nct]/F"     );
   m_tree->Branch("ct_etOut"    , &(m_treePtr->ct_etOut    ), "ct_etOut[nct]/F"     );
   
}

CaloTowersFromTrigPrimsAnalyzerTree* 
FillCaloTowersFromTrigPrimsAnalyzerTree::getTreePtr() {
   return m_treePtr;
}

void 
FillCaloTowersFromTrigPrimsAnalyzerTree::fill() {
  m_tree->Fill();
}

void 
FillCaloTowersFromTrigPrimsAnalyzerTree::finalize() {
   m_file->Write();
   m_file->Close();
}

