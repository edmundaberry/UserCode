#include "Analyzers/HcalDigiAnalyzer/interface/FillDigiTree.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TROOT.h"

FillDigiTree::FillDigiTree(){}

FillDigiTree::~FillDigiTree(){}

void 
FillDigiTree::init(std::string filename, DigiTree* treePtr) {
  
   m_treePtr = treePtr;
   
   m_file = new TFile( filename.c_str(), "RECREATE" );
   
   m_tree = new TTree( "dqmtree", "HLT Digi Info", 1 );   
   
   m_tree->Branch("run",        &(m_treePtr->run),        "run/I");
   m_tree->Branch("event",      &(m_treePtr->event),      "event/I");
   m_tree->Branch("nint",       &(m_treePtr->nint),       "nint/I");
   m_tree->Branch("nchn",       &(m_treePtr->nchn),       "nchn/I");
   m_tree->Branch("nhit",       &(m_treePtr->nhit),       "nhit/I");

   m_tree->Branch("h_id",       &(m_treePtr->h_id),       "h_id    [nhit]/I");
   m_tree->Branch("h_depth",    &(m_treePtr->h_depth),    "h_depth [nhit]/I");
   m_tree->Branch("h_iphi",     &(m_treePtr->h_iphi),     "h_iphi  [nhit]/I");
   m_tree->Branch("h_ieta",     &(m_treePtr->h_ieta),     "h_ieta  [nhit]/I");
   
   m_tree->Branch("h_psam",     &(m_treePtr->h_psam),     "h_psam  [83][73][4]/I");
   m_tree->Branch("h_size",     &(m_treePtr->h_size),     "h_size  [83][73][4]/I");
   m_tree->Branch("h_adc",      &(m_treePtr->h_adc),      "h_adc   [83][73][4][10]/I");
   m_tree->Branch("h_capid",    &(m_treePtr->h_capid),    "h_capid [83][73][4][10]/I");
   m_tree->Branch("h_fiber",    &(m_treePtr->h_fiber),    "h_fiber [83][73][4][10]/I");
   m_tree->Branch("h_fchan",    &(m_treePtr->h_fchan),    "h_fchan [83][73][4][10]/I");
   m_tree->Branch("h_fC",       &(m_treePtr->h_fC),       "h_fC    [83][73][4][10]/F");
   m_tree->Branch("h_ped",      &(m_treePtr->h_ped),      "h_ped   [83][73][4][10]/F");
   m_tree->Branch("h_pedc",     &(m_treePtr->h_pedc),     "h_pedc  [83][73][4][10]/F");
   m_tree->Branch("h_gain",     &(m_treePtr->h_gain),     "h_gain  [83][73][4][10]/F");

   m_tree->Branch("t_spike",    &(m_treePtr->t_spike),    "t_spike [83][73][4]/I");
   m_tree->Branch("t_ntp",      &(m_treePtr->t_ntp),      "t_ntp   [83][73][4]/I");
   m_tree->Branch("t_found",    &(m_treePtr->t_found),    "t_found [83][73][4][2]/I");
   m_tree->Branch("t_size",     &(m_treePtr->t_size),     "t_size  [83][73][4][2]/I");
   m_tree->Branch("t_psam",     &(m_treePtr->t_psam),     "t_psam  [83][73][4][2]/I");
   m_tree->Branch("t_ntpts",    &(m_treePtr->t_ntpts),    "t_ntpts [83][73][4][2]/I"); 
   m_tree->Branch("t_ntpts",    &(m_treePtr->t_ntpts),    "t_ntpts [83][73][4][2]/I"); 
   m_tree->Branch("t_ieta",     &(m_treePtr->t_ieta),     "t_ieta  [83][73][4][2]/I"); 
   m_tree->Branch("t_iphi",     &(m_treePtr->t_iphi),     "t_iphi  [83][73][4][2]/I"); 
   m_tree->Branch("t_subdet",   &(m_treePtr->t_subdet),   "t_subdet[83][73][4][2]/I"); 
   m_tree->Branch("t_cET",      &(m_treePtr->t_cET),      "t_cET   [83][73][4][2][10]/I");

   //m_tree->Print();
}

DigiTree* FillDigiTree::getTreePtr() {
  return m_treePtr;
}

void FillDigiTree::fill() {
   m_tree->Fill();
}

void FillDigiTree::finalize() {
   m_file->Write();
   m_file->Close();
}

