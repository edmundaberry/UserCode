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
   m_tree->Branch("h_psam",     &(m_treePtr->h_psam),     "h_psam [nhit]/I");
   m_tree->Branch("h_size",     &(m_treePtr->h_size),     "h_size [nhit]/I");
   m_tree->Branch("h_id",       &(m_treePtr->h_id),       "h_id   [nhit]/I");
   m_tree->Branch("h_depth",    &(m_treePtr->h_depth),    "h_depth[nhit]/I");
   m_tree->Branch("h_iphi",     &(m_treePtr->h_iphi),     "h_iphi [nhit]/I");
   m_tree->Branch("h_ieta",     &(m_treePtr->h_ieta),     "h_ieta [nhit]/I");
   m_tree->Branch("h_eta",      &(m_treePtr->h_eta),      "h_eta  [nhit]/F");
   m_tree->Branch("h_phi",      &(m_treePtr->h_phi),      "h_phi  [nhit]/F");

   m_tree->Branch("t_spike",    &(m_treePtr->t_spike),    "t_spike[nhit]/I");
   m_tree->Branch("t_ntp",      &(m_treePtr->t_ntp),      "t_ntp  [nhit]/I");
   m_tree->Branch("t_found",    &(m_treePtr->t_found),    "t_found[nhit][2]/I");
   m_tree->Branch("t_size",     &(m_treePtr->t_size),     "t_size [nhit][2]/I");
   m_tree->Branch("t_psam",     &(m_treePtr->t_psam),     "t_psam [nhit][2]/I");
   m_tree->Branch("t_ntpts",    &(m_treePtr->t_ntpts),    "t_ntpts[nhit][2]/I"); 
   m_tree->Branch("t_cET",      &(m_treePtr->t_cET),      "t_cET  [nhit][2][10]/I");

   m_tree->Branch("h_adc",      &(m_treePtr->h_adc),      "h_adc  [nhit][10]/F");
   m_tree->Branch("h_fC",       &(m_treePtr->h_fC),       "h_fC   [nhit][10]/F");
   m_tree->Branch("h_ped",      &(m_treePtr->h_ped),      "h_ped  [nhit][10]/F");
   m_tree->Branch("h_gain",     &(m_treePtr->h_gain),     "h_gain [nhit][10]/F");
   m_tree->Branch("h_capid",    &(m_treePtr->h_capid),    "h_capid[nhit][10]/I");
   m_tree->Branch("h_fiber",    &(m_treePtr->h_fiber),    "h_fiber[nhit][10]/I");
   m_tree->Branch("h_fchan",    &(m_treePtr->h_fchan),    "h_fchan[nhit][10]/I");

   m_tree->Branch("nadc",       &(m_treePtr->nadc),       "nadc/I");

   m_tree->Print();
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

