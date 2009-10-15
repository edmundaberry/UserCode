#include "Analyzers/CSCTFMuonPtAnalyzer/interface/MuonBinnedPlotStorage.h"
#include "TMath.h"
#include <iostream>

MuonBinnedPlotStorage::MuonBinnedPlotStorage(int nHitCombos, int nFRBCombos, 
					     int nPtBins, int nEtaBins, 
					     int nQualityBins,
					     const std::string& dphiFileName,
					     const std::string& detaFileName,
					     const std::string& etaFileName):
  m_nHitCombos(nHitCombos),
  m_nFRBCombos(nFRBCombos),
  m_nPtBins(nPtBins),  
  m_nEtaBins(nEtaBins),
  m_nQualityBins (nQualityBins),
  m_dphiFileName(dphiFileName),
  m_detaFileName(detaFileName),
  m_etaFileName(etaFileName),
  //  m_dphiHists ( std::vector < std::vector < std::vector <std::vector <std::vector< std::vector<TH1F*> > > > > > (nHitCombos+1, std::vector <std::vector< std::vector < std::vector < std::vector <TH1F*> > > > > (nFRBCombos+1, std::vector < std::vector< std::vector< std::vector <TH1F*> > > > (nPtBins +1 , std::vector < std::vector< std::vector <TH1F*> > > ( nEtaBins+1, std::vector < std::vector <TH1F*> > (m_nQualityBins, std::vector <TH1F*> ( 2, new TH1F()))))))),
  m_dphiHists ( std::vector < std::vector < std::vector < std::vector <std::vector <std::vector< std::vector<TH1F*> > > > > > > (nHitCombos+1, std::vector < std::vector <std::vector< std::vector < std::vector < std::vector <TH1F*> > > > > > (nFRBCombos+1, std::vector < std::vector < std::vector< std::vector< std::vector <TH1F*> > > > > (nPtBins +1 , std::vector < std::vector < std::vector< std::vector <TH1F*> > > > ( nEtaBins+1, std::vector < std::vector < std::vector <TH1F*> > > (2, std::vector < std::vector < TH1F* > > (m_nQualityBins, std::vector <TH1F*> ( 2, new TH1F())))))))),
  m_detaHists ( std::vector < std::vector < std::vector < std::vector <std::vector <std::vector< std::vector<TH1F*> > > > > > > (nHitCombos+1, std::vector < std::vector <std::vector< std::vector < std::vector < std::vector <TH1F*> > > > > > (nFRBCombos+1, std::vector < std::vector < std::vector< std::vector< std::vector <TH1F*> > > > > (nPtBins +1 , std::vector < std::vector < std::vector< std::vector <TH1F*> > > > ( nEtaBins+1, std::vector < std::vector < std::vector <TH1F*> > > (2, std::vector < std::vector < TH1F* > > (m_nQualityBins, std::vector <TH1F*> ( 2, new TH1F())))))))),
  m_etaHists ( std::vector < std::vector < std::vector < std::vector < TH1F* > > > > ( nHitCombos+1, std::vector < std::vector < std::vector < TH1F* > > > ( nPtBins + 1, std::vector < std::vector < TH1F* > > ( m_nQualityBins, std::vector < TH1F* > (2, new TH1F())))))

{

  char histName[200];

  for (int iHitCombo = 0; iHitCombo <= m_nHitCombos - 1; ++iHitCombo){
    for (int iPtBin = 0; iPtBin <= m_nPtBins; ++iPtBin){
      for (int iQuality = 1; iQuality < m_nQualityBins; ++iQuality){
	for (int isFirstCombo = 0; isFirstCombo <= 1; ++isFirstCombo){
	  
	  sprintf(histName,"Eta_HitCombo%d_PtBin%d_Quality%d_IsFirstCombo%d", iHitCombo,iPtBin, iQuality, isFirstCombo);	  
	  m_etaHists [iHitCombo][iPtBin][iQuality][isFirstCombo] = new TH1F (histName,"",m_nEtaBins + 1, -0.5, (float) m_nEtaBins + 0.5 );
	  for (int iS1Ring = 0; iS1Ring <= 1; ++iS1Ring){	    
	    for (int iFRBCombo = 0; iFRBCombo <= m_nFRBCombos; ++iFRBCombo){
	      for (int iEtaBin = 0; iEtaBin <= m_nEtaBins; iEtaBin++){
		
		sprintf(histName,"DPhi_HitCombo%d_FDBCombo%d_PtBin%d_EtaBin%d_S1Ring%d_Quality%d_IsFirstCombo%d",
			iHitCombo,iFRBCombo,iPtBin, iEtaBin, iS1Ring, iQuality, isFirstCombo);            
		
		m_dphiHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo] = new TH1F(histName,"",210,-525., 525. );
		
		sprintf(histName,"DEta_HitCombo%d_FDBCombo%d_PtBin%d_EtaBin%d_S1Ring%d_Quality%d_IsFirstCombo%d",
			iHitCombo,iFRBCombo,iPtBin, iEtaBin, iS1Ring, iQuality, isFirstCombo);            
		
		m_detaHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo] = new TH1F(histName,"",60,-30.,30.);
	      }
	    }
	  }
	}
      }
    }
  }
}

MuonBinnedPlotStorage::~MuonBinnedPlotStorage(){

}

void MuonBinnedPlotStorage::fillEtaBin ( int iHitCombo, int iQuality, int iPtBin, int isFirstCombo, int etaBin ){
  m_etaHists [iHitCombo][iPtBin][iQuality][isFirstCombo] -> Fill (etaBin);
}

void MuonBinnedPlotStorage::fillDPhi (int iHitCombo, int iFRBCombo, int iPtBin, int iEtaBin, int iS1Ring, int iQuality, int isFirstCombo, int deltaPhi) {
  m_dphiHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo] -> Fill(deltaPhi);
}

void MuonBinnedPlotStorage::fillDEta (int iHitCombo, int iFRBCombo, int iPtBin, int iEtaBin, int iS1Ring, int iQuality, int isFirstCombo, int deltaEta) {
  m_detaHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo] -> Fill(deltaEta);
}

void MuonBinnedPlotStorage::save (){

  TFile *dphiFile = new TFile(m_dphiFileName.c_str(),"RECREATE");
  TFile *detaFile = new TFile(m_detaFileName.c_str(),"RECREATE");
  TFile *etaFile  = new TFile(m_etaFileName.c_str() ,"RECREATE");
  
  for (int iHitCombo = 0; iHitCombo <= m_nHitCombos - 1; ++iHitCombo){
    for (int iPtBin = 0; iPtBin <= m_nPtBins; ++iPtBin){
      for (int iQuality = 1; iQuality < m_nQualityBins; ++iQuality){
	for (int isFirstCombo = 0; isFirstCombo <= 1; ++isFirstCombo){

	  Double_t eta_entries =  m_etaHists [iHitCombo][iPtBin][iQuality][isFirstCombo] -> GetEntries();

	  etaFile -> cd();
	  if (eta_entries > 0) m_etaHists [iHitCombo][iPtBin][iQuality][isFirstCombo] -> Write();
	  delete m_etaHists [iHitCombo][iPtBin][iQuality][isFirstCombo];

	  for (int iS1Ring = 0; iS1Ring <= 1; ++iS1Ring){	    
	    for (int iEtaBin = 0; iEtaBin <= m_nEtaBins; iEtaBin++){
	      for (int iFRBCombo = 0; iFRBCombo <= m_nFRBCombos; ++iFRBCombo){
		
		Double_t dphi_entries =  m_dphiHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo] -> GetEntries();
		Double_t deta_entries =  m_detaHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo] -> GetEntries();
		
		dphiFile -> cd();
		if (dphi_entries > 0.0) m_dphiHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo]-> Write();
		
		detaFile -> cd();
		if (deta_entries > 0.0) m_detaHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo]-> Write();
		
		delete m_dphiHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo];
		delete m_detaHists[iHitCombo][iFRBCombo][iPtBin][iEtaBin][iS1Ring][iQuality][isFirstCombo];
	      }
	    }
	  }
	}
      }
    }
  }

  dphiFile -> Close();
  detaFile -> Close();
  etaFile  -> Close();
}
