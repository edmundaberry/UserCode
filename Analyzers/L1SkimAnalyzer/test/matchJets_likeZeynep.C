#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TLeafF.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TList.h"
#include "TMatrixD.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TMath.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include <map>
#include <sstream>

//-------------------------------------------------------------
// These parameters must be set for the script to run
//-------------------------------------------------------------

const float matchDeltaR_cut = 0.7;
const int nEventsToExaminePerFile = -1;// 1000;
const float minHltJetPt = 25.0;

//-------------------------------------------------------------
// These parameters should not be changed
//-------------------------------------------------------------

const int nJobsPerThreshold = 10;
const int nEventsPerJob = 10000;

const int minJetThreshold = 3;
const int maxJetThreshold = 30;
const int anaJetThreshold = 15;

const int minPtBin = 6;
const int maxPtBin = 6;
const int nPtBins  = maxPtBin - minPtBin + 1;

const int maxNL1CenJet = 12;
const int maxNL1ForJet = 12;
const int maxNL1TauJet = 12;
const int maxNGctHt    = 12;

const int maxNHltJet   = 100;

const float gctLSB = 1.0;// 0.25;

const double ptBinBounds[14] = {
  0.0  , 15.0 , 20.0 , 30.0 , 50.0 , 80.0 , 120.0, 
  170.0, 230.0, 300.0, 380.0, 470.0, 600.0, 800.0
};

const double ptBinXSecInPb[14] = { 
  0.0,
  51562800000.0, 949441000.0, 400982000.0, 94702500.0, 12195900.0, 1617240.0,
  255987.0, 48325.0, 10623.2, 2634.94, 722.099, 240.983, 62.4923
};

const double lumiInInvPb = 0.0000008;

//-------------------------------------------------------------
// Declare histograms
//-------------------------------------------------------------

TH1F h_l1leaf_size("l1_leaf_size","l1 jets sizes",20,0,20);
TH1F h_rcleaf_size("rc_leaf_size","reco jets sizes",50,0,50);  

TH1F recoPt("recoPt","Pt of lead reco jet", 100,0,250);
TH1F recoPtinforward("recoPtinforward","Pt of lead reco jet in forward region", 70,0,250);
TH1F forrecoPtmatchedl1forward("forrecoPtmatchedl1forward","Pt of for. lead reco jet that matched L1 forward jets", 70,0,250);
TH1F forrecoEtamatchedl1forward("forrecoEtamatchedl1forward","Eta of for. lead reco jet that matched L1 forward jets", 20,-6,6);
TH1F forrecoPhimatchedl1forward("forrecoPhimatchedl1forward","Phi of for. lead reco jet that matched L1 forward jets", 20,-4,4);
TH2F  forrecoEtavsPhimatchedl1forward("forrecoEtavsPhimatchedl1forward", "Eta vs Phi for forward lead reco that matched L1 For. jets",20,-6,6,20,-4,4);

TH1F forrecoPtmatchedl1central("forrecoPtmatchedl1central","Pt of for. lead reco jet that matched L1 central jets", 70,0,250);
TH1F forrecoEtamatchedl1central("forrecoEtamatchedl1central","Eta of for. lead reco jet that matched L1 central jets", 20,-6,6);
TH1F forrecoPhimatchedl1central("forrecoPhimatchedl1central","Phi of for. lead reco jet that matched L1 central jets", 20,-4,4);

TH1F forrecoPtmatchedl1tau("forrecoPtmatchedl1tau","Pt of for. lead reco jet that matched L1 tau jets", 70,0,250);
TH1F forrecoEtamatchedl1tau("forrecoEtamatchedl1tau","Eta of for. lead reco jet that matched L1 tau jets", 20,-6,6);
TH1F forrecoPhimatchedl1tau("forrecoPhimatchedl1tau","Phi of for. lead reco jet that matched L1 tau jets", 20,-4,4);

TH1F forrecoPtNOTmatched("forrecoPtNOTmatched","Pt of for. lead reco jet that  did not match L1  jets", 70,0,250);
TH1F forrecoEtaNOTmatched("forrecoEtaNOTmatched","Eta of for. lead reco jet that  did not match L1 jets", 20,-6,6);
TH1F forrecoPhiNOTmatched("forrecoPhiNOTmatched","Phi of for. lead reco jet that did not match L1 jets", 20,-4,4);
TH2F  forrecoEtavsPhiNOTmatched("forrecoEtavsPhiNOTmatched", "Eta vs Phi for forward lead reco that did not match L1 jets",20,-6,6,20,-4,4);

TH2F h_l1leaf_rcleaf_pt("h_l1leaf_rcleaf_pt","l1 jets pt versus lead reco jet pt distribution",100,0,250,100,0,250);
TH2F h_l1for_vs_reco_pt("h_l1forward_vs_reco_pt","l1 forward jets pt versus lead reco jet pt",100,0,250,100,0,250);
TH2F h_l1cen_vs_reco_pt("h_l1central_vs_reco_pt","l1 central jets pt versus lead reco jet pt",100,0,250,100,0,250);
TH2F h_l1tau_vs_reco_pt("h_l1tau_vs_reco_pt","l1 tau jets pt versus lead reco jet pt",100,0,250,100,0,250);

TH2F h_l1forUnmatched_hltEtaVsPhi ("h_l1forUnmatched_etaVsPhi","",20,-6,6,20,-4,4);

TH2F h_htComparison    ("htComparison","htComparison",100,0,1000,100,0,1000);

TH1F h_eta_unmatched("h_eta_unmatched","eta dist. for unmatched Pt>50 lead reco jet",20,-6,6);
TH1F h_phi_unmatched("h_phi_unmatched","phi dist. for unmatched Pt>50 lead reco jet",20,-4,4);

//-------------------------------------------------------------
// These parameters will be filled by the script
//-------------------------------------------------------------

float  l1_gctHT           [12];

float  l1CenJet_pt        [maxNL1CenJet];
float  l1CenJet_eta       [maxNL1CenJet];
float  l1CenJet_phi       [maxNL1CenJet];
       		     
float  l1TauJet_pt        [maxNL1TauJet];
float  l1TauJet_eta       [maxNL1TauJet];
float  l1TauJet_phi       [maxNL1TauJet];
                      
float  l1ForJet_pt        [maxNL1ForJet];
float  l1ForJet_eta       [maxNL1ForJet];
float  l1ForJet_phi       [maxNL1ForJet];
       
float  hltJet_pt          [maxNHltJet];
float  hltJet_phi         [maxNHltJet];
float  hltJet_eta         [maxNHltJet];
       
float  gctHT              [maxNGctHt ];
float  hltHT;

int    nMatchCenJets      [maxPtBin+1];
int    nMatchTauJets      [maxPtBin+1];
int    nMatchForJets      [maxPtBin+1];
       
int    nTotalCenJets      [maxPtBin+1];
int    nTotalTauJets      [maxPtBin+1];
int    nTotalForJets      [maxPtBin+1];

double ptBinCenter        [maxPtBin+1];
double ptBinError         [maxPtBin+1];
				  
double ptBinMatchEffi_cen [maxPtBin+1];
double ptBinMatchEffi_tau [maxPtBin+1];
double ptBinMatchEffi_for [maxPtBin+1];
				  
double ptBinMatchEEffi_cen[maxPtBin+1];
double ptBinMatchEEffi_tau[maxPtBin+1];
double ptBinMatchEEffi_for[maxPtBin+1];

int   nHLTJetCands;					 
int   nL1GctEtHads;
int   nL1CenJet;
int   nL1TauJet;
int   nL1ForJet;

TChain *tree[maxPtBin+1][maxJetThreshold+1];

//-------------------------------------------------------------
// Geometry functions
//-------------------------------------------------------------

float deltaPhi (float phi1, float phi2) { 
  float ans = fabs(phi1 - phi2);
  while(ans > TMath::Pi())
    ans = fabs(ans - (2.0*(TMath::Pi())));
  return(ans);
}

float deltaR2 (float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1 - eta2;
  float dphi = deltaPhi (phi1, phi2);
  return deta*deta + dphi*dphi;
}

float deltaR (float eta1, float phi1, float eta2, float phi2) {
  return (float) TMath::Sqrt(deltaR2 (eta1, phi1, eta2, phi2));
}

//-------------------------------------------------------------
// Create the TChain
//-------------------------------------------------------------


void createChain(){

  char tempFileName[500];
  
  for (int iPtBin = minPtBin; iPtBin <= maxPtBin; iPtBin++){
    for (int iJetThreshold = minJetThreshold; iJetThreshold <= maxJetThreshold; iJetThreshold++){
      
      tree[iPtBin][iJetThreshold] = new TChain("l1SkimTree","");
      
      for (int iJob = 1; iJob <= nJobsPerThreshold; iJob++){
	
	sprintf(tempFileName,
		"/uscms/home/eberry/data/L1AndHLTOnMC_AllJetThresholds/L1SkimAnalyzerOutput_L1EmulatorOnMC_PtBin%d_%dGeV_job%d_%devents_withHLT.root",
		iPtBin,
		iJetThreshold,
		iJob,
		nEventsPerJob);
	
	tree[iPtBin][iJetThreshold] -> Add(tempFileName);
	
      }
    }
  }
}

//-------------------------------------------------------------
// Do the analysis
//-------------------------------------------------------------

void analyzeChain(){

  int nEvents;
  
  float l1CenJet_tempPt , l1TauJet_tempPt , l1ForJet_tempPt ;
  float l1CenJet_tempEta, l1TauJet_tempEta, l1ForJet_tempEta;
  float l1CenJet_tempPhi, l1TauJet_tempPhi, l1ForJet_tempPhi;
  float hltJet_tempPt   , hltJet_tempEta  , hltJet_tempPhi;
  float gctHT, hltHT;

  float matchDeltaR;

  double weight;

  for (int iPtBin = minPtBin; iPtBin <= maxPtBin; iPtBin++){

    //-------------------------------------------------------------
    // Set branch addresses
    //-------------------------------------------------------------

    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("gctHT",           l1_gctHT);

    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1CenJet_pt"    , l1CenJet_pt    );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1CenJet_eta"   , l1CenJet_eta   );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1CenJet_phi"   , l1CenJet_phi   );

    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1TauJet_pt"    , l1TauJet_pt    );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1TauJet_eta"   , l1TauJet_eta   );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1TauJet_phi"   , l1TauJet_phi   );

    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1ForJet_pt"    , l1ForJet_pt    );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1ForJet_eta"   , l1ForJet_eta   );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("l1ForJet_phi"   , l1ForJet_phi   );

    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("hltJet_pt"      , hltJet_pt      );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("hltJet_eta"     , hltJet_eta     );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("hltJet_phi"     , hltJet_phi     );

    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("nHLTJetCands"   ,&nHLTJetCands   );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("nL1CenJet"      ,&nL1CenJet      );
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("nL1TauJet"      ,&nL1TauJet      ); 
    tree[iPtBin][anaJetThreshold] -> SetBranchAddress("nL1ForJet"      ,&nL1ForJet      );

    //-------------------------------------------------------------
    // Determine how many events to look at 
    //-------------------------------------------------------------

    if (nEventsToExaminePerFile <= 0) nEvents = tree[iPtBin][anaJetThreshold] -> GetEntries();
    if (nEventsToExaminePerFile >  0) nEvents = nEventsToExaminePerFile;

    //-------------------------------------------------------------
    // Set counters to zero  
    //-------------------------------------------------------------

    nMatchCenJets[iPtBin] = 0;
    nMatchTauJets[iPtBin] = 0;
    nMatchForJets[iPtBin] = 0;

    nTotalCenJets[iPtBin] = 0;
    nTotalTauJets[iPtBin] = 0;
    nTotalForJets[iPtBin] = 0;
    
    //-------------------------------------------------------------
    // Let the user know where you are
    //-------------------------------------------------------------

    cout << "This is p_T bin " << iPtBin 
	 << " (" << minPtBin << "," << maxPtBin << "), with " 
	 << nEvents << " events" << endl;
    
    //-------------------------------------------------------------
    // Loop over events
    //-------------------------------------------------------------

    for (int iEvent = 1; iEvent <= nEvents; iEvent++){    

      tree[iPtBin][anaJetThreshold] -> GetEntry(iEvent);

      //-------------------------------------------------------------
      //
      //-------------------------------------------------------------

      gctHT = l1_gctHT[1] * (1.0 / gctLSB);

      hltHT = 0.0;
      
      for (int iHltJet = 0; iHltJet < nHLTJetCands; iHltJet++)
	if (hltJet_pt[iHltJet] > anaJetThreshold) hltHT += hltJet_pt[iHltJet];

      h_htComparison.Fill(gctHT, hltHT);
      
      //-------------------------------------------------------------
      // Print out the event number
      //-------------------------------------------------------------

      if (iEvent%10000 == 0) cout << "  Processing event " << iEvent << endl;
      
      //-------------------------------------------------------------
      // Get vectors for storing deltaR values
      //-------------------------------------------------------------

      std::vector<float>   dr_withfor, dr_withcen, dr_withtau; 
      std::vector<float>   Pt_for, Pt_cen, Pt_tau; 
      std::vector<float>   Eta_for, Eta_cen, Eta_tau; 
      std::vector<float>   Phi_for, Phi_cen, Phi_tau; 
      
      //-------------------------------------------------------------
      // Get lead HLT jet info for jet matching
      //-------------------------------------------------------------    
    
      hltJet_tempPt  = hltJet_pt [0];
      hltJet_tempEta = hltJet_eta[0];
      hltJet_tempPhi = hltJet_phi[0];

      if (hltJet_tempPt < minHltJetPt) continue;

      //-------------------------------------------------------------
      // Central Jet matching
      //-------------------------------------------------------------
            
      for (int iL1CenJet = 0; iL1CenJet < nL1CenJet; iL1CenJet++){
	
	l1CenJet_tempPt     = l1CenJet_pt    [iL1CenJet  ];
	l1CenJet_tempEta    = l1CenJet_eta   [iL1CenJet  ];
	l1CenJet_tempPhi    = l1CenJet_phi   [iL1CenJet  ];

	if ( fabs(l1CenJet_tempEta) > 5.0 || 
	     fabs(l1CenJet_tempPhi) > 3.3   ) continue;

	matchDeltaR = deltaR(hltJet_tempEta, hltJet_tempPhi,
			     l1CenJet_tempEta, l1CenJet_tempPhi );
	
	dr_withcen.push_back( matchDeltaR );
	Pt_cen.push_back (l1CenJet_tempPt);       
	Eta_cen.push_back(l1CenJet_tempEta);       
	Phi_cen.push_back(l1CenJet_tempPhi);       
	nTotalCenJets[iPtBin]++;
      }

      for (int iL1TauJet = 0; iL1TauJet < nL1TauJet; iL1TauJet++){
	
	l1TauJet_tempPt     = l1TauJet_pt    [iL1TauJet  ];
	l1TauJet_tempEta    = l1TauJet_eta   [iL1TauJet  ];
	l1TauJet_tempPhi    = l1TauJet_phi   [iL1TauJet  ];
	
	if ( fabs(l1TauJet_tempEta) > 5.0 || 
	     fabs(l1TauJet_tempPhi) > 3.3   ) continue;
	
	matchDeltaR = deltaR(hltJet_tempEta, hltJet_tempPhi,
			     l1TauJet_tempEta, l1TauJet_tempPhi );

	
	dr_withtau.push_back( matchDeltaR );
	Pt_tau.push_back(l1TauJet_tempPt);     
	Eta_tau.push_back(l1TauJet_tempEta);       
	Phi_tau.push_back(l1TauJet_tempPhi);         	
	nTotalTauJets[iPtBin]++;
      }
      
      for (int iL1ForJet = 0; iL1ForJet < nL1ForJet; iL1ForJet++){
	
	l1ForJet_tempPt     = l1ForJet_pt    [iL1ForJet  ];
	l1ForJet_tempEta    = l1ForJet_eta   [iL1ForJet  ];
	l1ForJet_tempPhi    = l1ForJet_phi   [iL1ForJet  ];
	
	if ( fabs(l1ForJet_tempEta) > 5.0 || 
	     fabs(l1ForJet_tempPhi) > 3.3   ) continue;

	matchDeltaR = deltaR(hltJet_tempEta, hltJet_tempPhi,
			     l1ForJet_tempEta, l1ForJet_tempPhi );

	dr_withfor.push_back( matchDeltaR );
	Pt_for.push_back(l1ForJet_tempPt);       	
	Eta_for.push_back(l1ForJet_tempEta);       
	Phi_for.push_back(l1ForJet_tempPhi);         	
	nTotalForJets[iPtBin]++;
      }

      //-------------------------------------------------------------
      // Which jet has the minimum?
      //-------------------------------------------------------------
      
      int ind_cen = 9999;
      int ind_tau = 9999;
      int ind_for = 9999;

      float mindr_withcenvalue = 9999.0;
      float mindr_withtauvalue = 9999.0;
      float mindr_withforvalue = 9999.0;

      if( dr_withcen.size()){
	std::vector <float>::iterator mindr_withcen = std::min_element( dr_withcen.begin(), dr_withcen.end() );
	ind_cen =  mindr_withcen  -  dr_withcen.begin();
	mindr_withcenvalue = *mindr_withcen;
      }

      if( dr_withtau.size()){
	std::vector <float>::iterator mindr_withtau = std::min_element( dr_withtau.begin(), dr_withtau.end() );
	ind_tau =  mindr_withtau  -  dr_withtau.begin();
	mindr_withtauvalue = *mindr_withtau;
      }
      
      if( dr_withfor.size()){
	std::vector <float>::iterator mindr_withfor = std::min_element( dr_withfor.begin(), dr_withfor.end() );
	ind_for =  mindr_withfor  -  dr_withfor.begin();
	mindr_withforvalue = *mindr_withfor;
      }

      double min_dr= 9999;
      int index_l1jet= 9999 ;
      int matchiswith=0;

      if( mindr_withforvalue < min_dr) {min_dr = mindr_withforvalue; index_l1jet = ind_for; matchiswith = 1;}
      if( mindr_withcenvalue < min_dr) {min_dr = mindr_withcenvalue; index_l1jet = ind_cen; matchiswith = 2;}
      if( mindr_withtauvalue < min_dr) {min_dr = mindr_withtauvalue; index_l1jet = ind_tau; matchiswith = 3;}

      //-------------------------------------------------------------
      // If the min deltaR is less than our cut value, fill the bins
      //-------------------------------------------------------------

      if ( min_dr < matchDeltaR_cut) {
	if(matchiswith == 1) {
	  h_l1leaf_rcleaf_pt.Fill(hltJet_tempPt,Pt_for[index_l1jet]); 
	  h_l1for_vs_reco_pt.Fill(hltJet_tempPt,Pt_for[index_l1jet]);
	  
	  if (fabs(hltJet_tempPt - Pt_for[index_l1jet]) > 10.0){
	    h_l1forUnmatched_hltEtaVsPhi.Fill(hltJet_tempEta,hltJet_tempPhi);
	  }

	}
	else if(matchiswith == 2) {
	  h_l1leaf_rcleaf_pt.Fill(hltJet_tempPt,Pt_cen[index_l1jet]);
	  h_l1cen_vs_reco_pt.Fill(hltJet_tempPt,Pt_cen[index_l1jet]);
	}
	else if(matchiswith == 3) {
	  h_l1leaf_rcleaf_pt.Fill(hltJet_tempPt,Pt_tau[index_l1jet]);
	  h_l1tau_vs_reco_pt.Fill(hltJet_tempPt,Pt_tau[index_l1jet]);
	}
      } 

      //-------------------------------------------------------------
      // Plot the distribution of unmatched hlt jets
      //-------------------------------------------------------------

      else {
	h_l1leaf_rcleaf_pt.Fill(-1,-1);
	if(hltJet_tempPt>50.0){
	  h_eta_unmatched.Fill(hltJet_tempEta);
	  h_phi_unmatched.Fill(hltJet_tempPhi);
	}	
      }

    } // end loop over events

    cout << "There are " << nTotalCenJets[iPtBin] << " L1 CenJets" << endl;
    cout << "There are " << nTotalTauJets[iPtBin] << " L1 TauJets" << endl;
    cout << "There are " << nTotalForJets[iPtBin] << " L1 ForJets" << endl;
    
  } // end loop over pt bins  

  char fileName[500];
  sprintf(fileName,"l1SkimFile.root");

  TFile *plotFile = new TFile(fileName,"RECREATE");
  plotFile -> cd();

  h_htComparison.Write(h_htComparison.GetName());
  h_l1leaf_rcleaf_pt.Write(h_l1leaf_rcleaf_pt.GetName());
  h_l1cen_vs_reco_pt.Write(h_l1cen_vs_reco_pt.GetName());
  h_l1tau_vs_reco_pt.Write(h_l1tau_vs_reco_pt.GetName());
  h_l1for_vs_reco_pt.Write(h_l1for_vs_reco_pt.GetName());
  h_l1forUnmatched_hltEtaVsPhi.Write(h_l1forUnmatched_hltEtaVsPhi.GetName());

  
} //end function

//-------------------------------------------------------------
// Main analysis function
//-------------------------------------------------------------

void matchJets_likeZeynep(){

  createChain();
  analyzeChain();
}
