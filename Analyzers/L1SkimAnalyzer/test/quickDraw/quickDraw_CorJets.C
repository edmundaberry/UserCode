
TChain *chain;

const float phiMin = -0.27;
const float phiMax = -0.15;
const float etaMin = -0.4;
const float etaMax = 0.0;

TH2F *hltLeadCorJet_etaVsPhi_all = new TH2F("hltLeadCorJet_etaVsPhi_all","",
					 82,-5,5,
					 72,-TMath::Pi(),TMath::Pi());
TH2F *hltLeadCorJet_etaVsPhi_outer = new TH2F("hltLeadCorJet_etaVsPhi_outer","",
					   82,-5,5,
					   72,-TMath::Pi(),TMath::Pi());

TH1F *hltLeadCorJet_phi_all = new TH1F("hltLeadCorJet_phi_all","",
				    72,-TMath::Pi(),TMath::Pi());

TH1F *hltLeadCorJet_et_phiStrip = new TH1F("hltLeadCorJet_et_phiStrip","",
					100,0,10);

TH1F *hltLeadCorJet_et_phiSpike = new TH1F("hltLeadCorJet_et_phiSpike","",
					100,0,10);

void quickDraw_CorJets(){
  
  chain = new TChain("l1SkimTree","");
    
  cout << "Hello world!!!" << endl;

  char fileName[300];
  
  for (int i = 1; i <= 10; i++){
     
     sprintf(fileName,"/uscms/home/eberry/3DayLifetime/223_COR/L1SkimAnalyzerOutput_L1EmulatorOnMC_PtBin1_20GeV_job%d_10000events_223_COR.root",i);
     
     chain -> Add(fileName);
    
  }
  
  
  chain -> Draw("hltCorJet_phi[0]:hltCorJet_eta[0]>>hltLeadCorJet_etaVsPhi_all",""); 
  chain -> Draw("hltCorJet_phi[0]:hltCorJet_eta[0]>>hltLeadCorJet_etaVsPhi_outer","abs(hltCorJet_et[0])>2.0");
  chain -> Draw("hltCorJet_et[0]>>hltLeadCorJet_et_phiStrip","hltCorJet_eta[0] < 0 && hltCorJet_eta[0] > -0.4 && hltCorJet_phi[0] < -0.15 && hltCorJet_phi[0] >  -0.27");
  chain -> Draw("hltCorJet_phi[0]:hltCorJet_eta[0]>>hltLeadCorJet_etaVsPhi_stripZoom","hltCorJet_et[0] > 1.5 && hltCorJet_eta[0] < 0 && hltCorJet_eta[0] > -0.4 && hltCorJet_phi[0] < -0.15 && hltCorJet_phi[0] >  -0.27");


  TFile *outFile = new TFile("hltLeadCorJet_etaVsPhi.root","RECREATE");
  outFile -> cd();
  
  hltLeadCorJet_phi_all        -> Write();
  hltLeadCorJet_etaVsPhi_all   -> Write();
  hltLeadCorJet_etaVsPhi_outer -> Write(); 
  hltLeadCorJet_et_phiStrip    -> Write();

  return;
  

}  

