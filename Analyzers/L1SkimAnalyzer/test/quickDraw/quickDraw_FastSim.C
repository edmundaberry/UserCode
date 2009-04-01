
TChain *chain;

const float phiMin = -0.27;
const float phiMax = -0.15;
const float etaMin = -0.4;
const float etaMax = 0.0;

TH2F *recoLeadJet_etaVsPhi_all = new TH2F("recoLeadJet_etaVsPhi_all","",
					 82,-5,5,
					 72,-TMath::Pi(),TMath::Pi());
TH2F *recoLeadJet_etaVsPhi_outer = new TH2F("recoLeadJet_etaVsPhi_outer","",
					   82,-5,5,
					   72,-TMath::Pi(),TMath::Pi());

TH1F *recoLeadJet_phi_all = new TH1F("recoLeadJet_phi_all","",
				    72,-TMath::Pi(),TMath::Pi());

TH1F *recoLeadJet_et_phiStrip = new TH1F("recoLeadJet_et_phiStrip","",
					100,0,10);

TH1F *recoLeadJet_et_phiSpike = new TH1F("recoLeadJet_et_phiSpike","",
					100,0,10);

void quickDraw_FastSim(){
  
  chain = new TChain("l1SkimTree","");
    
  cout << "Hello world!!!" << endl;

  char fileName[300];
  
  for (int i = 1; i <= 10; i++){
     
     sprintf(fileName,"/uscms/home/eberry/data/FastSim/L1SkimAnalyzerOutput_FastSim_PtBin1_job%d_10000events_fastSim.root",i);
     
     chain -> Add(fileName);
    
  }
  chain -> Draw("recoJet_phi[0]:recoJet_eta[0]>>recoLeadJet_etaVsPhi_all",""); 
  chain -> Draw("recoJet_phi[0]:recoJet_eta[0]>>recoLeadJet_etaVsPhi_outer","abs(recoJet_et[0])>2.0");
  chain -> Draw("recoJet_et[0]>>recoLeadJet_et_phiStrip","recoJet_eta[0] < 0 && recoJet_eta[0] > -0.4 && recoJet_phi[0] < -0.15 && recoJet_phi[0] >  -0.27");
  chain -> Draw("recoJet_phi[0]:recoJet_eta[0]>>recoLeadJet_etaVsPhi_stripZoom","recoJet_et[0] > 1.5 && recoJet_eta[0] < 0 && recoJet_eta[0] > -0.4 && recoJet_phi[0] < -0.15 && recoJet_phi[0] >  -0.27");

  TFile *outFile = new TFile("recoLeadJet_etaVsPhi.root","RECREATE");
  outFile -> cd();
  
  recoLeadJet_phi_all        -> Write();
  recoLeadJet_etaVsPhi_all   -> Write();
  recoLeadJet_etaVsPhi_outer -> Write(); 
  recoLeadJet_et_phiStrip    -> Write();


}  

