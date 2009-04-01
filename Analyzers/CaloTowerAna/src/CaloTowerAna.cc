#include "Analyzers/CaloTowerAna/interface/CaloTowerAna.h"

void CaloTowerAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //-----------------------------------------------------
  // Use the appropriate name spaces
  //-----------------------------------------------------

  using namespace edm;

  //-----------------------------------------------------
  // Initialize the ROOT tree
  //-----------------------------------------------------
  
  m_caloTowerTree.init();

  //-----------------------------------------------------
  // Pass on necessary pointers
  //-----------------------------------------------------

  m_event = &iEvent;
  m_setup = &iSetup;

  //-----------------------------------------------------
  // Get run and event information
  //-----------------------------------------------------
  
  int run   = (int) m_event -> id().run();
  int event = (int) m_event -> id().event();

  m_caloTowerTree.run   = run;
  m_caloTowerTree.event = event;

  //-----------------------------------------------------
  // Run helper analysis functions
  //-----------------------------------------------------

  //analyzeCaloTowers();
  analyzeJets();

  //-----------------------------------------------------
  // Fill ROOT tree
  //-----------------------------------------------------
  
  m_fillTree.fill();

}

template <class Digi>
void CaloTowerAna::analyzeDigi(Digi& digi, int nrjet, int rjet_nct, int rjet_ct_ndigi){
  
  //-----------------------------------------------------
  // Get digi location information
  //-----------------------------------------------------
  
  HcalDetId cell(digi.id());
  
  int h_ieta  = (int) cell.ieta();
  int h_iphi  = (int) cell.iphi();
  int h_depth = (int) cell.depth();

  m_caloTowerTree.rjet_ct_digi_ieta [nrjet][rjet_nct][rjet_ct_ndigi] = h_ieta;
  m_caloTowerTree.rjet_ct_digi_iphi [nrjet][rjet_nct][rjet_ct_ndigi] = h_iphi;
  m_caloTowerTree.rjet_ct_digi_depth[nrjet][rjet_nct][rjet_ct_ndigi] = h_depth;

  //-----------------------------------------------------
  // Get digi conditions for linearization
  //-----------------------------------------------------

  const HcalCalibrations& calibrations = m_conditions -> getHcalCalibrations(cell);
  const HcalQIECoder* channelCoder = m_conditions -> getHcalCoder(cell);
  HcalCoderDb coder (*channelCoder, *m_shape); 

  //-----------------------------------------------------
  // Declare a CaloSamples object (fC time samples),
  // and fill it
  //-----------------------------------------------------
  
  CaloSamples tool;
  coder.adc2fC(digi,tool);

  //-----------------------------------------------------
  // Get digi size information (number of [pre]samples)
  //-----------------------------------------------------

  int presamples = (int) digi.presamples();
  int size       = (int) digi.size();
  
  m_caloTowerTree.rjet_ct_digi_ps   [nrjet][rjet_nct][rjet_ct_ndigi] = presamples;
  m_caloTowerTree.rjet_ct_digi_size [nrjet][rjet_nct][rjet_ct_ndigi] = size;

  //-----------------------------------------------------
  // Loop over the time samples from this CaloSample
  //-----------------------------------------------------
  
  for( int ii=0; ii<tool.size(); ii++ ) { 
    
    int   adc    = (int)   digi[ii].adc();
    int   capid  = (int)   digi[ii].capid();
    float fC     = (float) tool[ii];
    float ped    = (float) calibrations.pedestal(capid);	
    float gain   = (float) calibrations.rawgain(capid);  
    float rcgain = (float) calibrations.respcorrgain(capid);
    
    m_caloTowerTree.rjet_ct_digi_adc   [nrjet][rjet_nct][rjet_ct_ndigi][ii] = adc;
    m_caloTowerTree.rjet_ct_digi_capid [nrjet][rjet_nct][rjet_ct_ndigi][ii] = capid;
    m_caloTowerTree.rjet_ct_digi_fC    [nrjet][rjet_nct][rjet_ct_ndigi][ii] = fC;
    m_caloTowerTree.rjet_ct_digi_ped   [nrjet][rjet_nct][rjet_ct_ndigi][ii] = ped;
    m_caloTowerTree.rjet_ct_digi_gain  [nrjet][rjet_nct][rjet_ct_ndigi][ii] = gain;
    m_caloTowerTree.rjet_ct_digi_rcgain[nrjet][rjet_nct][rjet_ct_ndigi][ii] = rcgain;

  } 
}

//-----------------------------------------------------
// Can we match an HcalDetId to a digi?
//-----------------------------------------------------

bool CaloTowerAna::matchHcalDetId(edm::Handle<HBHEDigiCollection> & hbhedigis,
				  edm::Handle<HODigiCollection>   & hodigis,
				  edm::Handle<HFDigiCollection>   & hfdigis,
				  HcalDetId hcalDetId){
  
  bool matched;

  if (hcalDetId.subdet() == HcalBarrel || hcalDetId.subdet() == HcalEndcap ){
    HBHEDigiCollection::const_iterator hbhedigi = hbhedigis -> find(hcalDetId);
    matched = (hbhedigi != hbhedigis -> end());
  }
  
  else if (hcalDetId.subdet() == HcalOuter ){
    HODigiCollection::const_iterator hodigi = hodigis -> find(hcalDetId);
    matched = (hodigi != hodigis -> end());
  }
  
  else if (hcalDetId.subdet() == HcalForward){
    HFDigiCollection::const_iterator hfdigi = hfdigis -> find(hcalDetId);
    matched = (hfdigi != hfdigis -> end());
  }
  
  else {
    cout << "hcalDetId.subdet() == " << hcalDetId.subdet() << endl;
    matched = false;
  }
  
  return matched;
  

}

void CaloTowerAna::processHcalDetId(edm::Handle<HBHEDigiCollection> & hbhedigis,
				    edm::Handle<HODigiCollection>   & hodigis,
				    edm::Handle<HFDigiCollection>   & hfdigis,
				    HcalDetId hcalDetId,
				    int nrjet, int rjet_nct, int rjet_ct_ndigi){

  if (hcalDetId.subdet() == HcalBarrel || hcalDetId.subdet() == HcalEndcap ){
    HBHEDigiCollection::const_iterator hbhedigi = hbhedigis -> find(hcalDetId);
    analyzeDigi(*hbhedigi,nrjet, rjet_nct, rjet_ct_ndigi);
  }
  
  else if (hcalDetId.subdet() == HcalOuter ){    
    HODigiCollection::const_iterator hodigi = hodigis -> find(hcalDetId);
    analyzeDigi(*hodigi,nrjet, rjet_nct, rjet_ct_ndigi);
  }
  
  else if (hcalDetId.subdet() == HcalForward){
    HFDigiCollection::const_iterator hfdigi = hfdigis -> find(hcalDetId);
    analyzeDigi(*hfdigi,nrjet, rjet_nct, rjet_ct_ndigi);
  }
  
  else {
    cout << "***hcalDetId.subdet() == " << hcalDetId.subdet() << endl;
  }

}

void CaloTowerAna::analyzeJets(){

  //-----------------------------------------------------
  // Appropriate namespaces
  //-----------------------------------------------------

  using namespace edm;
  using namespace reco;

  //-----------------------------------------------------
  // Declare handles
  //-----------------------------------------------------

  edm::Handle< View<Candidate> > jets;
  edm::Handle< HBHEDigiCollection > hbhedigis;
  edm::Handle< HODigiCollection > hodigis;
  edm::Handle< HFDigiCollection > hfdigis;

  //-----------------------------------------------------
  // Make sure your tags are right
  //-----------------------------------------------------

  bool gotJetHandle  = m_event -> getByLabel(m_jetTag     ,jets     );
  bool gotHBHEHandle = m_event -> getByLabel(m_hcalDigiTag,hbhedigis);
  bool gotHOHandle   = m_event -> getByLabel(m_hcalDigiTag,hodigis  );
  bool gotHFHandle   = m_event -> getByLabel(m_hcalDigiTag,hfdigis  );

  if (!gotJetHandle ) LogWarning("CaloTowerAna") << "Could not extract jets with tag: "       << m_jetTag      << endl;
  if (!gotHBHEHandle) LogWarning("CaloTowerAna") << "Could not extract HBHE digis with tag: " << m_hcalDigiTag << endl;
  if (!gotHOHandle  ) LogWarning("CaloTowerAna") << "Could not extract HO   digis with tag: " << m_hcalDigiTag << endl;
  if (!gotHFHandle  ) LogWarning("CaloTowerAna") << "Could not extract HF   digis with tag: " << m_hcalDigiTag << endl;
  
  //-----------------------------------------------------
  // Get setup information for decoding digis
  //-----------------------------------------------------

  m_setup -> get<HcalDbRecord>().get(m_conditions);
  m_shape = m_conditions -> getHcalShape();

  //-----------------------------------------------------
  // Declare useful values
  //-----------------------------------------------------

  // Reco jet information
  int nrjet = 0;
  float rjet_e, rjet_et, rjet_pt;
  float rjet_eta, rjet_phi;

  // CaloTower information
  int rjet_nct;
  int rjet_ct_ieta, rjet_ct_iphi;
  int rjet_ct_ndigi;
  float rjet_ct_hadE, rjet_ct_emE, rjet_ct_outE;
  
  // Digi matching information
  bool matched = false;

  //-----------------------------------------------------
  // Loop over the jets
  //-----------------------------------------------------

  if (gotJetHandle){
  
    for (unsigned int iJet = 0; iJet < jets -> size(); iJet++){
      
      //-----------------------------------------------------
      // Make sure your jet reference is OK 
      //-----------------------------------------------------
      
      CaloJetRef jetRef = ((*jets).refAt(iJet)).castTo<CaloJetRef>();
      
      if (!jetRef.isNull()){
	
	rjet_e   = (float) (*jetRef).energy();
	rjet_et  = (float) (*jetRef).et();
	rjet_pt  = (float) (*jetRef).pt(); 
	rjet_eta = (float) (*jetRef).eta();
	rjet_phi = (float) (*jetRef).phi();
	
	m_caloTowerTree.rjet_e   [nrjet] = rjet_e;
	m_caloTowerTree.rjet_et  [nrjet] = rjet_et;
	m_caloTowerTree.rjet_pt  [nrjet] = rjet_pt;
	m_caloTowerTree.rjet_eta [nrjet] = rjet_eta;
	m_caloTowerTree.rjet_phi [nrjet] = rjet_phi;
	
	//-----------------------------------------------------
	// Loop over the CaloTowers associated with this jet
	//-----------------------------------------------------
	
	rjet_nct = 0;
	
	for (unsigned int iCaloTower = 0; iCaloTower < (*jetRef).getCaloConstituents().size(); iCaloTower++){
	  
	  //-----------------------------------------------------
	  // Make sure this CaloTower pointer is OK
	  //-----------------------------------------------------
	  
	  CaloTowerPtr caloTowerPtr = (*jetRef).getCaloConstituent(iCaloTower);
	  
	  // We will want to count the number of matched, HCAL digis we find in the CaloTower
	  
	  if (!caloTowerPtr.isNull()){
	    
	    //-----------------------------------------------------
	    // Get CaloTower location information
	    //-----------------------------------------------------
	    
	    rjet_ct_ieta = (int) (*caloTowerPtr).id().ieta();
	    rjet_ct_iphi = (int) (*caloTowerPtr).id().iphi(); 

	    rjet_ct_emE  = (float) (*caloTowerPtr).emEnergy();
	    rjet_ct_hadE = (float) (*caloTowerPtr).hadEnergy();
	    rjet_ct_outE = (float) (*caloTowerPtr).outerEnergy();
	    
	    m_caloTowerTree.rjet_ct_ieta [nrjet][rjet_nct] = rjet_ct_ieta;
	    m_caloTowerTree.rjet_ct_iphi [nrjet][rjet_nct] = rjet_ct_iphi;
	    
	    m_caloTowerTree.rjet_ct_emE  [nrjet][rjet_nct] = rjet_ct_emE ;
	    m_caloTowerTree.rjet_ct_hadE [nrjet][rjet_nct] = rjet_ct_hadE;
	    m_caloTowerTree.rjet_ct_outE [nrjet][rjet_nct] = rjet_ct_outE;

	    //-----------------------------------------------------
	    // Loop over the DetId's in this calotower
	    // CaloTower constituent DetId's should have the same
	    // ieta & iphi as the CaloTower
	    // UNLESS: CaloTower |ieta| == 28 or 29
	    //-----------------------------------------------------
	    
	    rjet_ct_ndigi = 0;
	    
	    for (size_t iConstituent = 0; iConstituent < caloTowerPtr -> constituentsSize(); ++iConstituent){
	      
	      DetId detId = caloTowerPtr -> constituent(iConstituent);
	    
	      //-----------------------------------------------------
	      // Make sure this DetId is within the HCAL
	      //-----------------------------------------------------
	      
	      if (detId.det() == DetId::Hcal){
		
		//-----------------------------------------------------
		// Get the HcalDetId and find where it's located
		//-----------------------------------------------------
		
		HcalDetId hcalDetId = HcalDetId(detId);
		
		//-----------------------------------------------------
		// Can we match the HcalDetId to a digi?
		//-----------------------------------------------------	     	      
		
		if (gotHBHEHandle && gotHOHandle && gotHFHandle)		
		  matched = matchHcalDetId(hbhedigis,hodigis,hfdigis,hcalDetId);
		
		if (matched){
		  
		  processHcalDetId(hbhedigis,hodigis,hfdigis,hcalDetId,
				   nrjet,rjet_nct,rjet_ct_ndigi);
		  
		  //-----------------------------------------------------
		  // We only want to count matched HCAL digis
		  //-----------------------------------------------------
		  
		  rjet_ct_ndigi++;
		  
		}
		
	      } // end if (DetID is in the HCAL)
	    } // end loop over det id's in the calo towers
	    
	    //-----------------------------------------------------
	    // We only want to count non-null CaloTowers
	    //-----------------------------------------------------
	    
	    m_caloTowerTree.rjet_ct_ndigi[nrjet][rjet_nct] = rjet_ct_ndigi;
	    rjet_nct++;
	    
	  } // end if (caloTowerPtr.isNull()) 
	} // end loop over calo towers      
      
	m_caloTowerTree.rjet_nct[nrjet] = rjet_nct;
	nrjet++;            
      } // end if(jetRef.isNull())
      else {
	cout << "Jet Ref is NULL!" << endl;
      }
    } // end loop over candidates
  } // end if (gotJetHandle)
    
  m_caloTowerTree.nrjet = nrjet;
}

void CaloTowerAna::analyzeCaloTowers(){

  using namespace reco;
  using namespace std;
  
  //---------------------------------------------
  // Get the 'view' of CaloTowers from Candidate
  //---------------------------------------------

  edm::Handle<edm::View<Candidate> > towers;
  bool gotHandle =  m_event -> getByLabel(m_caloTowerTag,towers);
  if (!gotHandle) LogWarning("CaloTowerAna") << "Could not extract calo towers with tag: " << m_caloTowerTag << endl;

  //---------------------------------------------
  // Declare the iterator over the candidates
  //---------------------------------------------

  int nct = 0;
  int ct_ieta, ct_iphi, ct_nhcon;
  int h_ieta, h_iphi, h_depth;
  
  for (unsigned int iTower = 0; iTower < towers -> size(); iTower++){
    
    CaloTowerRef calotower = ((*towers).refAt(iTower)).castTo<CaloTowerRef>();
    if (!calotower.isNull()){
      
      ct_ieta = calotower->id().ieta  ();
      ct_iphi = calotower->id().iphi  ();

      //---------------------------------------------
      // Loop over the constituent detector id's
      //---------------------------------------------
      
      ct_nhcon = 0;
      
      for (size_t iConstituent = 0; iConstituent < calotower -> constituentsSize(); ++iConstituent){
	
	DetId detId = calotower -> constituent(iConstituent);
	
	if (detId.det() == DetId::Hcal){
	  
	  HcalDetId hcalDetId = HcalDetId(detId);
	  
	  h_ieta  = hcalDetId.ieta();
	  h_iphi  = hcalDetId.iphi();
	  h_depth = hcalDetId.depth();
	  
	  if ((h_ieta != ct_ieta || h_iphi != ct_iphi) && abs(h_ieta) != 28 && abs(h_ieta) != 29){
	    cout << "----------CaloTower-----------" << endl;
	    cout << "h_ieta  = " << h_ieta  << ", h_iphi  = " << h_iphi  << ", h_depth = " << h_depth << endl; ;
	    cout << "ct_ieta = " << ct_ieta << ", ct_iphi = " << ct_iphi << endl;	    
	  }
	  
	  // count number of HCAL constituents
	  ct_nhcon++;
	  
	} // end if (detector ID is in the HCAL)
	
      }// end loop over calo tower constituents
      
      nct++;
      
    } // end if(calotower)
    
  } // end loop over candidates
}

void CaloTowerAna::beginJob(const edm::EventSetup&) {}
void CaloTowerAna::endJob() {

  m_fillTree.finalize();

}

DEFINE_FWK_MODULE(CaloTowerAna);
