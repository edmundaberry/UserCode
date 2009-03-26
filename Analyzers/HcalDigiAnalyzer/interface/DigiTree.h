#ifndef Analyzers_HcalDigiAnalyzer_DigiTree_h
#define Analyzers_HcalDigiAnalyzer_DigiTree_h

class DigiTree
{

   public:

      enum { MAXNHITS  = 5000 };
      enum { MAXNTP    = 2    };
      enum { NQIE      = 10   };
      enum { MAX_IETA  = 83   };
      enum { MAX_IPHI  = 73   };
      enum { MAX_DEPTH = 4    };
      enum { MAX_NCT   = 5000 };

      DigiTree();
      virtual ~DigiTree();

      // Public member data
      int run;
      int event;

      // Information about digi hits
      int   nchn;
      int   nhit;
      int   nct;

      int   h_id    [MAXNHITS];
      int   h_depth [MAXNHITS];
      int   h_iphi  [MAXNHITS];
      int   h_ieta  [MAXNHITS];
		    
      int   h_psam  [MAX_IETA][MAX_IPHI][MAX_DEPTH];
      int   h_size  [MAX_IETA][MAX_IPHI][MAX_DEPTH];
      float h_rh_GeV_amp [MAX_IETA][MAX_IPHI][MAX_DEPTH];
      float h_rh_fC_amp  [MAX_IETA][MAX_IPHI][MAX_DEPTH];
      float h_correction [MAX_IETA][MAX_IPHI][MAX_DEPTH];
      float h_threshold  [MAX_IETA][MAX_IPHI][MAX_DEPTH];
      int   h_adc   [MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];
      float h_fC    [MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];
      float h_ped   [MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];
      float h_pedc  [MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];
      float h_gain  [MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];
      float h_rcgain[MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];
      int   h_capid [MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];
      int   h_fiber [MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];
      int   h_fchan [MAX_IETA][MAX_IPHI][MAX_DEPTH][NQIE];

      int   t_spike [MAX_IETA][MAX_IPHI][MAX_DEPTH];
      int   t_ntp   [MAX_IETA][MAX_IPHI][MAX_DEPTH];      
      int   t_subdet[MAX_IETA][MAX_IPHI][MAX_DEPTH][MAXNTP];
      int   t_ieta  [MAX_IETA][MAX_IPHI][MAX_DEPTH][MAXNTP];
      int   t_iphi  [MAX_IETA][MAX_IPHI][MAX_DEPTH][MAXNTP];
      int   t_found [MAX_IETA][MAX_IPHI][MAX_DEPTH][MAXNTP];
      int   t_size  [MAX_IETA][MAX_IPHI][MAX_DEPTH][MAXNTP];
      int   t_psam  [MAX_IETA][MAX_IPHI][MAX_DEPTH][MAXNTP];
      int   t_ntpts [MAX_IETA][MAX_IPHI][MAX_DEPTH][MAXNTP];
      int   t_cET   [MAX_IETA][MAX_IPHI][MAX_DEPTH][MAXNTP][NQIE];
      
      int   ct_ieta [MAX_NCT];
      int   ct_iphi [MAX_NCT];

      void init(){

	run   = -999;
	event = -999;
	nchn  = -999;
	nhit  = -999;
	nct   = -999;

	for (int i = 0; i < MAX_NCT; i++){
	  ct_ieta[i] = -999;
	  ct_iphi[i] = -999;
	}
	
	for (int i = 0; i < MAXNHITS; i++){

	  h_id   [i] = -999;
	  h_depth[i] = -999;
	  h_iphi [i] = -999;
	  h_ieta [i] = -999;

	}

	for (int ieta = 0; ieta < MAX_IETA; ieta++){
	  for (int iphi = 0; iphi < MAX_IPHI; iphi++){
	    for (int depth = 0; depth < MAX_DEPTH; depth++){
	      
	      h_threshold  [ieta][iphi][depth] = -999.0;
	      h_correction [ieta][iphi][depth] = -999.0;
	      h_rh_GeV_amp [ieta][iphi][depth] = -999.0;
	      h_rh_fC_amp  [ieta][iphi][depth] = -999.0;
	      
	      h_psam    [ieta][iphi][depth] = -999;
	      h_size    [ieta][iphi][depth] = -999;
	      t_spike   [ieta][iphi][depth] = -999;
	      t_ntp     [ieta][iphi][depth] = -999;

	      for (int iqie = 0; iqie < NQIE; iqie++){

		h_fC    [ieta][iphi][depth][iqie] = -999.0;
		h_ped   [ieta][iphi][depth][iqie] = -999.0;
		h_pedc  [ieta][iphi][depth][iqie] = -999.0;
		h_gain  [ieta][iphi][depth][iqie] = -999.0;
		
		h_adc   [ieta][iphi][depth][iqie] = -999;
		h_capid [ieta][iphi][depth][iqie] = -999;
		h_fiber [ieta][iphi][depth][iqie] = -999;
		h_fchan [ieta][iphi][depth][iqie] = -999;		 		  

	      }

	      for (int itp = 0; itp < MAXNTP; itp++){

		t_subdet[ieta][iphi][depth][itp]  = -999;
		t_ieta  [ieta][iphi][depth][itp]  = -999;
		t_iphi  [ieta][iphi][depth][itp]  = -999;
		t_found [ieta][iphi][depth][itp]  = -999;
		t_size  [ieta][iphi][depth][itp]  = -999;
		t_psam  [ieta][iphi][depth][itp]  = -999;
		t_ntpts [ieta][iphi][depth][itp]  = -999;

		for (int iqie2 = 0; iqie2 < NQIE; iqie2++)
		  t_cET [ieta][iphi][depth][itp][iqie2] = -999;

		
	      }
	    }
	  }
	}
      }

 private:
      DigiTree(const DigiTree&); // stop default
      
      const DigiTree& operator=(const DigiTree&); // stop default
      

};


#endif
