#ifndef Analyzers_HcalDigiAnalyzer_DigiTree_h
#define Analyzers_HcalDigiAnalyzer_DigiTree_h

class DigiTree
{

   public:

      enum { MAXNHITS = 5000 };
      enum { MAXNTP   = 2    };
      enum { NQIE     = 10   };
      enum { NADC     = 128  };

      DigiTree();
      virtual ~DigiTree();

      // Public member data
      int run;
      int event;
      int nint;

      // Information about digi hits
      int   nchn;
      int   nhit;
      int   h_id   [MAXNHITS];
      int   h_psam [MAXNHITS];
      int   h_size [MAXNHITS];
      int   h_depth[MAXNHITS];
      int   h_iphi [MAXNHITS];
      int   h_ieta [MAXNHITS];
      float h_eta  [MAXNHITS];
      float h_phi  [MAXNHITS];
      float h_adc  [MAXNHITS][NQIE];
      float h_fC   [MAXNHITS][NQIE];
      float h_ped  [MAXNHITS][NQIE];
      float h_pedc [MAXNHITS][NQIE];
      float h_gain [MAXNHITS][NQIE];
      int   h_capid[MAXNHITS][NQIE];
      int   h_fiber[MAXNHITS][NQIE];
      int   h_fchan[MAXNHITS][NQIE];

      int   t_spike[MAXNHITS];
      int   t_ntp  [MAXNHITS];
      int   t_found[MAXNHITS][MAXNTP];
      int   t_size [MAXNHITS][MAXNTP];
      int   t_psam [MAXNHITS][MAXNTP];
      int   t_ntpts[MAXNHITS][MAXNTP];
      int   t_cET  [MAXNHITS][MAXNTP][NQIE];

      int nadc;

      void init()
      {
	 run   = -999;
	 event = -999;
	 nint = -999;
	 nchn = -999;
	
	 nhit = 0;
	 for( int i=0; i<MAXNHITS; i++ )
	 {
	   
	   h_psam [i] = -999;
	   h_size [i] = -999;
	   h_id   [i] = -999;
	   h_depth[i] = -999;
	   h_ieta [i] = -999;
	   h_iphi [i] = -999;
	   h_eta  [i] = -999;
	   h_phi  [i] = -999;

	   t_spike[i] = -999;
	   t_ntp  [i] = -999;
	   
	   for (int k = 0; k < MAXNTP; k ++){
	   
	     t_found[i][k]    = -999;
	     t_size [i][k]    = -999;
	     t_psam [i][k]    = -999;
	     t_ntpts[i][k]    = -999;

	     for (int l = 0; l < NQIE; l++)
	       t_cET[i][k][l] = -999;
	     
	   }
	   
	   for( int j=0; j<NQIE; j++ )
	     {
	       h_adc  [i][j] = -999;
	       h_fC   [i][j] = -999;
	       h_ped  [i][j] = -999;
	       h_pedc [i][j] = -999;
	       h_gain [i][j] = -999;
	       h_capid[i][j] = -999;
	       h_fiber[i][j] = -999;
	       h_fchan[i][j] = -999;
	     }
	 }
	 
	 nadc = 0;
	 
      }

 private:
      DigiTree(const DigiTree&); // stop default
      
      const DigiTree& operator=(const DigiTree&); // stop default
      

};


#endif
