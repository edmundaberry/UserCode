#ifndef Analyzers_L1SkimAnalyzer_L1RegionLookup_h
#define Analyzers_L1SkimAnalyzer_L1RegionLookup_h

class L1RegionLookup
{

 public:
  
  L1RegionLookup();
  virtual ~L1RegionLookup();

  void SetArrays();
  int getMinHcalIeta(unsigned int etaIndex, bool isForward);
  int getMaxHcalIeta(unsigned int etaIndex, bool isForward);

  //-----------------------------------------------
  // A map from etaIndex to HCAL ieta range
  //-----------------------------------------------
    
 private:
  
  int m_minHcalIetaAbs_for[12];
  int m_maxHcalIetaAbs_for[12];
      
  int m_minHcalIetaAbs_cen[15];
  int m_maxHcalIetaAbs_cen[15];

  L1RegionLookup (const L1RegionLookup&);
  const L1RegionLookup& operator=(const L1RegionLookup&);


};

#endif
