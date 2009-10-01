#ifndef ANALYZERS_CSCTFMUONPTANALYZER_FILLCSCTFMUONTREE_h
#define ANALYZERS_CSCTFMUONPTANALYZER_FILLCSCTFMUONTREE_h

#include <string>

#include "Analyzers/CSCTFMuonPtAnalyzer/interface/CSCTFMuonTree.h"

class TFile;
class TTree;

class FillCSCTFMuonTree
{

   public:
      FillCSCTFMuonTree();
      virtual ~FillCSCTFMuonTree();

      void init(std::string filename, CSCTFMuonTree* treePtr);
      CSCTFMuonTree* getTreePtr();
      void fill();
      void finalize();

   private:
      FillCSCTFMuonTree(const FillCSCTFMuonTree&); // stop default

      const FillCSCTFMuonTree& operator=(const FillCSCTFMuonTree&); // stop default

      CSCTFMuonTree* m_treePtr;
      TFile* m_file;
      TTree* m_tree;

};


#endif


