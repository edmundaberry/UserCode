#ifndef Producers_CaloTowersFromTrigPrimsCreator_FillCaloTowersFromTrigPrimsAnalyzerTree_h
#define Producers_CaloTowersFromTrigPrimsCreator_FillCaloTowersFromTrigPrimsAnalyzerTree_h

#include <string>

#include "Producers/CaloTowersFromTrigPrimsCreator/interface/CaloTowersFromTrigPrimsAnalyzerTree.h"

class TFile;
class TTree;

class FillCaloTowersFromTrigPrimsAnalyzerTree
{

   public:
      FillCaloTowersFromTrigPrimsAnalyzerTree();
      virtual ~FillCaloTowersFromTrigPrimsAnalyzerTree();

      void init(std::string filename, CaloTowersFromTrigPrimsAnalyzerTree* treePtr);
      CaloTowersFromTrigPrimsAnalyzerTree* getTreePtr();
      void fill();
      void finalize();

   private:
      FillCaloTowersFromTrigPrimsAnalyzerTree(const FillCaloTowersFromTrigPrimsAnalyzerTree&); 

      const FillCaloTowersFromTrigPrimsAnalyzerTree& operator=(const FillCaloTowersFromTrigPrimsAnalyzerTree&);

      CaloTowersFromTrigPrimsAnalyzerTree* m_treePtr;
      TFile* m_file;
      TTree* m_tree;

};


#endif
