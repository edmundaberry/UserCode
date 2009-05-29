#ifndef Analyzers_CaloTowerComp_FillCaloTowerCompTree_h
#define Analyzers_CaloTowerComp_FillCaloTowerCompTree_h

#include <string>

#include "Analyzers/CaloTowerComp/interface/CaloTowerCompTree.h"

class TFile;
class TTree;

class FillCaloTowerCompTree
{

   public:
      FillCaloTowerCompTree();
      virtual ~FillCaloTowerCompTree();

      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void init(std::string filename, CaloTowerCompTree* treePtr);
      CaloTowerCompTree* getTreePtr();
      void fill();
      void finalize();

   private:
      FillCaloTowerCompTree(const FillCaloTowerCompTree&); // stop default

      const FillCaloTowerCompTree& operator=(const FillCaloTowerCompTree&); // stop default

      // ---------- member data --------------------------------
      CaloTowerCompTree* m_treePtr;
      TFile* m_file;
      TTree* m_tree;

};


#endif
