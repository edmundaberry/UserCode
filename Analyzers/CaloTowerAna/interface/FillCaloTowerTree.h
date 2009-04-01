#ifndef Analyzers_CaloTowerAna_FillCaloTowerTree_h
#define Analyzers_CaloTowerAna_FillCaloTowerTree_h

#include <string>

#include "Analyzers/CaloTowerAna/interface/CaloTowerTree.h"

class TFile;
class TTree;

class FillCaloTowerTree
{

   public:
      FillCaloTowerTree();
      virtual ~FillCaloTowerTree();

      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void init(std::string filename, CaloTowerTree* treePtr);
      CaloTowerTree* getTreePtr();
      void fill();
      void finalize();

   private:
      FillCaloTowerTree(const FillCaloTowerTree&); // stop default

      const FillCaloTowerTree& operator=(const FillCaloTowerTree&); // stop default

      // ---------- member data --------------------------------
      CaloTowerTree* m_treePtr;
      TFile* m_file;
      TTree* m_tree;

};


#endif
