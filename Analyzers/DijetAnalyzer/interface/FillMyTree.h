#ifndef Analyzers_DijetAnalyzer_FillMyTree_h
#define Analyzers_DijetAnalyzer_FillMyTree_h

#include <string>

#include "Analyzers/DijetAnalyzer/interface/MyTree.h"

class TFile;
class TTree;

class FillMyTree
{

   public:
      FillMyTree();
      virtual ~FillMyTree();

      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void init(std::string filename, MyTree* treePtr);
      MyTree* getTreePtr();
      void fill();
      void finalize();

   private:
      FillMyTree(const FillMyTree&); // stop default

      const FillMyTree& operator=(const FillMyTree&); // stop default

      // ---------- member data --------------------------------
      MyTree* m_treePtr;
      TFile* m_file;
      TTree* m_tree;

};


#endif
