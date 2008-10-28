#ifndef Analyzers_L1SkimAnalyzer_FillL1SkimTree_h
#define Analyzers_L1SkimAnalyzer_FillL1SkimTree_h

#include <string>
#include "Analyzers/L1SkimAnalyzer/interface/L1SkimTree.h"

class TFile;
class TTree;

class FillL1SkimTree
{

 public:

  FillL1SkimTree();
  virtual ~FillL1SkimTree();

  void init (std::string filename, L1SkimTree* treePtr);
  L1SkimTree* getTreePtr();
  void fill();
  void finalize();

 private:

  FillL1SkimTree(const FillL1SkimTree&);
  
  const FillL1SkimTree& operator=(const FillL1SkimTree&);

  L1SkimTree *m_treePtr;
  TFile *m_file;
  TTree *m_tree;

};

#endif
