#ifndef Analyzers_HcalDigiAnalyzer_FillDigiTree_h
#define Analyzers_HcalDigiAnalyzer_FillDigiTree_h
// -*- C++ -*-
//
// Package:     HcalDigiAnalyzer
// Class  :     FillDigiTree
// 
/**\class FillDigiTree FillDigiTree.h Analyzers/HcalDigiAnalyzer/interface/FillDigiTree.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  
//         Created:  Sun Nov 26 16:21:12 CST 2006
// $Id$
//

// system include files
#include <string>

// user include files
#include "Analyzers/HcalDigiAnalyzer/interface/DigiTree.h"

// forward declarations
class TFile;
class TTree;

class FillDigiTree
{

   public:
      FillDigiTree();
      virtual ~FillDigiTree();

      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void init(std::string filename, DigiTree* treePtr);
      DigiTree* getTreePtr();
      void fill();
      void finalize();

   private:
      FillDigiTree(const FillDigiTree&); // stop default

      const FillDigiTree& operator=(const FillDigiTree&); // stop default

      // ---------- member data --------------------------------
      DigiTree* m_treePtr;
      TFile* m_file;
      TTree* m_tree;

};


#endif
