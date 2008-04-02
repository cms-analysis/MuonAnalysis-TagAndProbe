#ifndef TagAndProbe_FillEvtTree_h
#define TagAndProbe_FillEvtTree_h
// -*- C++ -*-
//
// Package:     HLTPhysicsAnalyzer
// Class  :     FillEvtTree
// 
/**\class FillEvtTree FillEvtTree.h HLTPhysicsPerformance/HLTPhysicsAnalyzer/interface/FillEvtTree.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  
//         Created:  Tue Oct 31 10:23:41 CST 2006
// $Id$
//

// system include files
#include <string>

// user include files
#include "MuonAnalysis/TagAndProbe/interface/EvtTree.h"

// forward declarations
class TFile;
class TTree;

class FillEvtTree
{

   public:
      FillEvtTree();
      virtual ~FillEvtTree();

      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void init(std::string filename, EvtTree* treePtr);
      EvtTree* getTreePtr();
      void fill();
      void finalize();

   private:
      FillEvtTree(const FillEvtTree&); // stop default

      const FillEvtTree& operator=(const FillEvtTree&); // stop default

      // ---------- member data --------------------------------
      EvtTree* m_treePtr;
      TFile* m_file;
      TTree* m_tree;

};


#endif
