//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 17 22:57:15 2023 by ROOT version 6.14/04
// from TTree gem_hits/GEM TTree with single track hit info
// found on file: trd_singleTrackHits_Run_003200.root
//////////////////////////////////////////////////////////


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"


// Declaration of leaf types
Int_t           event_num;
float           ecal_energy;
float           presh_energy;
float           mult_energy;
//Int_t           parID;
vector<bool>     *parID;
Int_t           gem_nhit;
vector<int>     *xpos;
vector<int>     *ypos;
vector<float>   *zpos;
vector<float>   *dedx;
vector<float>   *zHist;

//Int_t           clu_nhit;
vector<float>   *xposc;
vector<float>   *yposc;
vector<float>   *zposc;
vector<float>   *dedxc;
vector<float>   *widthc;

// List of branches
TBranch        *b_event_num;   //!
TBranch        *b_ecal_energy;   //!
TBranch        *b_presh_energy;   //!
TBranch        *b_mult_energy;   //!
TBranch        *b_gem_nhit;   //!
TBranch        *b_xpos;   //!
TBranch        *b_ypos;   //!
TBranch        *b_zpos;   //!
TBranch        *b_dedx;   //!
TBranch        *b_parID;   //!
TBranch        *b_zHist;   //!
//
//TBranch        *b_clu_nhit;   //!
TBranch        *b_xposc;   //!
TBranch        *b_yposc;   //!
TBranch        *b_zposc;   //!
TBranch        *b_dedxc;   //!
TBranch        *b_widthc;   //!
