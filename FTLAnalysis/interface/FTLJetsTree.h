#ifndef FTL_JETS_TREE
#define FTL_JETS_TREE

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLJetsTree

#define DATA_TABLE                              \
    DATA(int, event)                            \
    DATA(int, lumi)                             \
    DATA(int, run)                              \
    DATA(float, mass)

#define DATA_CLASS_TABLE                                \
    DATA(vector<int>, idx)                              \
    DATA(vector<float>, pt)                             \
    DATA(vector<float>, eta)                            \
    DATA(vector<float>, phi)                            \
    DATA(vector<float>, energy)                         \
    DATA(vector<float>, sc_energy)                      \
    DATA(vector<float>, phoEfrac)                             

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
