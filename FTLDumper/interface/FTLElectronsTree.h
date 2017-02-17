#ifndef FTL_ELECTRONS_TREE
#define FTL_ELECTRONS_TREE

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLElectronsTree

#define DATA_TABLE                              \
    DATA(int, event)                            \
    DATA(int, lumi)                             \
    DATA(int, run)                              \
    DATA(float, mass)

#define DATA_CLASS_TABLE                        \
    DATA(vector<int>, idx)                      \
    DATA(vector<float>, pt)                     \
    DATA(vector<float>, eta)                    \
    DATA(vector<float>, phi)                    \
    DATA(vector<float>, energy)                 \
    DATA(vector<float>, sc_energy)              \
    DATA(vector<float>, r9)                     \
    DATA(vector<float>, sIeIe)                  

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
