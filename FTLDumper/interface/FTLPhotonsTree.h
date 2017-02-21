#ifndef FTL_PHOTONS_TREE
#define FTL_PHOTONS_TREE

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLPhotonsTree

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
    DATA(vector<float>, sIeIe)                  \
    DATA(vector<float>, convRadius)                  

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
