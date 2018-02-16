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

#define DATA_CLASS_TABLE                                \
    DATA(vector<int>, idx)                              \
    DATA(vector<float>, pt)                             \
    DATA(vector<float>, eta)                            \
    DATA(vector<float>, phi)                            \
    DATA(vector<float>, energy)                         \
    DATA(vector<float>, sc_energy)                      \
    DATA(vector<float>, sc_eop)                         \
    DATA(vector<float>, r9)                             \
    DATA(vector<float>, sIeIe)                          \
    DATA(vector<float>, nBrem)                          \
    DATA(vector<float>, fBrem)                          \
    DATA(vector<float>, firstBremRadius)                \
    DATA(vector<float>, genBrem)                        \
    DATA(vector<float>, genNBrem)                       \
    DATA(vector<float>, genPt)                          \
    DATA(vector<float>, genEta)                         \
    DATA(vector<float>, genPhi)                         \
    DATA(vector<float>, genEnergy)                      \
    DATA(vector<float>, ftlHitsEnergySum)               \
    DATA(vector<int>, ftlNHits)                         \
    DATA(vector<float>, ftlHits3x3Sum)                  \
    DATA(vector<int>, ftlNHits3x3)                      \
    DATA(vector<float>, ftlSieie)                       \
    DATA(vector<float>, ftlSipip)                       \
    DATA(vector, ftlHitsEleIdx, <vector<int> >)         \
    DATA(vector, ftlHitsEnergy, <vector<float> >)       \
    DATA(vector, ftlHitsTime, <vector<float> >)         \
    DATA(vector, ftlHitsEta, <vector<float> >)          \
    DATA(vector, ftlHitsPhi, <vector<float> >)          

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
