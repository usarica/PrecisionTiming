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
    DATA(float, mass)                           \
    DATA(int, ftlTotHits)

#define DATA_CLASS_TABLE                                \
    DATA(vector<int>, idx)                              \
    DATA(vector<float>, pt)                             \
    DATA(vector<float>, eta)                            \
    DATA(vector<float>, phi)                            \
    DATA(vector<float>, energy)                         \
    DATA(vector<float>, sc_energy)                      \
    DATA(vector<float>, r9)                             \
    DATA(vector<float>, sIeIe)                          \
    DATA(vector<float>, convRadius)                     \
    DATA(vector<float>, convPhi)                        \
    DATA(vector<float>, convZ)                          \
    DATA(vector<float>, genPt)                          \
    DATA(vector<float>, genEta)                         \
    DATA(vector<float>, genPhi)                         \
    DATA(vector<float>, genEnergy)                      \
    DATA(vector<float>, ftlHitsEnergySum)               \
    DATA(vector<int>, ftlNHits)                         \
    DATA(vector<float>, ftlHits3x3Sum)                  \
    DATA(vector<int>, ftlNHits3x3)                      \
    DATA(vector, ftlHitsPhoIdx, <vector<int> >)         \
    DATA(vector, ftlHitsEnergy, <vector<float> >)       \
    DATA(vector, ftlHitsTime, <vector<float> >)         \
    DATA(vector, ftlHitsEta, <vector<float> >)          \
    DATA(vector, ftlHitsPhi, <vector<float> >)          \
    DATA(vector, ftlHitsZ, <vector<float> >)              

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
