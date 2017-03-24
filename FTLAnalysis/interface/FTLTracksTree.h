#ifndef FTL_TRACKS_TREE
#define FTL_TRACKS_TREE

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLTracksTree

#define DATA_TABLE                              \
    DATA(int, event)                            \
    DATA(int, lumi)                             \
    DATA(int, run)                              \
    DATA(int, nTracks)                          \
    DATA(int, ftlTotHits)    

#define DATA_CLASS_TABLE                                \
    DATA(vector<int>, idx)                              \
    DATA(vector<int>, charge)                           \
    DATA(vector<float>, p)                              \
    DATA(vector<float>, pt)                             \
    DATA(vector<float>, eta)                            \
    DATA(vector<float>, phi)                            \
    DATA(vector<float>, out_pt)                         \
    DATA(vector<float>, tk_chi2)                        \
    DATA(vector<float>, trkZAtFTL)                      \
    DATA(vector<float>, trkEtaAtFTL)                    \
    DATA(vector<float>, trkPhiAtFTL)                    \
    DATA(vector<float>, genPt)                          \
    DATA(vector<float>, genEta)                         \
    DATA(vector<float>, genPhi)                         \
    DATA(vector<float>, genEnergy)                      \
    DATA(vector<float>, ftlHitsEnergySum)               \
    DATA(vector<int>, ftlNHits)                         \
    DATA(vector<float>, ftlHits3x3Sum)                  \
    DATA(vector<int>, ftlNHits3x3)                      \
    DATA(vector<int>, ftlSeedIdx)                       \
    DATA(vector<float>, ftlSeedEnergy)                  \
    DATA(vector<float>, ftlSeedTime)                    \
    DATA(vector<int>, ftlClusNHits)                     \
    DATA(vector<float>, ftlClusEnergy)                  \
    DATA(vector<float>, ftlClusTime)                    \
    DATA(vector, ftlHitsTrkIdx, <vector<int> >)         \
    DATA(vector, ftlHitsEnergy, <vector<float> >)       \
    DATA(vector, ftlHitsTime, <vector<float> >)         \
    DATA(vector, ftlHitsEta, <vector<float> >)          \
    DATA(vector, ftlHitsPhi, <vector<float> >)          \
    DATA(vector, ftlHitsZ, <vector<float> >)              

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
