#ifndef FTL_MUONISO_TREE
#define FTL_MUONISO_TREE

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLMuonIsoTree

#define DATA_TABLE                              \
    DATA(int, event)                            \
    DATA(int, lumi)                             \
    DATA(int, run)                              \
    DATA(float, pt)                             \
    DATA(float, eta)                            \
    DATA(float, phi)                            \
    DATA(float, px)                             \
    DATA(float, py)                             \
    DATA(float, pz)                             \
    DATA(float, muonZ)                          \
    DATA(bool, isLooseMuon)                     \
    DATA(bool, isMediumMuon)                    \
    DATA(bool, isTightMuon)                     \
    DATA(bool, genMatched)                      \
    DATA(bool, genMatchedPrompt)                \
    DATA(bool, genMatchedJet)                   \
    DATA(float, genPt)                          \
    DATA(float, genEta)                         \
    DATA(float, genPhi)                         \
    DATA(int, vtxIndex)                         \
    DATA(int, nTracksVtx)                       \
    DATA(int, nTracksVtxAssoc)                  \
    DATA(float, nVtxs)                          \
    DATA(float, vtxX)                           \
    DATA(float, vtxY)                           \
    DATA(float, vtxZ)                           \
    DATA(float, vtxT)                           \
    DATA(float, vtxNdof)                        \
    DATA(float, vtxChi2)                        \
    DATA(bool, vtxIsFake)

#define DATA_CLASS_TABLE                        \
    DATA(std::vector<float>, chIsoDR)           \
    DATA(std::vector<float>, chIsoZCut)         \
    DATA(std::vector<float>, chIsoZTCut_3sigma) \
    DATA(std::vector<float>, chIsoZTCut_4sigma) \
    DATA(std::vector<float>, chIsoZTCut_5sigma) 

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif

    

