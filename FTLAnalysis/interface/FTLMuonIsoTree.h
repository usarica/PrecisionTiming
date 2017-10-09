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
    DATA(int, iEvent)                           \
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
    DATA(int, vtx3DIdx)                         \
    DATA(int, vtx4DIdx)                         \
    DATA(int, nTracksVtx)                       \
    DATA(int, nTracksVtxAssoc)                  \
    DATA(float, nVtxs)                          \
    DATA(float, vtx3DX)                         \
    DATA(float, vtx3DY)                         \
    DATA(float, vtx3DZ)                         \
    DATA(float, vtx3DT)                         \
    DATA(float, vtx3DNdof)                      \
    DATA(float, vtx3DChi2)                      \
    DATA(bool, vtx3DIsFake)                     \
    DATA(float, vtx4DX)                         \
    DATA(float, vtx4DY)                         \
    DATA(float, vtx4DZ)                         \
    DATA(float, vtx4DT)                         \
    DATA(float, vtx4DNdof)                      \
    DATA(float, vtx4DChi2)                      \
    DATA(bool, vtx4DIsFake)                     

#define DATA_CLASS_TABLE                                \
    DATA(std::vector<float>, chIsoDR)                   \
    DATA(std::vector<float>, chIsoZCut)                 \
    DATA(std::vector<float>, chIsoZTCut_3sigma)         \
    DATA(std::vector<float>, chIsoZTCut_4sigma)         \
    DATA(std::vector<float>, chIsoZTCut_5sigma)         \
    DATA(std::vector<float>, chIsoZTCut_7sigma)         \
    DATA(std::vector<float>, chIsoZTCut_10sigma)

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif

    
