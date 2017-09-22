#ifndef _FTL_MUONISO_TREE_MINIAOD_
#define _FTL_MUONISO_TREE_MINIAOD_

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLMuonIsoTreeMINIAOD

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
    DATA(float, dz)                             \
    DATA(float, dxy)                            \
    DATA(bool, isLooseMuon)                     \
    DATA(bool, isMediumMuon)                    \
    DATA(bool, isTightMuon)                     \
    DATA(bool, genMatched)                      \
    DATA(bool, genMatchedPrompt)                \
    DATA(bool, genMatchedJet)                   \
    DATA(float, genPt)                          \
    DATA(float, genEta)                         \
    DATA(float, genPhi)                         \
    DATA(int, vtxIdx)                           \
    DATA(float, nVtxs)                          \
    DATA(float, vtxX)                           \
    DATA(float, vtxY)                           \
    DATA(float, vtxZ)                           \
    DATA(float, vtxT)                           \
    DATA(float, vtxNdof)                        \
    DATA(float, vtxChi2)                        \
    DATA(bool, vtxIsFake)                     

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

    

