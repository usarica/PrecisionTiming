#ifndef _FTL_MUONISO_TREE_
#define _FTL_MUONISO_TREE_

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLMuonIsoTree

#define DATA_TABLE                                              \
    DATA(int, event)                                            \
    DATA(int, lumi)                                             \
    DATA(int, run)                                              \
    DATA(int, iEvent)                                           \
    DATA(int, iMuon)                                            \
    DATA(float, pt)                                             \
    DATA(float, eta)                                            \
    DATA(float, phi)                                            \
    DATA(float, px)                                             \
    DATA(float, py)                                             \
    DATA(float, pz)                                             \
    DATA(float, dz)                                             \
    DATA(float, dxy)                                            \
    DATA(bool, isLooseMuon)                                     \
    DATA(bool, isMediumMuon)                                    \
    DATA(bool, isTightMuon)                                     \
    DATA(bool, genMatched)                                      \
    DATA(bool, genMatchedPrompt)                                \
    DATA(bool, genMatchedJet)                                   \
    DATA(float, genPt)                                          \
    DATA(float, genEta)                                         \
    DATA(float, genPhi)                                         \
    DATA(float, genIso)                                         \
    DATA(int, vtxIdx)                                           \
    DATA(float, nVtxs)                                          \
    DATA(float, vtxX)                                           \
    DATA(float, vtxY)                                           \
    DATA(float, vtxZ)                                           \
    DATA(float, vtxT)                                           \
    DATA(float, genVtxZ)                                        \
    DATA(float, genVtxT)                                        \
    DATA(float, vtxNdof)                                        \
    DATA(float, vtxChi2)                                        \
    DATA(int, nCloseVtxs)                                       \
    DATA(bool, vtxIsFake)                       
    

#define DATA_CLASS_TABLE                                \
    DATA(std::vector<float>, chIsoDR)                   \
    DATA(std::vector<float>, chIsoZCut)                 \
    DATA(std::vector<float>, chIsoZ1Cut)                \
    DATA(std::vector<float>, chIsoZ3Cut)                \
    DATA(std::vector<float>, chIsoZ5Cut)                \
    DATA(std::vector<float>, chIsoZ7Cut)                \
    DATA(std::vector<float>, chIsoZ10Cut)               \
    DATA(std::vector<float>, chIsoZTCut_3sigma)         \
    DATA(std::vector<float>, chIsoZTCut_4sigma)         \
    DATA(std::vector<float>, chIsoZTCut_5sigma)         \
    DATA(std::vector<float>, chIsoZTCut_7sigma)         \
    DATA(std::vector<float>, chIsoZTCut_10sigma)        \
    DATA(std::vector<float>, tracksInConePt)            \
    DATA(std::vector<float>, tracksRemovedPt)           \
    DATA(std::vector<float>, tracksPt)                  \
    DATA(std::vector<float>, tracksEta)                 \
    DATA(std::vector<float>, tracksPhi)                 \
    DATA(std::vector<float>, tracksDZ)                  \
    DATA(std::vector<float>, tracksDXY)                 \
    DATA(std::vector<float>, tracksT)                   \
    DATA(std::vector<bool>, tracksKeepZ)                \
    DATA(std::vector<bool>, tracksKeepT)               
    
#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif

    

