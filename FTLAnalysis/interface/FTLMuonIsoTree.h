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
    DATA(float, simPVX)                         \
    DATA(float, simPVY)                         \
    DATA(float, simPVZ)                         \
    DATA(float, simPVT)                         \
    DATA(unsigned int, nMuons)                  \
    DATA(unsigned int, nVtx3D)                  \
    DATA(unsigned int, nVtx4D)                  

#define DATA_CLASS_TABLE                        \
    DATA(std::vector<float>, muon_pt)                             \
    DATA(std::vector<float>, muon_eta)                            \
    DATA(std::vector<float>, muon_phi)                            \
    DATA(std::vector<float>, muon_px)                             \
    DATA(std::vector<float>, muon_py)                             \
    DATA(std::vector<float>, muon_pz)                             \
    DATA(std::vector<float>, muon_vx)                          \
    DATA(std::vector<float>, muon_vy)                          \
    DATA(std::vector<float>, muon_vz)                          \
    DATA(std::vector<float>, muon_t)                          \
    DATA(std::vector<float>, muon_terr)                          \
    DATA(std::vector<unsigned int>, isLooseMuon)                     \
    DATA(std::vector<unsigned int>, isMediumMuon)                    \
    DATA(std::vector<unsigned int>, isTightMuon)                     \
    DATA(std::vector<unsigned int>, muonGenMatched)                      \
    DATA(std::vector<unsigned int>, muonGenMatchedPrompt)                \
    DATA(std::vector<unsigned int>, muonGenMatchedJet)                   \
    DATA(std::vector<float>, muonGenPt)                          \
    DATA(std::vector<float>, muonGenEta)                         \
    DATA(std::vector<float>, muonGenPhi)                         \
    DATA(std::vector<float>, muonGenJetE)                           \
    DATA(std::vector<float>, muonGenJetPt)                          \
    DATA(std::vector<float>, muonGenJetEta)                         \
    DATA(std::vector<float>, muonGenJetPhi)                         \
    DATA(std::vector<int>, muonTrkId)                               \
    DATA(std::vector<unsigned int>, muonVtx3DId)                \
    DATA(std::vector<float>, muonIP3DVtx3D)                      \
    DATA(std::vector<float>, muondIP3DVtx3D)                     \
    DATA(std::vector<unsigned int>, muonVtx4DId)                \
    DATA(std::vector<float>, muonIP3DVtx4D)                      \
    DATA(std::vector<float>, muondIP3DVtx4D)                     \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_unassociated)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_associationrank_0)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_associationrank_1)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_nodzcut_unassociated)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_sipcut_unassociated)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_sipcut_associationrank_0)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx3D_sipcut_associationrank_1)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_unassociated)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_associationrank_0)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_associationrank_1)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_associationrank_2)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_nodzcut_unassociated)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_nodzcut_associationrank_0)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_nodzcut_associationrank_1)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_nodzcut_associationrank_2)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_sipcut_unassociated)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_sipcut_associationrank_0)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_sipcut_associationrank_1)          \
    DATA(std::vector<float>, muon_isosumtrackpt_vtx4D_sipcut_associationrank_2)          \
    DATA(std::vector<float>, track_pt)                             \
    DATA(std::vector<float>, track_eta)                            \
    DATA(std::vector<float>, track_phi)                            \
    DATA(std::vector<float>, track_px)                             \
    DATA(std::vector<float>, track_py)                             \
    DATA(std::vector<float>, track_pz)                             \
    DATA(std::vector<float>, track_vx)                          \
    DATA(std::vector<float>, track_vy)                          \
    DATA(std::vector<float>, track_vz)                          \
    DATA(std::vector<float>, track_t)                          \
    DATA(std::vector<float>, track_terr)                          \
    DATA(std::vector<unsigned int>, trackVtx3DId)                \
    DATA(std::vector<float>, trackIP3DVtx3D)                      \
    DATA(std::vector<float>, trackdIP3DVtx3D)                     \
    DATA(std::vector<int>, trackVtx3DAssociationRank)   \
    DATA(std::vector<unsigned int>, trackVtx4DId)                \
    DATA(std::vector<float>, trackIP3DVtx4D)                      \
    DATA(std::vector<float>, trackdIP3DVtx4D)                     \
    DATA(std::vector<int>, trackVtx4DAssociationRank)   \
    DATA(std::vector<float>, vtx3D_vx)                             \
    DATA(std::vector<float>, vtx3D_vy)                             \
    DATA(std::vector<float>, vtx3D_vz)                             \
    DATA(std::vector<unsigned int>, vtx3D_ntrks)                \
    DATA(std::vector<float>, vtx3D_ndof)                \
    DATA(std::vector<float>, vtx3D_chisq)                \
    DATA(std::vector<float>, vtx4D_vx)                             \
    DATA(std::vector<float>, vtx4D_vy)                             \
    DATA(std::vector<float>, vtx4D_vz)                             \
    DATA(std::vector<float>, vtx4D_t)                          \
    DATA(std::vector<float>, vtx4D_terr)                          \
    DATA(std::vector<unsigned int>, vtx4D_ntrks)                \
    DATA(std::vector<float>, vtx4D_ndof)                \
    DATA(std::vector<float>, vtx4D_chisq)                


#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif

    

