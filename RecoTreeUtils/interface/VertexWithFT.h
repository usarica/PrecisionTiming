#ifndef VertexWithFT_H
#define VertexWithFT_H

#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"

using namespace std;

class PFCandidateWithFT;

class VertexWithFT: public reco::Vertex
{
public:
    //---ctors---
    VertexWithFT();
    VertexWithFT(const reco::Vertex* recoVtx);
    //---dtor---
    ~VertexWithFT() {};

    //---relation operators---
    friend bool operator< (const VertexWithFT& v1, const VertexWithFT& v2);
    friend bool operator> (const VertexWithFT& v1, const VertexWithFT& v2);
    
    //---getters---
    inline int                               GetGenVtxId() {return genVtxId_;};
    inline const SimVertex*                  GetGenVtxRef() {return genVtxRef_;};
    inline const reco::Vertex*               GetRecoVtxRef() {return recoVtxRef_;};
    PFCandidateWithFT*                       GetSeedRef();
    vector<pair<PFCandidateWithFT*, float> > GetParticlesWithDZ();
    vector<PFCandidateWithFT*>               GetParticles();
    int                                      GetNPart() {return n_time_tracks_;};
    int                                      GetNPartEE() {return n_time_tracks_EE_;};
    int                                      GetNPartEB() {return n_time_tracks_EB_;};
    float                                    ComputeTime(int particle_type,
                                                         float pt_cut=2,
                                                         float pz2_cut=20,
                                                         float smearing=0);
    float                                    ComputeTimeBottomUp(int particle_type,
                                                                 float pt_cut=2,
                                                                 float smearing=0);

    //---setters---
    inline void   SetGhostTime(float time){time_ = time; isGhost_ = true; hasSeed_ = false;};
    void          SetGenVtxRef(const SimVertex* genVtx, int id);
    void          SetGenVtxRef(math::XYZTLorentzVector vp, int id);
    void          SetSeed(PFCandidateWithFT* seed);
    
    //---utils---
    inline bool   isGhost() {return isGhost_;};
    inline bool   hasSeed() {return hasSeed_;};
    inline float  GetTimeEE() {return time_EE_;};
    inline float  GetTimeEB() {return time_EB_;};
    float         sumPtSquared(float dz_cut=0.2, float pt_cut=2, int EB=-1) const;
    void          AddParticle(PFCandidateWithFT* particle, float dz=-1000);
    void          RemoveParticle(PFCandidateWithFT* particle);
    void          FixVtxRefs();

private:

    vector<pair<PFCandidateWithFT*, float> > particles_;
    bool                                     hasSeed_;
    bool                                     isGhost_;
    float                                    time_;
    float                                    time_EB_;
    float                                    time_EE_;
    int                                      n_time_tracks_;
    int                                      n_time_tracks_EE_;
    int                                      n_time_tracks_EB_;
    int                                      genVtxId_;
    const SimVertex*                         genVtxRef_;
    const reco::Vertex*                      recoVtxRef_;    
};
   
#endif 
