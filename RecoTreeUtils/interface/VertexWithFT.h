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
    float                                    ComputeTime(int particle_type,
                                                         float pt_cut=2, float smearing=0);
    float                                    ComputeTimeBottomUp(int particle_type,
                                                                 float pt_cut=2,
                                                                 float smearing=0);

    //---setters---
    void          SetGenVtxRef(const SimVertex* genVtx, int id);
    void          SetGenVtxRef(math::XYZTLorentzVector vp, int id);
    void          SetSeed(PFCandidateWithFT* seed);
    
    //---utils---
    float         sumPtSquared(float dz_cut=0.2, float pt_cut=1) const;
    void          AddParticle(PFCandidateWithFT* particle, float dz=-1000);
    void          RemoveParticle(PFCandidateWithFT* particle);
    void          FixVtxRefs();
    bool          hasSeed() {return hasSeed_;};

private:

    vector<pair<PFCandidateWithFT*, float> > particles_;
    bool                                     hasSeed_;
    float                                    time_;
    int                                      n_time_tracks_;
    int                                      genVtxId_;
    const SimVertex*                         genVtxRef_;
    const reco::Vertex*                      recoVtxRef_;    
};
   
#endif 
