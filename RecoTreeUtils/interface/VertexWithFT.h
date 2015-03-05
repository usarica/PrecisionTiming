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
    inline const reco::Vertex*              GetRecoVtxRef() {return recoVtxRef_;};
    PFCandidateWithFT*                      GetSeedRef();
    vector<pair<PFCandidateWithFT, float> > GetParticles();
    int                                     GetNPart(float dz_cut=0.1, float smearing=0);
    float                                   ComputeTime(float pt_cut=2, float smearing=0);

    //---setters---
    void          SetSeed(PFCandidateWithFT seed);
    
    //---utils---
    void          AddParticle(PFCandidateWithFT particle, float dz=-1000);
    void          RemoveParticle(PFCandidateWithFT particle);
    float         sumPtSquared(float dz_cut=0.2, float pt_cut=1) const;
    bool          hasSeed() {return hasSeed_;};

private:

    vector<pair<PFCandidateWithFT, float> > particles_;
    bool                                    hasSeed_;
    float                                   time_;
    int                                     n_time_tracks_;
    const reco::Vertex*                     recoVtxRef_;
};
   
#endif 
