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
    vector<pair<PFCandidateWithFT, float> > GetParticles();
    int                                     GetNPart(float dz_cut=0.1, float smearing=0);
    float                                   ComputeWightedTime(float dz_cut=0.1, float smearing=0);
    
    //---utils---
    void          AddParticle(PFCandidateWithFT particle, float dz=-1000);
    void          RemoveParticle(PFCandidateWithFT particle);
    float         sumPtSquared(float dz_cut=0.2, float pt_cut=1) const;

private:

    vector<pair<PFCandidateWithFT, float> > particles_;    
    float                                   time_;
    int                                     n_time_tracks_;
    const reco::Vertex*                     recoVtxRef_;
};
   
#endif 
