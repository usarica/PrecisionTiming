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

    //---getters---
    inline const reco::Vertex*             GetRecoVtxRef() {return recoVtxRef_;};
    vector<pair<PFCandidateWithFT*, int> > GetParticles();
    inline int                             GetNPart() {return n_time_tracks_;};
    float                                  ComputeWightedTime(float dz=0.1, float smearing=0);
    
    //---utils---
    void          AddParticle(PFCandidateWithFT* particle, float dz=-1000);
    void          RemoveParticle(PFCandidateWithFT* particle);

private:

    vector<pair<PFCandidateWithFT*, int> > particles_;    
    float                                  time_;
    int                                    n_time_tracks_;
    const reco::Vertex*                    recoVtxRef_;
};
   
#endif 
