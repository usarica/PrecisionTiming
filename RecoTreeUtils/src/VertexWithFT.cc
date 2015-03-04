#include "FastTiming/RecoTreeUtils/interface/VertexWithFT.h"

//----------ctors-------------------------------------------------------------------------

VertexWithFT::VertexWithFT()
{}

VertexWithFT::VertexWithFT(const reco::Vertex* recoVtx):
    reco::Vertex(recoVtx->position(), recoVtx->error(), recoVtx->chi2(), recoVtx->ndof(), recoVtx->tracksSize()),
    time_(-1000), n_time_tracks_(0)
{
    recoVtxRef_ = recoVtx;
}

//----------Get particles container-------------------------------------------------------
vector<pair<PFCandidateWithFT*, int> > VertexWithFT::GetParticles()
{
    return particles_;
}
    
//----------Add particle to the vertex----------------------------------------------------

void VertexWithFT::AddParticle(PFCandidateWithFT* particle, float dz)
{
    float dz_tmp=dz;
    if(dz_tmp == -1000)
        dz_tmp = particle->GetTrack()->dz(this->position());

    particles_.push_back(make_pair(particle, dz_tmp));

    return;
}

//----------Remove particle to the vertex-------------------------------------------------
void VertexWithFT::RemoveParticle(PFCandidateWithFT* particle)
{
}
        
//----------Compute vertex time-----------------------------------------------------------
float VertexWithFT::ComputeWightedTime(float dz_cut, float smearing)
{
    time_ = 0;
    n_time_tracks_ = 0;
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
        if(fabs(particles_.at(iPart).second) < dz_cut &&
           particles_.at(iPart).first->hasTime() &&
           particles_.at(iPart).first->GetTrack())
        {
            cout << iPart << "  " << particles_.at(iPart).first << endl;
            time_ += particles_.at(iPart).first->GetVtxTime(smearing);
            ++n_time_tracks_;
        }
    }
    if(n_time_tracks_ == 0)
        return -100;
    
    return time_ = time_ / n_time_tracks_;
}
