#include "FastTiming/RecoTreeUtils/interface/VertexWithFT.h"

//----------ctors-------------------------------------------------------------------------

VertexWithFT::VertexWithFT()
{}

VertexWithFT::VertexWithFT(const reco::Vertex* recoVtx):
    reco::Vertex(recoVtx->position(), recoVtx->error(), recoVtx->chi2(), recoVtx->ndof(), recoVtx->tracksSize()),
    time_(-1000), n_time_tracks_(-1)
{
    recoVtxRef_ = recoVtx;
}

//----------Get particles container-------------------------------------------------------
vector<pair<PFCandidateWithFT, float> > VertexWithFT::GetParticles()
{
    return particles_;
}
    
//----------Add particle to the vertex----------------------------------------------------

void VertexWithFT::AddParticle(PFCandidateWithFT particle, float dz)
{
    float dz_tmp=dz;
    if(dz_tmp == -1000)
        dz_tmp = particle.GetTrack()->dz(this->position());

    particles_.push_back(make_pair(particle, dz_tmp));

    return;
}

//----------Remove particle to the vertex-------------------------------------------------
void VertexWithFT::RemoveParticle(PFCandidateWithFT particle)
{
}

//----------Get the number of particles used for the time reconstruction------------------
int VertexWithFT::GetNPart(float dz_cut, float smearing)
{
    if(n_time_tracks_ == -1)
        ComputeWightedTime(dz_cut, smearing);

    return n_time_tracks_;
}

//----------Compute vertex time-----------------------------------------------------------
float VertexWithFT::ComputeWightedTime(float dz_cut, float smearing)
{
    time_ = 0;
    n_time_tracks_ = 0;
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
        float dz_tmp = particles_.at(iPart).first.GetPFCandidate()->pt();
        if(dz_tmp < dz_cut && particles_.at(iPart).first.hasTime())
        {
            time_ += particles_.at(iPart).first.GetVtxTime(smearing);
            ++n_time_tracks_;
        }
    }
    if(n_time_tracks_ == 0)
        return -100;
    
    return time_ = time_ / n_time_tracks_;
}

//----------compute sumpt2 using all the particles related to the vtx---------------------

float VertexWithFT::sumPtSquared(float dz_cut, float pt_cut) const
{
    double sum = 0.;
    double pT;
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
        pT = particles_.at(iPart).first.GetPFCandidate()->pt();
        if(fabs(particles_.at(iPart).second) < dz_cut && pT > pt_cut)           
            sum += pT*pT;
    }
    return sum;
}
    
//----------relation operators------------------------------------------------------------

bool operator< (const VertexWithFT& v1, const VertexWithFT& v2)
{
    if(v1.sumPtSquared() < v2.sumPtSquared())
        return true;

    return false;
}

bool operator> (const VertexWithFT& v1, const VertexWithFT& v2)
{
    return v2 < v1;
}
