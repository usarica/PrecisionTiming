#include "FastTiming/RecoTreeUtils/interface/VertexWithFT.h"

//----------ctors-------------------------------------------------------------------------

VertexWithFT::VertexWithFT()
{}

VertexWithFT::VertexWithFT(const reco::Vertex* recoVtx):
    reco::Vertex(recoVtx->position(), recoVtx->error(), recoVtx->chi2(), recoVtx->ndof(), recoVtx->tracksSize()),
    hasSeed_(false), time_(-1000), n_time_tracks_(-1)
{
    recoVtxRef_ = recoVtx;
}

//----------Get particles container with dz info------------------------------------------
vector<pair<PFCandidateWithFT*, float> > VertexWithFT::GetParticlesWithDZ()
{
    return particles_;
}

//----------Get particles list------------------------------------------------------------
vector<PFCandidateWithFT*> VertexWithFT::GetParticles()
{
    vector<PFCandidateWithFT*> particles;
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
        particles.push_back(particles_.at(iPart).first);
    
    return particles;
}


//----------Set the seed for the combined space-time vtx reco-----------------------------
void VertexWithFT::SetSeed(PFCandidateWithFT* seed)
{
    if(!seed->GetTrack())
        return;
    if(hasSeed_)
        particles_.clear();

    hasSeed_=true;
    AddParticle(seed);
    return;
}
    
//----------Add particle to the vertex----------------------------------------------------
void VertexWithFT::AddParticle(PFCandidateWithFT* particle, float dz)
{
    if(!hasSeed_)
        SetSeed(particle);
    else
    {
        float dz_tmp=dz;
        if(dz_tmp == -1000)
            dz_tmp = particle->GetTrack()->dz(this->position());

        particles_.push_back(make_pair(particle, dz_tmp));
    }
    return;
}

//----------Remove particle to the vertex-------------------------------------------------
void VertexWithFT::RemoveParticle(PFCandidateWithFT* particle)
{
}

//----------Get the number of particles used for the time reconstruction------------------
int VertexWithFT::GetNPart(float dz_cut, float smearing)
{
    if(n_time_tracks_ == -1)
        ComputeTime(dz_cut, smearing);

    return n_time_tracks_;
}

//----------Get a reference to the seed particle------------------------------------------
PFCandidateWithFT* VertexWithFT::GetSeedRef()
{
    if(hasSeed())
        return particles_.at(0).first;

    return NULL;
}

//----------Compute vertex time-----------------------------------------------------------
float VertexWithFT::ComputeTime(float pt_cut, float smearing)
{
    time_ = 0;
    n_time_tracks_ = 0;
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
        float pt_tmp = particles_.at(iPart).first->pt();
        if(pt_tmp > pt_cut && particles_.at(iPart).first->hasTime())
        {
            time_ += particles_.at(iPart).first->GetVtxTime(smearing);
            ++n_time_tracks_;
        }
    }
    if(n_time_tracks_ == 0)
        return -100;

    time_ = time_ / n_time_tracks_;
    
    return time_;
}

//----------compute sumpt2 using all the particles related to the vtx---------------------
float VertexWithFT::sumPtSquared(float dz_cut, float pt_cut) const
{
    double sum = 0.;
    double pT;
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
        pT = particles_.at(iPart).first->pt();
        if(fabs(particles_.at(iPart).second) < dz_cut && pT > pt_cut)           
            sum += pT*pT;
    }
    return sum;
}

//**********Utils*************************************************************************

//----------update associated particles vtx reference-------------------------------------
void VertexWithFT::FixVtxRefs()
{
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
        if(particles_.at(iPart).first->GetRecoVtx() != this)
            particles_.at(iPart).first->SetRecoVtx(this);
    }

    return;
}        

//**********relation operators************************************************************

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
