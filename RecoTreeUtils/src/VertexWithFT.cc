#include "FastTiming/RecoTreeUtils/interface/VertexWithFT.h"

//----------ctors-------------------------------------------------------------------------

VertexWithFT::VertexWithFT()
{}

VertexWithFT::VertexWithFT(const reco::Vertex* recoVtx):
    reco::Vertex(recoVtx->position(), recoVtx->error(), recoVtx->chi2(), recoVtx->ndof(), recoVtx->tracksSize()),
    hasSeed_(false), isGhost_(false), time_(-1000), time_EB_(-1000), time_EE_(-1000) ,
    nCharged_(0), n_time_tracks_(-1), n_time_tracks_EE_(-1), n_time_tracks_EB_(-1) 
{
    genVtxRef_ = NULL;
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

//----------Set gen vtx reference and id--------------------------------------------------
//--- id=0 => signal
//--- id>0 => PU
void VertexWithFT::SetGenVtxRef(const SimVertex* genVtx, int id)
{
    genVtxRef_ = genVtx; 
    genVtxId_ = id;

    return;
}

//----------Create gen vtx from position and set id---------------------------------------
void VertexWithFT::SetGenVtxRef(math::XYZTLorentzVector vp, int id)
{
    if(genVtxRef_)
        delete genVtxRef_;
    
    genVtxRef_ = new SimVertex(math::XYZVector(vp.x(), vp.y(), vp.z()), vp.t()*1E9);
    genVtxId_ = id;

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
        if(particle->particleId() >= 4)
            dz_tmp = 0;
        if(dz_tmp == -1000)
            dz_tmp = particle->GetTrack()->dz(this->position());

        particles_.push_back(make_pair(particle, dz_tmp));
        if(particle->particleId() < 4 && particle->pt() > 2)
            ++nCharged_;
    }
    return;
}

//----------Remove particle to the vertex-------------------------------------------------
void VertexWithFT::RemoveParticle(PFCandidateWithFT* particle)
{
}

//----------Get a reference to the seed particle------------------------------------------
PFCandidateWithFT* VertexWithFT::GetSeedRef()
{
    if(hasSeed() && !isGhost())
        return particles_.at(0).first;

    return NULL;
}

//----------Compute vertex time-----------------------------------------------------------
//---particle_type: 0 --> all
//---               1 --> charged particle
//---               2 --> neutral (photons)
float VertexWithFT::ComputeTime(int particle_type, float pt_cut, float pz2_cut, float smearing)
{
    //---return gen time if ghost
    if(isGhost())
        return time_;
            
    float seed_time=0;
    float seed_time_EB_=0;
    float seed_time_EE_=0;
    time_ = 0;
    time_EB_ = 0;
    time_EE_ = 0;
    n_time_tracks_ = 0;
    n_time_tracks_EB_ = 0;
    n_time_tracks_EE_ = 0;

    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
        //---select the particles type used to compute the vtx time
        if(particle_type == 1 && particles_.at(iPart).first->particleId() > 3)
            continue;
        if(particle_type == 2 && particles_.at(iPart).first->particleId() < 4)
            continue;

        //---loop over the selected particles
        float pt_tmp = particles_.at(iPart).first->pt();
        float pz2_tmp = pow(particles_.at(iPart).first->pz(), 2);
        if(pt_tmp > pt_cut && pz2_tmp > pz2_cut && particles_.at(iPart).first->hasTime())
        {
            if(n_time_tracks_ == 0)
            {
                seed_time = particles_.at(iPart).first->GetVtxTime(smearing, false);
                time_ += seed_time;
                ++n_time_tracks_;
            }
            else if(fabs(particles_.at(iPart).first->GetVtxTime(smearing, false) - seed_time) < smearing*2)
            {
                time_ += particles_.at(iPart).first->GetVtxTime(smearing, false);
                ++n_time_tracks_;
            }
	    //EB
	    if(fabs(particles_.at(iPart).first->eta()) < 1.47 ){
	      if(n_time_tracks_EB_ == 0)
		{
		  seed_time_EB_ = particles_.at(iPart).first->GetVtxTime(smearing, false);
		  time_EB_ += seed_time;
		  ++n_time_tracks_EB_;
		}
	      else if(fabs(particles_.at(iPart).first->GetVtxTime(smearing, false) - seed_time_EB_) < smearing*2)
		{
		  time_EB_ += particles_.at(iPart).first->GetVtxTime(smearing, false);
		  ++n_time_tracks_EB_;
		}
	    }//EB
	    else{
	      if(n_time_tracks_EE_ == 0)
		{
		  seed_time_EE_ = particles_.at(iPart).first->GetVtxTime(smearing, false);
		  time_EE_ += seed_time;
		  ++n_time_tracks_EE_;
		}
	      else if(fabs(particles_.at(iPart).first->GetVtxTime(smearing, false) - seed_time_EE_) < smearing*2)
		{
		  time_EE_ += particles_.at(iPart).first->GetVtxTime(smearing, false);
		  ++n_time_tracks_EE_;
		}
	    }//EE
        }
    }
    if(n_time_tracks_ == 0)
        return -100;

    time_ = time_ / n_time_tracks_;
    time_EB_ = time_EB_ / n_time_tracks_EB_;
    time_EE_ = time_EE_ / n_time_tracks_EE_;
    
    return time_;
}

float VertexWithFT::ComputeTimeBottomUp(int particle_type, float pt_cut, float smearing)
{
    time_ = 0;
    n_time_tracks_ = 0;
    vector<float> used_times;
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
        //---select the particles type used to compute the vtx time
        if(particle_type == 1 && particles_.at(iPart).first->particleId() > 3)
            continue;
        if(particle_type == 2 && particles_.at(iPart).first->particleId() < 4)
            continue;
        //---loop over the selected particles
        float pt_tmp = particles_.at(iPart).first->pt();
        if(pt_tmp > pt_cut && particles_.at(iPart).first->hasTime())
        {
            float tmp_time = particles_.at(iPart).first->GetVtxTime(smearing);
            used_times.push_back(tmp_time);
            time_ += tmp_time;
            ++n_time_tracks_;
        }
    }
    time_ = time_ / n_time_tracks_;

    if(n_time_tracks_ < 3)
        return time_;

    // int bad_particle=-1;
    // float new_time=time_;            
    // do
    // {
    //     time_ = new_time;
    //     float minDeltaRMS=100;
    //     for(unsigned int i=0; i<used_index.size(); ++i)
    //     {
    //         float tRMS=0;
    //         float sRMS=0;
    //         float tmp_time=((time_*n_time_tracks_) -
    //                         particles_.at(used_index.at(i)).first->GetVtxTime(smearing))/(n_time_tracks_-1);
    //         for(unsigned int j=0; j<used_index.size(); ++j)
    //         {
    //             tRMS += pow(time_ - particles_.at(used_index.at(j)).first->GetVtxTime(smearing), 2);
    //             if(j != i)
    //                 sRMS += pow(tmp_time - particles_.at(used_index.at(j)).first->GetVtxTime(smearing), 2);
    //         }
    //         tRMS = sqrt(1/(n_time_tracks_-1)*tRMS);
    //         sRMS = sqrt(1/(n_time_tracks_-2)*sRMS);
    //         float tmpDeltaRMS=tRMS*(sqrt(n_time_tracks_/(n_time_tracks_-1))) - sRMS;
    //         if(tmpDeltaRMS > 0 && tmpDeltaRMS < minDeltaRMS)
    //         {
    //             new_time = tmp_time;
    //             bad_particle = i;
    //         }
    //     }
    //     if(bad_particle != -1)
    //     {
    //         used_index.erase(used_index.begin()+bad_particle);
    //         --n_time_tracks_;
    //     }
    // }
    // while(n_time_tracks_ > 2 && bad_particle!=-1);

    return time_;
}                      

//----------compute sumpt2 using all the particles related to the vtx---------------------
float VertexWithFT::sumPtSquared(float dz_cut, float pt_cut, int EB) const
{
    double sum = 0.;
    double pT;
    for(unsigned int iPart=0; iPart<particles_.size(); ++iPart)
    {
      if((EB == 1 && fabs(particles_.at(iPart).first->eta()) < 1.47) ||
	 (EB == 0 && fabs(particles_.at(iPart).first->eta()) >= 1.47) ||
	 EB == -1){
        pT = particles_.at(iPart).first->pt();
        if(fabs(particles_.at(iPart).second) < dz_cut && pT > pt_cut)           
	  sum += pT*pT;
      }
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
