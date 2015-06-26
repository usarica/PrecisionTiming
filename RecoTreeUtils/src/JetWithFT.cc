#include "FastTiming/RecoTreeUtils/interface/JetWithFT.h"

JetWithFT::JetWithFT()
{}

JetWithFT::JetWithFT(const reco::PFJet* pfJet, VertexWithFT* pV, float smearing,
                     vector<EcalRecHit>* recVectEK, const CaloGeometry* skGeometry,
                     const MagneticField* magField):
  reco::PFJet(*pfJet), tRes_(smearing)
{
    primaryVtx_=pV;
    corrMomentum_=TLorentzVector(0,0,0,0);

    for(auto& constituent : getPFConstituents())
    {
        if(constituent->particleId() > 4)
            continue;
        PFCandidateWithFT particle(constituent.get(), recVectEK, skGeometry, magField, NULL, primaryVtx_);
        if(!particle.hasTime() ||  
           fabs(particle.GetVtxTime(tRes_) - primaryVtx_->GetSeedRef()->GetVtxTime(tRes_))<tRes_*3)        
            corrMomentum_ += TLorentzVector(particle.px(), particle.py(), particle.pz(), particle.energy());        
        else
        {
            if(particle.particleId() == 4)
                ++nNeutralRej_;
            else
                ++nChargedRej_;
        }
    }
}
