#ifndef JetWithFT_H
#define JetWithFT_H

#include "TLorentzVector.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"

class JetWithFT : public reco::PFJet
{
public:
    //---ctors---
    JetWithFT();
    JetWithFT(const reco::PFJet* pfJet, VertexWithFT* pV, float smearing,
              vector<EcalRecHit>* ecalRecHits, const CaloGeometry* skGeometry,
              const MagneticField* magField);
    //---dtor---
    ~JetWithFT(){};

    bool operator>(const JetWithFT& j) const {return corrMomentum_.Pt()>j.GetCorrMomentum()->Pt();};
    bool operator<(const JetWithFT& j) const {return corrMomentum_.Pt()<j.GetCorrMomentum()->Pt();};

    //---getters---
    const TLorentzVector* GetCorrMomentum() const {return &corrMomentum_;};

private:
    TLorentzVector corrMomentum_;
    VertexWithFT*  primaryVtx_;
    int            nChargedRej_;
    int            nNeutralRej_;
    float          tRes_;
};

#endif
