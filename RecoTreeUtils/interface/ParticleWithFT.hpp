#ifndef ParticleWithFT_H
#define ParticleWithFT_H

#include <stdlib>
#include <iostream>

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"

using namespace std;

class ParticleWithFT
{
public:
    //---ctors---
    ParticleWithFT();
    ParticleWithFT(PFCandidate* PFCand, vector<EcalRecHit>* ecalRecHits);
    //---dtor---
    ~ParticleWithFT();
    //---getters---
    inline float GetTime() {return time_;};
    inline PFCluster* GetPFCluster() {return pfCluster_;};
    inline PFCandidate* GetPFCandidate() {return pfCand_;};
    pair<float, float> GetRecHitTimeMaxE();
    vector<pair<float, float> > GetRecHitsTimeE();
    
private:
    PFCandidate        pfCand_;
    PFCluster*         pfCluster_;
    vector<EcalRecHit> recHitColl_;
    float              clusterE_;
    float              maxRecHitE_;
    float              time_;   
};

#endif
