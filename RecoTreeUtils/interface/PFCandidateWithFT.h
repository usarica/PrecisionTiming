#ifndef PFCandidateWithFT_H
#define PFCandidateWithFT_H

#include <vector>

#include "TMath.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "Math/GenVector/VectorUtil.h"

#include "FastTiming/RecoTreeUtils/interface/Utils.h"

using namespace std;

//****************************************************************************************
class PFCandidateWithFT: public reco::PFCandidate
{
public:
    //---ctors---
    PFCandidateWithFT();
    PFCandidateWithFT(const reco::PFCandidate* PFCand, vector<EcalRecHit>* ecalRecHits,
                      const SimVertex* primaryVtx);
    //---dtor---
    ~PFCandidateWithFT();
    //---getters---
    inline float GetTime() {return time_;};
    inline const reco::PFCluster* GetPFCluster() {return pfCluster_;};
    inline const reco::PFCandidate* GetPFCandidate() {return pfCand_;};    
    inline float GetTrackR() { return trackPt_ / 0.3 / 3.8; };
    inline float GetDrTrackCluster() { return drTrackCluster_; };   
    inline float GetTrackTOF() { return trackL_/3E8; };
    inline pair<float, float> GetRecHitTimeMaxE() {return GetRecHitTimeE(ecalSeed_);};
    float GetTrackLength();
    void GetTrackInfo(float& alpha, float& trackR, float& secant, int& charge);
    pair<float, float> GetRecHitTimeE(DetId id);
    vector<pair<float, float> > GetRecHitsTimeE();
    //---utils---
    inline bool hasTime() {if(pfCluster_) return true; return false;};
    void TrackReconstruction();
    

private:
    const reco::PFCandidate* pfCand_;
    const reco::PFCluster*   pfCluster_;
    const SimVertex*         primaryVtx_;
    vector<EcalRecHit>       recHitColl_;
    DetId                    ecalSeed_;
    float                    clusterE_;
    float                    maxRecHitE_;
    float                    time_;
    float                    vtxTime_;
    const reco::Track*       recoTrack_;
    math::XYZVector          vtxPos_;
    math::XYZVector          secant_;
    float                    alpha_;
    float                    trackPt_;
    float                    trackR_;
    float                    trackL_;
    float                    drTrackCluster_;
};

#endif
