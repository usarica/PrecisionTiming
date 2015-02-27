#ifndef PFCandidateWithFT_H
#define PFCandidateWithFT_H

//****************************************************************************************
// lengths are in cm, times are in ns
//
//
//
//****************************************************************************************

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
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

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
                      const SimVertex* genVtx, const reco::Vertex* recoVtx,
                      const CaloGeometry* skGeometry=NULL, const MagneticField* magField=NULL);
    //---dtor---
    ~PFCandidateWithFT();
    //---getters---
    inline float                    GetTime() {return absTime_;};
    inline const reco::PFCluster*   GetPFCluster() {return pfCluster_;};
    inline const reco::PFCandidate* GetPFCandidate() {return pfCand_;};
    inline const reco::Track*       GetTrack() {return recoTrack_;};
    inline math::XYZVector          GetRecoVtxPos() {return recoVtxPos_;};
    inline float                    GetDrTrackCluster() { return drTrackCluster_; };   
    inline float                    GetTOF() { return GetPropagatedTrackLength()/3E10*1E9; };
    inline pair<float, float>       GetRecHitTimeMaxE() {return GetRecHitTimeE(ecalSeed_);};
    float                           GetTrackLength();
    float                           GetPropagatedTrackLength();
    float                           GetGenTOF();
    void                            GetTrackInfo(float& alpha, float& trackR, float& secant, int& charge);
    pair<float, float>              GetRecHitTimeE(DetId id);
    vector<pair<float, float> >     GetRecHitsTimeE();
    //---utils---
    inline bool                     hasTime() {if(pfCluster_) return true; return false;};
    void                            TrackReconstruction();
    

private:
    const reco::PFCandidate* pfCand_;
    const reco::PFCluster*   pfCluster_;
    const reco::Vertex*      recoVtx_;
    const MagneticField*     magField_;
    const CaloGeometry*      skGeometry_;
    const SimVertex*         genVtx_;
    vector<EcalRecHit>*      recHitColl_;
    DetId                    ecalSeed_;
    float                    clusterE_;
    float                    maxRecHitE_;
    float                    absTime_;
    float                    vtxTime_;
    const reco::Track*       recoTrack_;
    math::XYZVector          ecalPos_;
    math::XYZVector          genVtxPos_;
    math::XYZVector          recoVtxPos_;
    math::XYZVector          secant_;
    float                    alpha_;
    float                    trackPt_;
    float                    trackR_;
    float                    trackL_;
    float                    propagatedTrackL_;
    float                    drTrackCluster_;
};

#endif
