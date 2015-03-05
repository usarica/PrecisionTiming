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
#include "TRandom.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToBase.h" 
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

#include "FastTiming/RecoTreeUtils/interface/VertexWithFT.h"

typedef edm::RefToBase<reco::Track> TrackBaseRef;
typedef ROOT::Math::PositionVector3D<ROOT::Math::CylindricalEta3D<Double32_t> > REPPoint;

using namespace std;

class VertexWithFT;

//****************************************************************************************

class PFCandidateWithFT
{
public:
    //---ctors---
    PFCandidateWithFT();
    PFCandidateWithFT(const reco::PFCandidate* PFCand, vector<EcalRecHit>* ecalRecHits,
                      const SimVertex* genVtx, VertexWithFT* recoVtx=NULL,
                      const CaloGeometry* skGeometry=NULL, const MagneticField* magField=NULL);
    //---dtor---
    ~PFCandidateWithFT();
    //---getters---
    inline const reco::PFCluster*   GetPFCluster() {return pfCluster_;};
    inline const reco::PFCandidate* GetPFCandidate() const {return pfCand_;}; 
    inline VertexWithFT*            GetRecoVtx() {return recoVtx_;};
    inline math::XYZVector          GetRecoVtxPos() {return recoVtxPos_;};
    inline float                    GetDrTrackCluster() { return drTrackCluster_; };
    inline float                    GetRawTime() {return rawTime_;};
    inline float                    GetTOF() { return GetPropagatedTrackLength()/3E10*1E9; };
    inline pair<float, float>       GetRecHitTimeMaxE() {return GetRecHitTimeE(ecalSeed_);};
    const reco::Track*              GetTrack();
    float                           GetTrackLength();
    float                           GetPropagatedTrackLength();
    float                           GetGenTOF();
    float                           GetECALTime(float smearing=0);
    float                           GetVtxTime(float smearing=0);
    pair<float, float>              GetRecHitTimeE(DetId id);
    vector<pair<float, float> >     GetRecHitsTimeE();
    //---setters---
    void                            SetRecoVtx(VertexWithFT* recoVtx);
    //---utils---
    inline bool                     hasTime() {if(pfCluster_) return true; return false;};
    void                            TrackReconstruction();
    

private:
    const reco::PFCandidate* pfCand_;
    const reco::PFCluster*   pfCluster_;
    const MagneticField*     magField_;
    const CaloGeometry*      skGeometry_;
    const SimVertex*         genVtx_;
    VertexWithFT*            recoVtx_;
    vector<EcalRecHit>*      recHitColl_;
    DetId                    ecalSeed_;
    float                    clusterE_;
    float                    maxRecHitE_;
    float                    rawTime_;
    float                    tSmearing_;
    float                    smearedRawTime_;
    float                    vtxTime_;
    const reco::Track*       recoTrack_;
    math::XYZVector          ecalPos_;
    math::XYZVector          genVtxPos_;
    math::XYZVector          recoVtxPos_;
    float                    trackPt_;
    float                    trackL_;
    float                    propagatedTrackL_;
    float                    drTrackCluster_;
};

#endif
