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
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "FastTiming/RecoTreeUtils/interface/VertexWithFT.h"
#include "FastTiming/RecoTreeUtils/interface/MVATimeComputer.h"

typedef edm::RefToBase<reco::Track> TrackBaseRef;
typedef ROOT::Math::PositionVector3D<ROOT::Math::CylindricalEta3D<Double32_t> > REPPoint;

using namespace std;

class VertexWithFT;

typedef enum tof_algo{
    propagatedTrack,
    pzTOF,
} tof_algo;

//****************************************************************************************

class PFCandidateWithFT : public reco::PFCandidate
{
public:
    //---ctors---
    PFCandidateWithFT();
    PFCandidateWithFT(const reco::PFCandidate* PFCand, vector<EcalRecHit>* ecalRecHits,
                      const CaloGeometry* skGeometry_, const MagneticField* magField,
                      const SimVertex* genVtx=NULL, VertexWithFT* recoVtx=NULL);
    //---dtor---
    ~PFCandidateWithFT();
    //---getters---
    inline const reco::PFCandidate* GetPFCandidateRef() const {return pfCand_;};
    inline const reco::PFCluster*   GetPFCluster() {return pfCluster_;};
    inline REPPoint                 GetPFClusterPos() {return pfClusterPos_;};                 
    inline VertexWithFT*            GetRecoVtx() {return recoVtx_;};
    inline math::XYZVector          GetRecoVtxPos() {return recoVtxPos_;};
    inline float                    GetDrTrackCluster() { return drTrackCluster_; };
    inline float                    GetRawTime() {return rawTime_;};
    inline pair<float, float>       GetRecHitTimeMaxE() {return GetRecHitTimeE(ecalSeed_);};
    const reco::Track*              GetTrack();
    float                           GetTrackLength();
    float                           GetPropagatedTrackLength();
    float                           GetGenTOF();
    float                           GetTOF(tof_algo method);
    float                           GetECALTime(float smearing=0);
    float                           GetECALTimeMVA(float smearing=0);
    float                           GetVtxTime(float smearing=0, bool mva=false, tof_algo tof_method=propagatedTrack);
    pair<float, float>              GetRecHitTimeE(DetId id);
    vector<FTEcalRecHit>*           GetRecHits();
    //---setters---
    void                            SetRecoVtx(VertexWithFT* recoVtx);
    void                            SetMVAComputer(MVATimeComputer* mvaComputer) {mvaComputer_ = mvaComputer;};
    //---utils---
    inline bool                     hasTime() {return hasTime_;};
    DetId                           FindEcalSeed();
    void                            TrackReconstruction();    

private:
    const reco::PFCandidate* pfCand_;
    const reco::PFCluster*   pfCluster_;    
    const MagneticField*     magField_;
    const CaloGeometry*      skGeometry_;
    const SimVertex*         genVtx_;
    VertexWithFT*            recoVtx_;
    vector<EcalRecHit>*      recHitColl_;
    vector<FTEcalRecHit>     ftRecHits_;
    REPPoint                 pfClusterPos_;    
    DetId                    ecalSeed_;
    bool                     hasTime_;
    float                    clusterE_;
    float                    rawTime_;
    float                    tSmearing_;
    float                    tSmearingMVA_;
    float                    smearedRawTime_;
    float                    mvaRawTime_;
    const reco::Track*       recoTrack_;
    math::XYZVector          ecalPos_;
    math::XYZVector          recoVtxPos_;
    float                    trackPt_;
    float                    trackL_;
    float                    propagatedTrackL_;
    float                    drTrackCluster_;
    MVATimeComputer*         mvaComputer_;
};

#endif
