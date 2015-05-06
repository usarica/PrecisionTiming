#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"

typedef std::vector<reco::TrackBaseRef >::const_iterator trackRef_iterator;

//**********Ctors*************************************************************************

PFCandidateWithFT::PFCandidateWithFT()
{}

PFCandidateWithFT::PFCandidateWithFT(const reco::PFCandidate* PFCand, vector<EcalRecHit>* ecalRecHits,
                                     const CaloGeometry* skGeometry, const MagneticField* magField,
                                     const SimVertex* genVtx, VertexWithFT* recoVtx):
    reco::PFCandidate(*PFCand), clusterE_(0), rawTime_(0), ecalPos_(0,0,0),
    recoVtxPos_(0,0,0), trackL_(-1), propagatedTrackL_(-1), drTrackCluster_(-1)
{
    hasTime_ = false;
    pfCand_ = PFCand;
    magField_ = magField;
    skGeometry_ = skGeometry;
    recHitColl_ = ecalRecHits;   
    genVtx_ = genVtx;
    recoVtx_ = recoVtx;
    pfCluster_ = NULL;
    recoTrack_ = NULL;
    smearedRawTime_=0;
    mvaRawTime_=0;
    tSmearing_=-1;
    tSmearingMVA_=-1;
    //---reco vtx info---
    if(recoVtx_)
        recoVtxPos_ = math::XYZVector(recoVtx_->position().x(),
                                      recoVtx_->position().y(),
                                      recoVtx_->position().z());
    //---get the right ecal cluster---
    float min_dist_cluster = 100;
    for(auto& blockPair : elementsInBlocks())
    {        
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
	if(blockElement.type() == 4 && blockElement.clusterRef().isAvailable())
        {
            reco::PFCluster tmpCluster = *blockElement.clusterRef().get();	 
            tmpCluster.calculatePositionREP();
            pfClusterPos_ = tmpCluster.positionREP();
            float tmp_dist=deltaR(pfClusterPos_.Eta(), pfClusterPos_.Phi(),
                                  positionAtECALEntrance().eta(), positionAtECALEntrance().phi());
            if(tmp_dist < min_dist_cluster)
            {
                min_dist_cluster = tmp_dist;
                pfCluster_ = blockElement.clusterRef().get();
	    }
	}
    }
    if(pfCluster_)
    {
        FindEcalSeed();
        const CaloCellGeometry* cell=skGeometry_->getGeometry(ecalSeed_);
        // if(pfCluster_->isEB()
        //     ecalPos_ = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(3.5);
        // else
            ecalPos_ = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(10*0.4-0.075-0.25);
        rawTime_ = GetRecHitTimeE(ecalSeed_).first + GetGenTOF();
        if(GetRecHitTimeMaxE().second != -1)
            hasTime_ = true;
        recoVtx_ = NULL;
    }
}

//**********Dtor**************************************************************************

PFCandidateWithFT::~PFCandidateWithFT()
{}

//**********Getters***********************************************************************

//----------Return <time, enegy> of a specific ecal cluster recHit------------------------
pair<float, float> PFCandidateWithFT::GetRecHitTimeE(DetId id)
{
    //---search for the right recHit---
    for(unsigned int iRec=0; iRec<recHitColl_->size(); iRec++)
    {
        if(recHitColl_->at(iRec).id() == id)
        {
            return make_pair(recHitColl_->at(iRec).time(), recHitColl_->at(iRec).energy());
        }
    }
    //---if not found return time=-1, energy=-1---
    return make_pair(-1, -1);
}

//----------Return <time, enegy> of all the ecal cluster recHits--------------------------
vector<FTEcalRecHit>* PFCandidateWithFT::GetRecHits()
{
    if(ftRecHits_.size() == 0)
    {
        vector<pair<DetId, float> > detIdMap = pfCluster_->hitsAndFractions();
        //---Sort cluster rechits---
        vector<DetId> sortedDetId;
        sortedDetId.push_back(detIdMap.at(0).first);
        for(unsigned int i=1; i<detIdMap.size(); i++)        
        {
            bool inserted=false;
            for(unsigned int j=0; j<sortedDetId.size(); j++)
            {
                if(detIdMap.at(i).first.rawId() < sortedDetId.at(j).rawId())
                {
                    sortedDetId.insert(sortedDetId.begin()+j, detIdMap.at(i).first);
                    inserted=true;
                    break;
                }
            }
            if(!inserted)
                sortedDetId.push_back(detIdMap.at(i).first);
        }
        //---search for the right rechits
        unsigned int rh_start=0;
        for(unsigned int iDet=0; iDet<sortedDetId.size(); iDet++)
        {
            for(unsigned int iRec=rh_start; iRec<recHitColl_->size(); iRec++)
            {
                if(sortedDetId.at(iDet) == recHitColl_->at(iRec).id())
                {
                    rh_start=iRec+1;
                    const CaloCellGeometry* cell=skGeometry_->getGeometry(sortedDetId.at(iDet));
                    GlobalPoint recHitPos = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0);
                    FTEcalRecHit tmp(sortedDetId.at(iDet).rawId(),
                                     EKDetId(sortedDetId.at(iDet)).ix(),
                                     EKDetId(sortedDetId.at(iDet)).iy(),
                                     recHitPos.z(),
                                     recHitColl_->at(iRec).time()+recHitPos.mag()/30,
                                     recHitColl_->at(iRec).energy());
                    ftRecHits_.push_back(tmp);
                    break;
                }
            }
        }
    }
    //---prepare the FTRecHits for the MVA
    sort(ftRecHits_.begin(), ftRecHits_.end());
    reverse(ftRecHits_.begin(), ftRecHits_.end());

    return &ftRecHits_;
}

//----------TOF wrt nominal IP------------------------------------------------------------
//---NOTE: the position wrt the TOF is calculated in CMSSW is wrong!
//---      3+0.5 is for the standard EE not SK
float PFCandidateWithFT::GetGenTOF()
{
    const CaloCellGeometry* cell=skGeometry_->getGeometry(ecalSeed_);
    math::XYZVector ecalBadPos(dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(3+0.5));
    
    //---(distance/c)*1E9
    return ecalBadPos.R()/30;
}

float PFCandidateWithFT::GetTOF(tof_algo method)
{
    if(method == pzTOF)
        return p()/(fabs(pz()))*fabs(ecalPos_.z()-recoVtxPos_.z())/30;
    
    return GetPropagatedTrackLength()/3E10*1E9;
};

//----------Apply a time resolution smearing to the raw time------------------------------
//---NOTE: smearing must be in ns
float PFCandidateWithFT::GetECALTime(float smearing)
{
    //---do not re-smear the raw time
    if(tSmearing_ == -1)
    {
        //---unique rndm seed
        TRandom rndm(0);
        tSmearing_ = smearing;
        smearedRawTime_ = rndm.Gaus(rawTime_, smearing);
    }
    
    return smearedRawTime_;
}

//----------Get the time smeared and evaluated from the mva-------------------------------
float PFCandidateWithFT::GetECALTimeMVA(float smearing)
{
    //---compute ecal time from the mva
    if(mvaComputer_ && tSmearingMVA_ == -1)
    {
        tSmearingMVA_ = smearing;
        mvaRawTime_ = mvaComputer_->GetMVATime(GetRecHits(), energy(), pt(), pz(), smearing);
    }
    
    return mvaRawTime_;
}

//----------Get TOF corrected vtx time----------------------------------------------------
float PFCandidateWithFT::GetVtxTime(float smearing, bool mva, tof_algo tof_method)
{
    if(mva)
    {
        GetECALTimeMVA(smearing);
        float tof_from_face = p()/(fabs(pz()))*fabs(GetRecHits()->at(0).z-recoVtxPos_.z())/30;
        return mvaRawTime_ - tof_from_face;
    }
    else
        return GetECALTime(smearing) - GetTOF(tof_method);
}

//----------Get the right track ref from the PFBlock--------------------------------------
const reco::Track* PFCandidateWithFT::GetTrack()
{
    if(recoTrack_)
        return recoTrack_;
    float min_dist_track = 100;
    //---get the track from the PFBlock
    for(auto& blockPair : elementsInBlocks())
    {
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
        if(blockElement.type() == 1 && blockElement.trackRef().isAvailable())
        {
            recoTrack_ = blockElement.trackRef().get();	 
            float tmp_dist = deltaR(recoTrack_->eta(), recoTrack_->phi(), eta(), phi());
            if(tmp_dist < min_dist_track)
                min_dist_track = tmp_dist;
        }
    }
    if(recoTrack_ && pfCluster_)
        drTrackCluster_ = deltaR(recoTrack_->outerEta(), recoTrack_->outerPhi(),
                                 pfClusterPos_.Eta(), pfClusterPos_.Phi());
    return recoTrack_;
}

//----------Get track length--------------------------------------------------------------
float PFCandidateWithFT::GetTrackLength()
{
    if(trackL_ == -1)
        TrackReconstruction();

    return trackL_;
}

//----------Get track length using CMSSW tools--------------------------------------------

float PFCandidateWithFT::GetPropagatedTrackLength()
{
    //---if charged compute the track length exploiting HelixPropagator
    if(particleId() < 4)
    {
        if(!recoTrack_)
            GetTrack();
        if(recoTrack_)
        {
            GlobalPoint startingPoint(recoVtxPos_.x(), recoVtxPos_.y(), recoVtxPos_.z());
            GlobalVector startingMomentum(px(), py(), pz());
            GlobalPoint endPoint(ecalPos_.x(), ecalPos_.y(), ecalPos_.z());
            FreeTrajectoryState trajectory(startingPoint, startingMomentum, recoTrack_->charge(), magField_);
            SteppingHelixPropagator propagator(magField_);
            propagatedTrackL_ = propagator.propagateWithPath(trajectory, endPoint).second;            
        }
    }
    //---if neutral return the lenght of the straight line
    else
        propagatedTrackL_ = GetTrackLength();
    
    return propagatedTrackL_;
}

//**********Setters***********************************************************************

//----------Set primary reco vtx reference------------------------------------------------
void PFCandidateWithFT::SetRecoVtx(VertexWithFT* recoVtx)
{
    if(!recoVtx)
        return;
    recoVtx_ = recoVtx;
    recoVtxPos_ = math::XYZVector(recoVtx_->position().x(),
                                  recoVtx_->position().y(),
                                  recoVtx_->position().z());
    trackL_ = -1;
    propagatedTrackL_ = -1;
}

//**********Utils*************************************************************************

//----------Search for the impact RecHit--------------------------------------------------
DetId PFCandidateWithFT::FindEcalSeed()
{
    ecalSeed_ = pfCluster_->seed();
    float maxE=0;
    vector<pair<DetId, float> > detIdMap = pfCluster_->hitsAndFractions();
    
    for(vector<pair<DetId, float> >::iterator it=detIdMap.begin(); it!=detIdMap.end(); ++it)
    {
        float tmpE = GetRecHitTimeE(it->first).second;
        if(tmpE > maxE) 
        {
            maxE = tmpE;
            ecalSeed_ = it->first;
        }
    }

    return ecalSeed_;    
}

//----------Track reconstruction----------------------------------------------------------
void PFCandidateWithFT::TrackReconstruction()
{
    if(particleId() < 4)
    {
        if(!recoTrack_)
            GetTrack();
    
        math::XYZVector secant(ecalPos_.x()-recoVtxPos_.x(), ecalPos_.y()-recoVtxPos_.y(), 0);              
        trackPt_ = recoTrack_->pt();
        float trackR = trackPt_*100 / 0.3 / 3.8;
        float alpha = asin(secant.R() / (2*trackR));
        trackL_ = sqrt(pow(2*alpha*trackR, 2) + pow(ecalPos_.z()-recoVtxPos_.z(), 2));
    }
    else
        trackL_ = (ecalPos_ - recoVtxPos_).R();
}

