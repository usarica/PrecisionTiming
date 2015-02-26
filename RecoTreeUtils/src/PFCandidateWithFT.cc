#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"

//**********Ctors*************************************************************************

PFCandidateWithFT::PFCandidateWithFT():
    clusterE_(0), maxRecHitE_(0), time_(0), vtxTime_(0)
{}

PFCandidateWithFT::PFCandidateWithFT(const reco::PFCandidate* PFCand,
                               vector<EcalRecHit>* ecalRecHits, const SimVertex* primaryVtx):
  reco::PFCandidate(*PFCand), clusterE_(0), maxRecHitE_(0), time_(0), vtxTime_(primaryVtx->position().t()), 
  vtxPos_(0,0,0), secant_(0,0,0), alpha_(0), trackR_(0), trackL_(-1)
{
    pfCand_ = PFCand;
    primaryVtx_ = primaryVtx;
    recHitColl_ = *ecalRecHits;   
    pfCluster_ = NULL;
    recoTrack_ = NULL;

    //---get the right ecal cluster---
    float min_dist_cluster = 100;
    for(auto& blockPair : pfCand_->elementsInBlocks())
    {        
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
	if(blockElement.type() == 4 && blockElement.clusterRef().isAvailable())
        {
            reco::PFCluster tmpCluster = *blockElement.clusterRef().get();	 
            tmpCluster.calculatePositionREP();
            REPPoint pfClusterPos = tmpCluster.positionREP();
            float tmp_dist=DeltaR(pfClusterPos.Eta(), pfCand_->eta(),
                                  pfClusterPos.Phi(), pfCand_->phi());
            if(tmp_dist < min_dist_cluster)
	    {
                min_dist_cluster = tmp_dist;
                pfCluster_ = blockElement.clusterRef().get();
                ecalSeed_ = pfCluster_->seed();
                time_ = GetRecHitTimeMaxE().first + vtxTime_ + GetGenTOF();
	    }
	}
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
    for(unsigned int iRec=0; iRec<recHitColl_.size(); iRec++)
    {
        if(recHitColl_.at(iRec).id() == id)
        {
                return make_pair(recHitColl_.at(iRec).time() - vtxTime_ * 1E9,
                                 recHitColl_.at(iRec).energy());
        }
    }
    //---if not found return time=0, energy=-1---
    return make_pair(0, -1);
}

//----------Return <time, enegy> of all the ecal cluster recHits--------------------------
vector<pair<float, float> > PFCandidateWithFT::GetRecHitsTimeE()
{
    vector<pair<float, float> > TandE_vect;
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
        for(unsigned int iRec=rh_start; iRec<recHitColl_.size(); iRec++)
        {
	    if(sortedDetId.at(iDet) == recHitColl_.at(iRec).id())
	    {
                rh_start=iRec+1;
                TandE_vect.push_back(make_pair(
                                         recHitColl_.at(iRec).time() - vtxTime_ * 1E9,
                                         recHitColl_.at(iRec).energy()));
                break;
	    }
	}
    }
    return TandE_vect;
}

//----------TOF wrt sim vertex------------------------------------------------------------
float PFCandidateWithFT::GetGenTOF()
{
    return 1;
}

//----------Get track length--------------------------------------------------------------
float PFCandidateWithFT::GetTrackLength()
{
    if(pfCand_->particleId() < 4)
    {
        if(trackL_ == -1)
            TrackReconstruction();
        return trackL_;

    }
    return trackL_;
}

//----------Simple track info getter------------------------------------------------------
void PFCandidateWithFT::GetTrackInfo(float& alpha, float& trackR, float& secant, int& charge)
{
    alpha = alpha_;
    trackR = trackR_;
    secant = secant;
    charge = pfCand_->charge();
    return;
}

//----------Track reconstruction----------------------------------------------------------
void PFCandidateWithFT::TrackReconstruction()
{
    float min_dist_track = 100;
    for(auto& blockPair : pfCand_->elementsInBlocks())
    {
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
        if(blockElement.type() == 1 && blockElement.trackRef().isAvailable())
        {
            reco::Track tmpTrack = *blockElement.trackRef().get();	 
            //---Characteristics of the track
            const reco::PFBlockElementTrack& elementTrack = dynamic_cast<const reco::PFBlockElementTrack &>(blockElement);
            const math::XYZPoint pfTrackPos(elementTrack.positionAtECALEntrance().x(), 
                                            elementTrack.positionAtECALEntrance().y(),
                                            elementTrack.positionAtECALEntrance().z());

            float tmp_dist = DeltaR(pfTrackPos.eta(), pfCand_->eta(),
                                    pfTrackPos.phi(), pfCand_->phi());
            if(tmp_dist < min_dist_track)
            {
                min_dist_track = tmp_dist;

                recoTrack_ = &tmpTrack;	 
                vtxPos_ = math::XYZVector(primaryVtx_->position().x(),
                                          primaryVtx_->position().y(),
                                          primaryVtx_->position().z());
                secant_ = math::XYZVector(pfTrackPos.x()-vtxPos_.x(), pfTrackPos.y()-vtxPos_.y(), 0);              
                trackPt_ = elementTrack.trackRef()->pt();
                trackR_ = trackPt_ / 0.3 / 3.8;
                alpha_ = asin(secant_.R() / (100*2*trackR_));
                trackL_ = sqrt(pow(2*alpha_*trackR_, 2) + pow((pfTrackPos.z()-vtxPos_.z())/100, 2));

                // DEBUG: controllare (con pietro magari) la differenza tra le diverse approssimazioni
                //        utili per calcolare la traccia.
                //trackR_ = (trackL_ - (fabs(pfTrackPos.z()-vtxPos_.z()))/100/cos(asin(pfCand_->pt()/pfCand_->p())))/3E8;
                //trackR_ = (2*alpha_*trackR_ - secant_.R()/100)/3E8;
            }
        }
    }
    return;
}
