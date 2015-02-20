#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"

//**********Ctors*************************************************************************

PFCandidateWithFT::PFCandidateWithFT():
    clusterE_(0), maxRecHitE_(0), time_(0), vtxTime_(0)
{}

PFCandidateWithFT::PFCandidateWithFT(const reco::PFCandidate* PFCand,
                               vector<EcalRecHit>* ecalRecHits, float vtxTime):
  reco::PFCandidate(*PFCand), clusterE_(0), maxRecHitE_(0), time_(0), vtxTime_(vtxTime), 
  innerP_(0,0,0), outerP_(0,0,0), secant_(0,0,0), alpha_(0), trackPt_(0), trackR_(0), trackL_(0), drTrackCluster_(0)
{
    pfCand_ = PFCand;
    pfCluster_ = NULL;
    recHitColl_ = *ecalRecHits;   
    recoTrack_ = NULL;
    //    innerMomentum_(0), outerMomentum_(0)

    //---get the right ecal cluster---
    float min_dist_cluster = 100;
    float min_dist_track = 100;
    for(auto& blockPair : pfCand_->elementsInBlocks())
    {        
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
	//Cluster
        // if(blockElement.type() != 4 || !blockElement.clusterRef().isAvailable())
        //     continue;
	if(blockElement.type() == 4 && blockElement.clusterRef().isAvailable()){
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
	    }
	}

	//Track
        if(blockElement.type() == 1 && blockElement.trackRef().isAvailable()){
	  
	  reco::Track tmpTrack = *blockElement.trackRef().get();	 
	  // Characteristics of the track
	  const reco::PFBlockElementTrack& elementTrack = dynamic_cast<const reco::PFBlockElementTrack &>(blockElement);
	  const math::XYZPoint pfTrackPos(elementTrack.positionAtECALEntrance().x(), 
					  elementTrack.positionAtECALEntrance().y(), 
					  elementTrack.positionAtECALEntrance().z());
	  float tmp_dist = DeltaR(pfTrackPos.eta(), pfCand_->eta(),
				  pfTrackPos.phi(), pfCand_->phi());
	  drTrackCluster_ = tmp_dist;
	  if(tmp_dist < min_dist_track)
	    {
	      min_dist_track = tmp_dist;

	      recoTrack_ = &tmpTrack;	 
	      innerP_ = math::XYZVector(tmpTrack.innerMomentum().x(), tmpTrack.innerMomentum().y(), tmpTrack.innerMomentum().z());
	      outerP_ = math::XYZVector(tmpTrack.outerMomentum().x(), tmpTrack.outerMomentum().y(), tmpTrack.outerMomentum().z());  // if ok => can be removed
	      secant_ = math::XYZVector(pfTrackPos.x()-innerP_.x(), pfTrackPos.y()-innerP_.y(), pfTrackPos.z()-innerP_.z());
	      if(pfCand_->charge() == 1) alpha_ = ROOT::Math::VectorUtil::Angle(innerP_, secant_);
	      if(pfCand_->charge() == -1) alpha_ = ROOT::Math::VectorUtil::Angle(innerP_, innerP_);

	      trackPt_ = elementTrack.trackRef()->pt();
	      trackR_ = trackPt_ / 0.3 / 3.8;
	      trackL_ = 2 * alpha_ * trackR_;
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

//----------     --------------------------
void PFCandidateWithFT::GetTrackInfo(float& phiIn, float& phiOut, float& alpha, float& trackR, float& secant, int& charge)
{
    phiIn = innerP_.phi();
    phiOut = outerP_.phi();
    alpha = alpha_;
    trackR = trackR_;
    secant = secant;
    charge = pfCand_->charge();
    return;
}

