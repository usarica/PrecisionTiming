#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"

//**********Ctors*************************************************************************

PFCandidateWithFT::PFCandidateWithFT():
    clusterE_(0), maxRecHitE_(0), time_(0), vtxTime_(0)
{}

PFCandidateWithFT::PFCandidateWithFT(const reco::PFCandidate* PFCand,
                               vector<EcalRecHit>* ecalRecHits, float vtxTime):
    reco::PFCandidate(*PFCand), clusterE_(0), maxRecHitE_(0), time_(0), vtxTime_(vtxTime)
{
    pfCand_ = PFCand;
    pfCluster_ = NULL;
    recHitColl_ = *ecalRecHits;   
    //---get the right ecal cluster---
    float min_dist_cluster=100;
    for(auto& blockPair : pfCand_->elementsInBlocks())
    {        
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
        if(blockElement.type() != 4 || !blockElement.clusterRef().isAvailable())
            continue;
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
    for(unsigned int iRec=0; iRec<recHitColl_.size(); iRec++)
    {
	for(unsigned int iDet=0; iDet<detIdMap.size(); iDet++)
	{
	    if(detIdMap.at(iDet).first == recHitColl_.at(iRec).id())
	    {
                    TandE_vect.push_back(make_pair(
                                           recHitColl_.at(iRec).time() - vtxTime_ * 1E9,
                                           recHitColl_.at(iRec).energy()));
	    }
	}
    }
    return TandE_vect;
}
