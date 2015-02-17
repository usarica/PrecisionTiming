#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"

//**********Ctors*************************************************************************

PFCandidateWithFT::PFCandidateWithFT():
    clusterE_(0), maxRecHitE_(0), time_(0), vtxTime_(0)
{}

PFCandidateWithFT::PFCandidateWithFT(const reco::PFCandidate* PFCand,
                               vector<EcalRecHit>* ecalRecHits, float vtxTime):
  reco::PFCandidate(*PFCand), clusterE_(0), maxRecHitE_(0), time_(0), vtxTime_(vtxTime), innerMomentum_(0,0,0,0), outerMomentum_(0,0,0,0), Momentum_(0)
{
    pfCand_ = PFCand;
    pfCluster_ = NULL;
    recHitColl_ = *ecalRecHits;   
    Track_ = NULL;
    //    innerMomentum_(0), outerMomentum_(0)

    //---get the right ecal cluster---
    float min_dist_cluster=100;
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

	std::cout << " costruttore " << std::endl;
	//Track
	std::cout <<"blockElement.type() = " << blockElement.type() << std::endl;
	std::cout <<"blockElement.trackRefPF().isAvailable() = " << blockElement.trackRef().isAvailable() << std::endl;
        if(blockElement.type() == 1 && blockElement.trackRef().isAvailable()){
	  
	  //	  reco::Track tmpTrack = *blockElement.trackRef().get();
	  reco::Track tmpTrack = *blockElement.trackRef().get();	 
	  // Characteristics of the track
	  const reco::PFBlockElementTrack& elementTrack = dynamic_cast<const reco::PFBlockElementTrack &>(blockElement);
	  TLorentzVector pfTrackPos(elementTrack.positionAtECALEntrance().x(), elementTrack.positionAtECALEntrance().y(), elementTrack.positionAtECALEntrance().z(), 0.);
	  float tmp_dist=DeltaR(pfTrackPos.Eta(), pfCand_->eta(),
				pfTrackPos.Phi(), pfCand_->phi());
	  std::cout << " tmp_dist = " << tmp_dist << std::endl;
	  if(tmp_dist < min_dist_cluster)
	    {
	      /*
	      double p = et.trackRef()->p();  
	      double pt = et.trackRef()->pt(); 
	      double eta = et.trackRef()->eta();
	      double phi = et.trackRef()->phi();
	      */

	      Track_ = &tmpTrack;	 
	      innerMomentum_ = TLorentzVector(tmpTrack.innerMomentum().x(), tmpTrack.innerMomentum().y(), tmpTrack.innerMomentum().z(), 0.);
	      outerMomentum_ = TLorentzVector(tmpTrack.outerMomentum().x(), tmpTrack.outerMomentum().y(), tmpTrack.outerMomentum().z(), 0.);
	      Momentum_ = elementTrack.trackRef()->pt();
	      min_dist_cluster = tmp_dist;
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


//----------     --------------------------
float PFCandidateWithFT::GetTrackAlpha(){
  std::cout << "innerMomentum_.Phi() = " << innerMomentum_.Phi() << std::endl;
  std::cout << "outerMomentum_.Phi() = " << outerMomentum_.Phi() << std::endl;
  // float deltaPhi = outerPhi
  // float innerPhi = innerMomentum_.Phi();
  // float outerPhi = innerMomentum_.Phi();

  return innerMomentum_.Phi();
}


