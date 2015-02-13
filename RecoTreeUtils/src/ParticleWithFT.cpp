#include "FastTiming/RecoTreeUtils/interface/ParticleWithFT.hpp"

//**********Ctors*************************************************************************

ParticleWithFT::ParticleWithFT():
    clusterE_(0), maxRecHitE_(0), time_(0)
{}

ParticleWithFT::ParticleWithFT(PFCandidate* PFCand, vector<EcalRecHit>* ecalRecHits):
    clusterE_(0), maxRecHitE_(0), time_(0), pfCand_(PFCand)
{
    float min_dist_cluster=100;
    for(auto& blockPair : pfCand_.elementsInBlocks())
    {        
	unsigned int pos = blockPair.second;
	const reco::PFBlockElement& blockElement = blockPair.first->elements()[pos];
	//{
    	// if(block.elements()[iEle].type() != 4)
    	//     continue;
        // sub_det_elements clusters = GetLinks(&block, iEle);        
        // for(unsigned int iCls=0; iCls<clusters.size(); iCls++)
        // {
        //     if(clusters.at(iCls).second->type() != 4)
        //         continue;
            if(blockElement.type() != 4)
                continue;
            reco::PFCluster tmpCluster = *blockElement.clusterRef().get();	 
            //reco::PFCluster tmpCluster = *clusters.at(iCls).second->clusterRef().get();
            tmpCluster.calculatePositionREP();
            REPPoint pfClusterPos = tmpCluster.positionREP();
            float tmp_dist=DeltaR(pfClusterPos.Eta(), photon->eta(),
                                  pfClusterPos.Phi(), photon->phi());
            if(tmp_dist < min_dist_cluster && blockElement.clusterRef().isAvailable())
            {
                min_dist_cluster = tmp_dist;
                pfCluster = blockElement.clusterRef().get();//clusters.at(iCls).second->clusterRef().get();                
                cout << "DeltaR        : " << min_dist_cluster << endl
                     << "Cluster Energy: " << pfCluster->energy() << endl;
            }
            //}        
    }
    if(pfCluster)
    {
        if(clusterList->find(pfCluster) == clusterList->end() ||
           clusterList->at(pfCluster).second > min_dist_cluster)
            (*clusterList)[pfCluster] = pair<const reco::PFCandidate*, float>(photon, min_dist_cluster);
        return 0;
    }
    return -1;
