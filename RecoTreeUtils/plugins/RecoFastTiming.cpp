#ifndef FastTiming_RecoTreeUtils_RecoFastTiming_h
#define FastTiming_RecoTreeUtils_RecoFastTiming_h

#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"

#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"
#include "FastTiming/RecoTreeUtils/interface/FTTree.h"

using namespace std;

//****************************************************************************************

class RecoFastTiming : public edm::EDAnalyzer
{
public:
    explicit RecoFastTiming(const edm::ParameterSet&) {};
    ~RecoFastTiming() {};

private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    int iEvent;
    const CaloGeometry* skGeometry;
    //---output file---
    edm::Service<TFileService> fs;
    TFile* outFile;
    FTTree* outTree;   
    //---objects interfaces---
    edm::Handle<vector<SimVertex> > genVtxHandle;
    edm::Handle<vector<reco::PFCandidate> > candHandle;
    // edm::Handle<vector<reco::PFJet> > jetsHandle;
    // edm::Handle<vector<reco::GenJet> > genJetsHandle;
    edm::Handle<edm::SortedCollection<EcalRecHit, 
                                      edm::StrictWeakOrdering<EcalRecHit > > > recSort;    
};

void RecoFastTiming::beginJob()
{
    iEvent=1;
    outFile = &fs->file();
    if(outFile)
    {
        outFile->cd();
        outTree = new FTTree();
    }
}

void RecoFastTiming::endJob()
{
    outFile->cd();
    outTree->Write("fast_timing");
}

void RecoFastTiming::analyze(const edm::Event& Event, const edm::EventSetup& Setup)
{
    if(!outTree)
    {
        cout << "WARNING: output tree is NULL" << endl;
        return;
    }
    outTree->event_n = iEvent;
    iEvent++;            
    //---get gen vertex time---
    const SimVertex* primaryVtx=NULL;
    Event.getByLabel("g4SimHits", genVtxHandle);
    if(genVtxHandle.product()->size() == 0 || genVtxHandle.product()->at(0).vertexId() != 0)
        return;
    primaryVtx = &genVtxHandle.product()->at(0);
    //---fill gen vtx infos 
    outTree->gen_vtx_z = primaryVtx->position().z();
    outTree->gen_vtx_t = primaryVtx->position().t()*1E9;                
    //---get EK detailed time RecHits---
    Event.getByLabel(edm::InputTag("ecalDetailedTimeRecHit", "EcalRecHitsEK", "RECO"),
                      recSort);
    if(!recSort.isValid())
        return;
    vector<EcalRecHit>* recVect = (vector<EcalRecHit>*)recSort.product();
    //---loop over all particles---
    Event.getByLabel("particleFlow", candHandle);
    for(unsigned int iCand=0; iCand<candHandle.product()->size(); iCand++)
    {
        PFCandidateWithFT particle(&candHandle.product()->at(iCand), recVect, primaryVtx);
        if(particle.particleId() > 4 || !particle.GetPFCluster())
            continue;
        outTree->particle_n = iCand;
        outTree->particle_type = particle.particleId();
        outTree->particle_E = particle.energy();
        outTree->particle_pt = particle.pt();
        outTree->particle_eta = particle.eta();
        outTree->particle_phi = particle.phi();
        outTree->maxE_time = particle.GetRecHitTimeMaxE().first;
        outTree->maxE_energy = particle.GetRecHitTimeMaxE().second;                
        outTree->reco_vtx_time = particle.GetTOF();
        outTree->all_time.clear();
        outTree->all_energy.clear();
        outTree->track_length = particle.GetTrackLength();
        outTree->track_radius = particle.GetTrackR();
        particle.GetTrackInfo(outTree->track_alpha, outTree->track_radius,
                              outTree->track_secant, outTree->track_charge);
        outTree->trackCluster_dr = particle.GetDrTrackCluster();                
        vector<pair<float, float> > TandE = particle.GetRecHitsTimeE();
        if(TandE.size() == 0)
            continue;
        for(unsigned int iRec=0; iRec<TandE.size(); iRec++)
        {
            outTree->all_time.push_back(TandE.at(iRec).first);
            outTree->all_energy.push_back(TandE.at(iRec).second);
        }
        outTree->Fill();
    }
}

#endif
//define this as a plug-in
DEFINE_FWK_MODULE(RecoFastTiming);

