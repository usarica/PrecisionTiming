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
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Common/interface/SortedCollection.h"

#include "FastTiming/RecoTreeUtils/interface/VertexWithFT.h"
#include "FastTiming/RecoTreeUtils/interface/FTFile.h"

typedef std::vector<reco::TrackBaseRef >::const_iterator trackRef_iterator;

using namespace std;

//****************************************************************************************

class RecoFastTiming : public edm::EDAnalyzer
{
public:
    explicit RecoFastTiming(const edm::ParameterSet&) {};
    ~RecoFastTiming() {};

    //---utils---
    int          FindPrimaryVtx(PFCandidateWithFT* particle);
    void         AssignParticleToVertices();
    
private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    int iEvent;
    const CaloGeometry* skGeometry;
    const MagneticField* magField;
    //---output file---
    edm::Service<TFileService> fs;
    FTFile* outFile;
    //---objects interfaces---
    edm::ESHandle<MagneticField> magFieldHandle;             
    edm::ESHandle<CaloGeometry> geoHandle;
    edm::Handle<vector<SimVertex> > genVtxHandle;
    edm::Handle<vector<reco::Vertex> > recoVtxHandle;
    edm::Handle<vector<reco::PFCandidate> > candHandle;
    // edm::Handle<vector<reco::PFJet> > jetsHandle;
    // edm::Handle<vector<reco::GenJet> > genJetsHandle;
    edm::Handle<edm::SortedCollection<EcalRecHit, 
                                      edm::StrictWeakOrdering<EcalRecHit > > > recSort;
    vector<VertexWithFT> recoVtxCollection;
    vector<PFCandidateWithFT> particles;
};

void RecoFastTiming::beginJob()
{
    iEvent=1;
    outFile = new FTFile(&fs->file());
}

void RecoFastTiming::endJob()
{
    outFile->cd();
    outFile->particlesTree.Write("ft_particles");
    outFile->verticesTree.Write("ft_vertices");
}

void RecoFastTiming::analyze(const edm::Event& Event, const edm::EventSetup& Setup)
{
    if(!outFile)
    {
        cout << "WARNING: output file is NULL" << endl;
        return;
    }
    outFile->particlesTree.event_n = iEvent;
    ++iEvent;
    //---get the magnetic field---
    Setup.get<IdealMagneticFieldRecord>().get(magFieldHandle);
    magField = magFieldHandle.product();
    //---get the geometry---
    Setup.get<CaloGeometryRecord>().get(geoHandle);
    skGeometry = geoHandle.product();
    //---get gen vertex time---
    const SimVertex* genVtx=NULL;
    Event.getByLabel("g4SimHits", genVtxHandle);
    if(genVtxHandle.product()->size() == 0 || genVtxHandle.product()->at(0).vertexId() != 0)
        return;
    genVtx = &genVtxHandle.product()->at(0);
    //---get reco primary vtxs---
    VertexWithFT* recoVtx=NULL;
    Event.getByLabel("offlinePrimaryVertices", recoVtxHandle);
    recoVtxCollection.clear();
    for(unsigned int iVtx=0; iVtx<recoVtxHandle.product()->size(); ++iVtx)
    {
        if(recoVtxHandle.product()->at(iVtx).isValid() && !recoVtxHandle.product()->at(iVtx).isFake())
            recoVtxCollection.push_back(VertexWithFT(&recoVtxHandle.product()->at(iVtx)));
    }
    if(recoVtxCollection.size() == 0)
        return;
    recoVtx = &recoVtxCollection.at(0);
    //---get EK detailed time RecHits---
    Event.getByLabel(edm::InputTag("ecalDetailedTimeRecHit", "EcalRecHitsEK", "RECO"),
                      recSort);
    if(!recSort.isValid())
        return;
    vector<EcalRecHit>* recVect = (vector<EcalRecHit>*)recSort.product();
    //---convert all particles---
    Event.getByLabel("particleFlow", candHandle);
    PFCandidateWithFT particle;
    particles.clear();
    for(unsigned int iCand=0; iCand<candHandle.product()->size(); ++iCand)
    {
        particle = PFCandidateWithFT(&candHandle.product()->at(iCand), 
                                     recVect, genVtx, recoVtx, skGeometry, magField);        
        if(particle.GetPFCandidate()->particleId() < 4 &&
           particle.GetPFCandidate()->pt() > 0.5 && particle.hasTime())
        {
            FindPrimaryVtx(&particle);
            if(!particle.GetRecoVtx()->hasSeed() ||
               particle.GetRecoVtx()->GetSeedRef()->GetPFCandidate()->pt() < particle.GetPFCandidate()->pt())
            {
                //---store the old seed in di particles collection
                //---redundant since di PFCandidate collection is pt ordered
                if(particle.GetRecoVtx()->hasSeed())
                    particles.push_back(*particle.GetRecoVtx()->GetSeedRef());
                //---set the new seed
                particle.GetRecoVtx()->SetSeed(particle);
            }
            else
                particles.push_back(particle);
        }
        else
            particles.push_back(particle);        
    }
    AssignParticleToVertices();
    //---loop over all particles---
    for(unsigned int iCand=0; iCand<candHandle.product()->size(); ++iCand)
    {
        continue;
        particle = PFCandidateWithFT(&candHandle.product()->at(iCand), 
                                     recVect, genVtx, recoVtx, skGeometry, magField);
        if(particle.GetPFCandidate()->particleId() >= 4)
            continue;        
        //---assign the right primary vtx to the track
        int index_tmp = FindPrimaryVtx(&particle);
        //---fill output tree---
        if(particle.hasTime())
        {
            //---fill gen vtx infos 
            outFile->particlesTree.gen_vtx_z = genVtx->position().z();
            outFile->particlesTree.gen_vtx_t = genVtx->position().t()*1E9;                
            //---particle variables
            outFile->particlesTree.particle_n = iCand;
            outFile->particlesTree.particle_type = particle.GetPFCandidate()->particleId();
            outFile->particlesTree.particle_E = particle.GetPFCandidate()->energy();
            outFile->particlesTree.particle_pt = particle.GetPFCandidate()->pt();
            outFile->particlesTree.particle_eta = particle.GetPFCandidate()->eta();
            outFile->particlesTree.particle_phi = particle.GetPFCandidate()->phi();
            //---ecal time variables
            outFile->particlesTree.maxE_time = particle.GetRecHitTimeMaxE().first;
            outFile->particlesTree.maxE_energy = particle.GetRecHitTimeMaxE().second;
            outFile->particlesTree.all_time.clear();
            outFile->particlesTree.all_energy.clear();
            //---vertex reco info
            outFile->particlesTree.reco_vtx_index = index_tmp;
            outFile->particlesTree.reco_vtx_t = particle.GetRawTime()-particle.GetTOF();
            outFile->particlesTree.reco_vtx_z = particle.GetRecoVtxPos().z();
            //---track info
            outFile->particlesTree.track_length = particle.GetTrackLength();
            outFile->particlesTree.track_length_helix = particle.GetPropagatedTrackLength();
            outFile->particlesTree.track_dz = particle.GetTrack()->dz(particle.GetRecoVtx()->position());
            outFile->particlesTree.track_dxy = particle.GetTrack()->dxy(particle.GetRecoVtx()->position());
            outFile->particlesTree.trackCluster_dr = particle.GetDrTrackCluster();                
            vector<pair<float, float> > TandE = particle.GetRecHitsTimeE();
            if(TandE.size() == 0)
                continue;
            for(unsigned int iRec=0; iRec<TandE.size(); ++iRec)
            {
                outFile->particlesTree.all_time.push_back(TandE.at(iRec).first);
                outFile->particlesTree.all_energy.push_back(TandE.at(iRec).second);
            }
            outFile->particlesTree.Fill();
        }
    }
    sort(recoVtxCollection.begin(), recoVtxCollection.end());
    for(unsigned int iVtx=0; iVtx<recoVtxCollection.size(); ++iVtx)
    {
        outFile->verticesTree.event_n = iEvent;
        outFile->verticesTree.gen_vtx_z = genVtx->position().z();
        outFile->verticesTree.gen_vtx_t = genVtx->position().t()*1E9;
        outFile->verticesTree.reco_vtx_index = iVtx;
        outFile->verticesTree.reco_vtx_sumpt2 = recoVtxCollection.at(iVtx).sumPtSquared();
        outFile->verticesTree.reco_vtx_z = recoVtxCollection.at(iVtx).z();
        outFile->verticesTree.reco_vtx_t =
            recoVtxCollection.at(iVtx).ComputeTime(2, 0.03);
        outFile->verticesTree.reco_vtx_n_part = recoVtxCollection.at(iVtx).GetNPart(2, 0.03);
        outFile->verticesTree.Fill();
    }
}

int RecoFastTiming::FindPrimaryVtx(PFCandidateWithFT* particle)
{
    if(!particle->GetTrack())
        return -1;
    
    int goodVtx=0;
    float dz_min = 100;    
    for(unsigned int iVtx=0; iVtx<recoVtxCollection.size(); ++iVtx)
    {            
        if(!recoVtxCollection.at(iVtx).isFake())
        {                
            float dz_tmp = fabs(particle->GetTrack()->dz(recoVtxCollection.at(iVtx).position()));
            if(dz_tmp < dz_min)
            {
                dz_min = dz_tmp;
                goodVtx = iVtx;
            }
        }
    }
    particle->SetRecoVtx(&recoVtxCollection.at(goodVtx));
    return goodVtx;
}

void RecoFastTiming::AssignParticleToVertices()
{
    vector<PFCandidateWithFT>::iterator it;
    while(particles.size() != 0)
    {
        it=particles.end();
        --it;
        if(!it->GetTrack())
        {
            particles.pop_back();
            continue;
        }
        int goodVtx=-1;
        float dz_min=100, dt_min=10;
        for(unsigned int iVtx=0; iVtx<recoVtxCollection.size(); ++iVtx)
        {
            float dz_tmp = fabs(it->GetTrack()->dz(recoVtxCollection.at(iVtx).position()));
            if(recoVtxCollection.at(iVtx).hasSeed() && dz_tmp < 0.2)
            {
                it->SetRecoVtx(&recoVtxCollection.at(iVtx));
                float dt_tmp = fabs(it->GetVtxTime()-recoVtxCollection.at(iVtx).ComputeTime());
                if(dz_tmp < dz_min && dt_tmp < 0.03*2 && dt_tmp < dt_min)
                    goodVtx=iVtx;
            }
        }
        if(goodVtx != -1)
        {
            recoVtxCollection.at(goodVtx).AddParticle(*it);
            it->SetRecoVtx(&recoVtxCollection.at(goodVtx));
        }
        particles.pop_back();        
    }
}

//define this as a plugin
DEFINE_FWK_MODULE(RecoFastTiming);

