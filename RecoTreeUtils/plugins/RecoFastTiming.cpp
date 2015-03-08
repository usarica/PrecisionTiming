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
    explicit RecoFastTiming(const edm::ParameterSet&);
    ~RecoFastTiming() {};

    //---utils---
    int          FindVtxSeedParticle(PFCandidateWithFT* particle);
    void         AssignChargedToVtxs(vector<PFCandidateWithFT*>* charged_cand);
    void         AssignNeutralToVtxs(vector<PFCandidateWithFT*>* neutral_cand);
    void         ProcessVertices(const SimVertex* genVtx);
    void         ProcessParticles(vector<PFCandidateWithFT*>* particles, int iVtx);
        
private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    int iEvent_;
    const CaloGeometry* skGeometry_;
    const MagneticField* magField_;
    //---output file---
    edm::Service<TFileService> fs_;
    FTFile* outFile_;
    //---objects interfaces---
    edm::ESHandle<MagneticField> magFieldHandle_;             
    edm::ESHandle<CaloGeometry> geoHandle_;
    edm::Handle<vector<SimVertex> > genVtxHandle_;
    edm::Handle<vector<reco::Vertex> > recoVtxHandle_;
    edm::Handle<vector<reco::PFCandidate> > candHandle_;
    // edm::Handle<vector<reco::PFJet> > jetsHandle;
    // edm::Handle<vector<reco::GenJet> > genJetsHandle;
    edm::Handle<edm::SortedCollection<EcalRecHit, 
                                      edm::StrictWeakOrdering<EcalRecHit > > > recSort_;
    //---FT objects---
    vector<VertexWithFT> recoVtxCollection_;
    vector<PFCandidateWithFT> particlesCollection_;
    //---Options---
    float tRes_;
    float dzCut_;       
    float ptCut_;       
    bool saveParticles_;
    bool saveAllRecHits_;
    bool saveVertices_;
};

RecoFastTiming::RecoFastTiming(const edm::ParameterSet& Config)
{
    tRes_ = Config.getUntrackedParameter<double>("timeResSmearing");
    dzCut_ = Config.getUntrackedParameter<double>("dzCut");
    ptCut_ = Config.getUntrackedParameter<double>("ptCut");
    saveParticles_ = Config.getUntrackedParameter<bool>("saveParticles");
    saveAllRecHits_ = Config.getUntrackedParameter<bool>("saveAllRecHits");
    saveVertices_ = Config.getUntrackedParameter<bool>("saveVertices");
}

void RecoFastTiming::beginJob()
{
    iEvent_=1;
    outFile_ = new FTFile(&fs_->file());
}

void RecoFastTiming::endJob()
{
    outFile_->cd();
    outFile_->particlesTree.Write("ft_particles");
    outFile_->verticesTree.Write("ft_vertices");
}

void RecoFastTiming::analyze(const edm::Event& Event, const edm::EventSetup& Setup)
{
    if(!outFile_)
    {
        cout << "WARNING: output file is NULL" << endl;
        return;
    }
    outFile_->particlesTree.event_n = iEvent_;
    ++iEvent_;
    //---get the magnetic field---
    Setup.get<IdealMagneticFieldRecord>().get(magFieldHandle_);
    magField_ = magFieldHandle_.product();
    //---get the geometry---
    Setup.get<CaloGeometryRecord>().get(geoHandle_);
    skGeometry_ = geoHandle_.product();
    //---get gen vertex time---
    const SimVertex* genVtx=NULL;
    Event.getByLabel("g4SimHits", genVtxHandle_);
    if(genVtxHandle_.product()->size() == 0 || genVtxHandle_.product()->at(0).vertexId() != 0)
        return;
    genVtx = &genVtxHandle_.product()->at(0);
    //---get reco primary vtxs---
    VertexWithFT* recoVtx=NULL;
    Event.getByLabel("offlinePrimaryVertices", recoVtxHandle_);
    recoVtxCollection_.clear();
    for(unsigned int iVtx=0; iVtx<recoVtxHandle_.product()->size(); ++iVtx)
    {
        if(recoVtxHandle_.product()->at(iVtx).isValid() &&
           !recoVtxHandle_.product()->at(iVtx).isFake() &&
           recoVtxHandle_.product()->at(iVtx).chi2()/recoVtxHandle_.product()->at(iVtx).ndof() < 2 &&
           recoVtxHandle_.product()->at(iVtx).chi2()/recoVtxHandle_.product()->at(iVtx).ndof() > 0.5)
        {
            recoVtxCollection_.push_back(VertexWithFT(&recoVtxHandle_.product()->at(iVtx)));
        }
    }
    if(recoVtxCollection_.size() == 0)
        return;
    recoVtx = &recoVtxCollection_.at(0);
    //---get EK detailed time RecHits---
    Event.getByLabel(edm::InputTag("ecalDetailedTimeRecHit", "EcalRecHitsEK", "RECO"),
                      recSort_);
    if(!recSort_.isValid())
        return;
    vector<EcalRecHit>* recVect = (vector<EcalRecHit>*)recSort_.product();
    //---convert all particles---
    Event.getByLabel("particleFlow", candHandle_);
    PFCandidateWithFT particle;
    particlesCollection_.clear();
    for(unsigned int iCand=0; iCand<candHandle_.product()->size(); ++iCand)
    {
        particle = PFCandidateWithFT(&candHandle_.product()->at(iCand), 
                                     recVect, genVtx, recoVtx, skGeometry_, magField_);
        particlesCollection_.push_back(particle);
    }
    //---search for seeds particles---
    vector<PFCandidateWithFT*> chargedRefs;
    vector<PFCandidateWithFT*> neutralRefs;
    for(unsigned int iPart=0; iPart<particlesCollection_.size(); ++iPart)
    {
        PFCandidateWithFT* particleRef = &particlesCollection_.at(iPart);
        if(particleRef->particleId() < 4)
        {
            FindVtxSeedParticle(particleRef);
            if(particleRef->hasTime() && particleRef->GetRecoVtx() &&
               (!particleRef->GetRecoVtx()->hasSeed() ||
                particleRef->GetRecoVtx()->GetSeedRef()->pt() < particleRef->pt()))
            {
                //---store the old seed in di particles collection
                //---redundant since di particles collection is pt ordered
                if(particleRef->GetRecoVtx()->hasSeed())
                    chargedRefs.push_back(particleRef->GetRecoVtx()->GetSeedRef());
                //---set the new seed
                particleRef->GetRecoVtx()->SetSeed(particleRef);
            }
            else
                chargedRefs.push_back(particleRef);
        }
        else 
            neutralRefs.push_back(particleRef);        
    }
    //---link the remaining charged particles to vtxs---
    AssignChargedToVtxs(&chargedRefs);
    //---link the neutral particles to vtxs---
    //AssignNeutralToVtxs(&neutralRefs);
    //---sort the vtxs by sumpt2---
    sort(recoVtxCollection_.begin(), recoVtxCollection_.end());    
    //---compute vtxs times and fill the output tree---
    if(saveVertices_)
        ProcessVertices(genVtx);
    else if(saveParticles_)
    {
        vector<PFCandidateWithFT*> allParticlesRefs;
        for(unsigned int iPart=0; iPart<particlesCollection_.size(); ++iPart)
            allParticlesRefs.push_back(&particlesCollection_.at(iPart));

        ProcessParticles(&allParticlesRefs, -1);
    }

    return;
}

int RecoFastTiming::FindVtxSeedParticle(PFCandidateWithFT* particle)
{
    if(!particle->GetTrack())
        return -1;
    
    int goodVtx=-1;
    float dz_min = 100;    
    for(unsigned int iVtx=0; iVtx<recoVtxCollection_.size(); ++iVtx)
    {            
        float dz_tmp = fabs(particle->GetTrack()->dz(recoVtxCollection_.at(iVtx).position()));
        if(dz_tmp < dzCut_ && dz_tmp < dz_min)
        {
            dz_min = dz_tmp;
            goodVtx = iVtx;
        }
    }
    if(goodVtx != -1)
        particle->SetRecoVtx(&recoVtxCollection_.at(goodVtx));
    
    return goodVtx;
}

void RecoFastTiming::AssignChargedToVtxs(vector<PFCandidateWithFT*>* charged_cand)
{
    vector<PFCandidateWithFT*>::iterator it;
    while(charged_cand->size() != 0)
    {
        it=charged_cand->end();
        --it;
        if(!(*it)->GetTrack())
        {
            charged_cand->pop_back();
            continue;
        }
        int goodVtx=-1;
        float dz_min=100, dt_min=10;
        for(unsigned int iVtx=0; iVtx<recoVtxCollection_.size(); ++iVtx)
        {
            float dz_tmp = fabs((*it)->GetTrack()->dz(recoVtxCollection_.at(iVtx).position()));
            if(recoVtxCollection_.at(iVtx).hasSeed() && dz_tmp < dzCut_)
            {
                (*it)->SetRecoVtx(&recoVtxCollection_.at(iVtx));
                float dt_tmp=0;
                dt_tmp = fabs((*it)->GetVtxTime() - recoVtxCollection_.at(iVtx).ComputeTime(ptCut_, tRes_));                
                if(dz_tmp < dz_min && dt_tmp < tRes_*2 && dt_tmp < dt_min)
                {
                    goodVtx=iVtx;
                    // dz_min=dz_tmp;
                    dt_min=dt_tmp;
                }
            }
        }
        if(goodVtx != -1)
        {
            recoVtxCollection_.at(goodVtx).AddParticle(*it);
            (*it)->SetRecoVtx(&recoVtxCollection_.at(goodVtx));
        }
        charged_cand->pop_back();        
    }
}

void RecoFastTiming::AssignNeutralToVtxs(vector<PFCandidateWithFT*>* neutral_cand)
{
    vector<PFCandidateWithFT*>::iterator it;
    while(neutral_cand->size() != 0)
    {
        it=neutral_cand->end();
        --it;
        if(!(*it)->hasTime())
        {
            neutral_cand->pop_back();
            continue;
        }
        int goodVtx=-1;
        float dt_min=10;
        for(unsigned int iVtx=0; iVtx<recoVtxCollection_.size(); ++iVtx)
        {
            float dz_tmp = fabs((*it)->GetTrack()->dz(recoVtxCollection_.at(iVtx).position()));
            if(recoVtxCollection_.at(iVtx).hasSeed() && dz_tmp < dzCut_)
            {
                (*it)->SetRecoVtx(&recoVtxCollection_.at(iVtx));
                float dt_tmp=0;
                //if(recoVtxCollection_.at(iVtx).hasSeed())
                dt_tmp = fabs((*it)->GetVtxTime() - recoVtxCollection_.at(iVtx).ComputeTime(ptCut_, tRes_));                
                if(dt_tmp < tRes_*2 && dt_tmp < dt_min)
                {                    
                    goodVtx=iVtx;
                    dt_min=dt_tmp;
                }
            }
        }
        if(goodVtx != -1)
        {
            recoVtxCollection_.at(goodVtx).AddParticle(*it);
            (*it)->SetRecoVtx(&recoVtxCollection_.at(goodVtx));
        }
        neutral_cand->pop_back();        
    }
}



void RecoFastTiming::ProcessVertices(const SimVertex* genVtx)
{
    for(unsigned int iVtx=0; iVtx<recoVtxCollection_.size(); ++iVtx)
    {
        //---fix references after sort
        recoVtxCollection_.at(iVtx).FixVtxRefs();
        outFile_->verticesTree.event_n = iEvent_;
        outFile_->verticesTree.gen_vtx_z = genVtx->position().z();
        outFile_->verticesTree.gen_vtx_t = genVtx->position().t()*1E9;
        outFile_->verticesTree.reco_vtx_index = iVtx;
        outFile_->verticesTree.reco_vtx_ndof = recoVtxCollection_.at(iVtx).GetRecoVtxRef()->ndof();
        outFile_->verticesTree.reco_vtx_chi2 = recoVtxCollection_.at(iVtx).GetRecoVtxRef()->chi2();
        outFile_->verticesTree.reco_vtx_sumpt2 = recoVtxCollection_.at(iVtx).sumPtSquared();
        outFile_->verticesTree.reco_vtx_seed_pt = recoVtxCollection_.at(iVtx).GetSeedRef() ?
            recoVtxCollection_.at(iVtx).GetSeedRef()->pt() : -1;
        outFile_->verticesTree.reco_vtx_seed_t = recoVtxCollection_.at(iVtx).GetSeedRef() ?
            recoVtxCollection_.at(iVtx).GetSeedRef()->GetVtxTime(tRes_) : -100;
        outFile_->verticesTree.reco_vtx_z = recoVtxCollection_.at(iVtx).z();
        outFile_->verticesTree.reco_vtx_t =
            recoVtxCollection_.at(iVtx).ComputeTime(ptCut_, tRes_);
        outFile_->verticesTree.reco_vtx_n_part = recoVtxCollection_.at(iVtx).GetNPart(ptCut_, tRes_);
        outFile_->verticesTree.Fill();
        if(saveParticles_)
        {
            vector<PFCandidateWithFT*> vtxParticles = recoVtxCollection_.at(iVtx).GetParticles();
            ProcessParticles(&vtxParticles, iVtx);
        }
    }
}

void RecoFastTiming::ProcessParticles(vector<PFCandidateWithFT*>* particles, int iVtx)
{
    int index=-1;
    //---loop over all particles---
    for(vector<PFCandidateWithFT*>::iterator it=particles->begin(); it!=particles->end(); ++it)
    {
        ++index;
        //---fill output tree---
        if((*it)->hasTime())
        {            
            //---fill gen vtx infos 
            // outFile_->particlesTree.gen_vtx_z = genVtx->position().z();
            // outFile_->particlesTree.gen_vtx_t = genVtx->position().t()*1E9;                
            //---particle variables
            outFile_->particlesTree.particle_n = index;
            outFile_->particlesTree.particle_type = (*it)->particleId();
            outFile_->particlesTree.particle_E = (*it)->energy();
            outFile_->particlesTree.particle_pt = (*it)->pt();
            outFile_->particlesTree.particle_eta = (*it)->eta();
            outFile_->particlesTree.particle_phi = (*it)->phi();
            //---ecal time variables
            outFile_->particlesTree.maxE_time = (*it)->GetRecHitTimeMaxE().first;
            outFile_->particlesTree.maxE_energy = (*it)->GetRecHitTimeMaxE().second;
            outFile_->particlesTree.all_time.clear();
            outFile_->particlesTree.all_energy.clear();
            //---vertex reco info
            outFile_->particlesTree.reco_vtx_index = iVtx;
            outFile_->particlesTree.reco_vtx_t = (*it)->GetRawTime()-(*it)->GetTOF();
            outFile_->particlesTree.reco_vtx_z = (*it)->GetRecoVtxPos().z();
            //---track info
            outFile_->particlesTree.track_length = (*it)->GetTrackLength();
            outFile_->particlesTree.track_length_helix = (*it)->GetPropagatedTrackLength();
            outFile_->particlesTree.track_dz = (*it)->GetTrack()->dz((*it)->GetRecoVtx()->position());
            outFile_->particlesTree.track_dxy = (*it)->GetTrack()->dxy((*it)->GetRecoVtx()->position());
            outFile_->particlesTree.trackCluster_dr = (*it)->GetDrTrackCluster();
            vector<pair<float, float> > TandE;
            if(saveAllRecHits_)
                TandE = (*it)->GetRecHitsTimeE();
            else
                TandE.push_back(make_pair(0, 0));
            if(TandE.size() != 0)
            {
                for(unsigned int iRec=0; iRec<TandE.size(); ++iRec)
                {
                    outFile_->particlesTree.all_time.push_back(TandE.at(iRec).first);
                    outFile_->particlesTree.all_energy.push_back(TandE.at(iRec).second);
                }
            }
            outFile_->particlesTree.Fill();
        }
    }
}

//define this as a plugin
DEFINE_FWK_MODULE(RecoFastTiming);

