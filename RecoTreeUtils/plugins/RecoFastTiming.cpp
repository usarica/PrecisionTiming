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
    pair<VertexWithFT*, int> FindPrimaryVtx(const reco::Track* track);
    double                   sumPtSquared(const reco::Vertex* v);
    
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
    vector<VertexWithFT*> recoVtxCollection;
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
            recoVtxCollection.push_back(new VertexWithFT(&recoVtxHandle.product()->at(iVtx)));
    }
    if(recoVtxCollection.size() == 0)
        return;
    recoVtx = recoVtxCollection.at(0);
    //---get EK detailed time RecHits---
    Event.getByLabel(edm::InputTag("ecalDetailedTimeRecHit", "EcalRecHitsEK", "RECO"),
                      recSort);
    if(!recSort.isValid())
        return;
    vector<EcalRecHit>* recVect = (vector<EcalRecHit>*)recSort.product();
    //---loop over all particles---
    Event.getByLabel("particleFlow", candHandle);
    particles.clear();
    for(unsigned int iCand=0; iCand<candHandle.product()->size(); ++iCand)
    {
        PFCandidateWithFT particle(&candHandle.product()->at(iCand), 
                                   recVect, genVtx, recoVtx, skGeometry, magField);
        if(particle.particleId() >= 4 || !particle.GetPFCluster())
            continue;
        //---assign the right primary vtx to the track
        pair<VertexWithFT*, int> vertex_info;
        vertex_info = FindPrimaryVtx(particle.GetTrack());
        particle.SetRecoVtx(vertex_info.first);
        //---fill output tree---
        //---fill gen vtx infos 
        outFile->particlesTree.gen_vtx_z = genVtx->position().z();
        outFile->particlesTree.gen_vtx_t = genVtx->position().t()*1E9;                
        //---particle variables
        outFile->particlesTree.particle_n = iCand;
        outFile->particlesTree.particle_type = particle.particleId();
        outFile->particlesTree.particle_E = particle.energy();
        outFile->particlesTree.particle_pt = particle.pt();
        outFile->particlesTree.particle_eta = particle.eta();
        outFile->particlesTree.particle_phi = particle.phi();
        //---ecal time variables
        outFile->particlesTree.maxE_time = particle.GetRecHitTimeMaxE().first;
        outFile->particlesTree.maxE_energy = particle.GetRecHitTimeMaxE().second;
        outFile->particlesTree.all_time.clear();
        outFile->particlesTree.all_energy.clear();
        //---vertex reco info
        outFile->particlesTree.reco_vtx_index = vertex_info.second;
        outFile->particlesTree.reco_vtx_t = particle.GetRawTime()-particle.GetTOF();
        outFile->particlesTree.reco_vtx_z = particle.GetRecoVtxPos().z();
        //---track info
        outFile->particlesTree.track_length = particle.GetTrackLength();
        outFile->particlesTree.track_length_helix = particle.GetPropagatedTrackLength();
        outFile->particlesTree.track_dz = particle.GetTrack()->dz(particle.GetRecoVtx()->position());
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
        //---add the new particle to the right vtx---
        particles.push_back(particle);
        vertex_info.first->AddParticle(&particles.at(particles.size()-1));
    }
    for(unsigned int iVtx=0; iVtx<recoVtxCollection.size(); ++iVtx)
    {
        outFile->verticesTree.event_n = iEvent;
        outFile->verticesTree.gen_vtx_z = genVtx->position().z();
        outFile->verticesTree.gen_vtx_t = genVtx->position().t()*1E9;
        outFile->verticesTree.reco_vtx_index = iVtx;
        outFile->verticesTree.reco_vtx_n_part = recoVtxCollection.at(iVtx)->GetNPart();
        outFile->verticesTree.reco_vtx_sumpt2 =
            sumPtSquared(recoVtxCollection.at(iVtx)->GetRecoVtxRef());
        outFile->verticesTree.reco_vtx_z = recoVtxCollection.at(iVtx)->z();
        outFile->verticesTree.reco_vtx_t =
            recoVtxCollection.at(iVtx)->ComputeWightedTime(0.2, 0.03);        
        outFile->verticesTree.Fill();
    }
}

pair<VertexWithFT*, int> RecoFastTiming::FindPrimaryVtx(const reco::Track* track)
{
    if(!track)
    {
        VertexWithFT* fake=NULL;
        return make_pair(fake, -1);
    }
    int goodVtx=0;
    float dz_min = 100;    
    for(unsigned int iVtx=0; iVtx<recoVtxCollection.size(); ++iVtx)
    {            
        if(!recoVtxCollection.at(iVtx)->isFake())
        {                
            float dz_tmp = fabs(track->dz(recoVtxCollection.at(iVtx)->position()));
            if(dz_tmp < dz_min)
            {
                dz_min = dz_tmp;
                goodVtx = iVtx;
            }
        }
    }
    return make_pair(recoVtxCollection.at(goodVtx), goodVtx);
}

double RecoFastTiming::sumPtSquared(const reco::Vertex* v) 
{
    double sum = 0.;
    double pT;
    for (reco::Vertex::trackRef_iterator it = v->tracks_begin(); it != v->tracks_end(); ++it)
    {
        pT = (**it).pt();
        double epT=(**it).ptError();
        pT=pT>epT ? pT-epT : 0;
        
        sum += pT*pT;
    }
    return sum;
}

//define this as a plugin
DEFINE_FWK_MODULE(RecoFastTiming);

