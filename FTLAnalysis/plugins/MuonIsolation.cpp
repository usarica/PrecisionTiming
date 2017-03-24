#ifndef _FTL_MUON_ISOLATION_
#define _FTL_MUON_ISOLATION_

// system include files
#include <cmath>
#include <cstdlib>
#include <memory>
#include <unordered_map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "TTree.h"
#include "TRandom.h"

#include "PrecisionTiming/FTLAnalysis/interface/FTLMuonIsoTree.h"

//
// class declaration
//

class FTLMuonIsolation : public edm::EDAnalyzer {
public:
  
    typedef edm::Association<reco::VertexCollection> CandToVertex;
    typedef edm::ValueMap<int> CandToVertexQuality;


    explicit FTLMuonIsolation(const edm::ParameterSet&);
    ~FTLMuonIsolation() {};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override {};
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override {};

    //---member data      
    edm::EDGetTokenT<reco::VertexCollection>            vtxToken_;
    edm::Handle<reco::VertexCollection>                 vtxHandle_;    
    edm::EDGetTokenT<reco::MuonCollection>              muonsToken_;
    edm::Handle<reco::MuonCollection>                   muonsHandle_;    
    edm::EDGetTokenT<edm::View<reco::Track> >           tracksToken_;
    edm::Handle<edm::View<reco::Track> >                tracksHandle_;    
    edm::EDGetTokenT<edm::ValueMap<float> >             timeToken_;
    edm::Handle<edm::ValueMap<float> >                  timeHandle_;    
    edm::EDGetTokenT<edm::ValueMap<float> >             timeResToken_;
    edm::Handle<edm::ValueMap<float> >                  timeResHandle_;        
    // edm::EDGetTokenT<std::vector<std::vector<float> > > ebtimeToken_;
    // edm::Handle<std::vector<std::vector<float> > >      ebtimeHandle_;    
    edm::EDGetTokenT<reco::GenParticleCollection>       genPartToken_;
    edm::Handle<reco::GenParticleCollection>            genPartHandle_;    
    edm::EDGetTokenT<vector<reco::GenJet> >             genJetToken_;
    edm::Handle<vector<reco::GenJet> >                  genJetHandle_;    

    //---I/O
    edm::Service<TFileService> fs;
    FTLMuonIsoTree outTree_;

    //---options
    vector<double> isoConeSizes_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FTLMuonIsolation::FTLMuonIsolation(const edm::ParameterSet& pSet) :
    vtxToken_(consumes<std::vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtxTag"))),    
    muonsToken_(consumes<reco::MuonCollection>(pSet.getUntrackedParameter<edm::InputTag>("muonsTag"))),
    tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("tracksTag"))),
    timeToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("timeTag"))),
    timeResToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("timeResTag"))),
    genPartToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genPartTag"))),
    genJetToken_(consumes<std::vector<reco::GenJet> >(pSet.getUntrackedParameter<edm::InputTag>("genJetsTag")))
{
    outTree_ = FTLMuonIsoTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "Muon tree for FTL studies");
    isoConeSizes_ = pSet.getUntrackedParameter<vector<double> >("isoConeSizes");
}

//
// member functions
//

// ------------ method called for each event  ------------
void
FTLMuonIsolation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //---get input collections
    iEvent.getByToken(muonsToken_, muonsHandle_);
    iEvent.getByToken(tracksToken_, tracksHandle_);
    iEvent.getByToken(timeToken_, timeHandle_);
    iEvent.getByToken(timeResToken_, timeResHandle_);
    iEvent.getByToken(vtxToken_, vtxHandle_);
    iEvent.getByToken(genPartToken_, genPartHandle_);
    iEvent.getByToken(genJetToken_, genJetHandle_);

    //---skip events with no "good" muons
    unsigned int nmuons = 0;
    for (const reco::Muon &muon : *muonsHandle_)
        if (muon.pt() > 5.)
            ++nmuons;
  
    if(!nmuons || vtxHandle_->size()<2)
        return;
          
    //---make a map of vertices to track refs within cuts
    std::unordered_multimap<unsigned,reco::TrackBaseRef>
        vertices_to_tracks_z,
        vertices_to_tracks_zt3,
        vertices_to_tracks_zt4,
        vertices_to_tracks_zt5;
  
    for(unsigned i = 0; i < tracksHandle_->size(); ++i)
    {
        auto ref = tracksHandle_->refAt(i);
        const float time = (*timeHandle_)[ref];
        const float timeReso = (*timeResHandle_)[ref] != 0.f ? (*timeResHandle_)[ref] : 0.170f;
        for(int ivtx = 0; ivtx < (int)vtxHandle_->size(); ++ivtx)
        {
            const auto& thevtx = (*vtxHandle_)[ivtx];
            const float dz = std::abs(ref->dz(thevtx.position()));
            const float dt = std::abs(time - thevtx.t());
            const bool useTime = (thevtx.t() != 0.);

            const float base_cut = std::sqrt(thevtx.tError()*thevtx.tError()
                                             + timeReso*timeReso);

            const float time_cut3 = 3.f*base_cut;
            const float time_cut4 = 4.f*base_cut;
            const float time_cut5 = 5.f*base_cut;

            const bool keepz = ( dz < 0.2f );
            const bool keept3 = (!useTime || std::isnan(dt) || dt < time_cut3);
            const bool keept4 = (!useTime || std::isnan(dt) || dt < time_cut4);
            const bool keept5 = (!useTime || std::isnan(dt) || dt < time_cut5);
      
            if( ref->quality(reco::TrackBase::highPurity) && keepz ) {
                vertices_to_tracks_z.emplace(ivtx, ref);        
                if( keept5 ) {
                    vertices_to_tracks_zt5.emplace(ivtx, ref);
                }
                if( keept4 ) {
                    vertices_to_tracks_zt4.emplace(ivtx, ref);
                }
                if( keept3 ) {
                    vertices_to_tracks_zt3.emplace(ivtx, ref);
                }	
            }      
        }
    }
  
    for(auto &muon : *muonsHandle_)
    {
        //---reset output
        outTree_.Reset();

        //---fill global info
        outTree_.event = iEvent.id().event();
        outTree_.lumi = iEvent.id().luminosityBlock();
        outTree_.run = iEvent.id().run();    
        
        //---basic check
        if(muon.track().isNull() || muon.pt() < 10)
            continue;
    
        int vtx_index = -1;

        //---find the 4D vertex this muon is best associated to..
        float max_weight = 0.f;
        for( unsigned i = 0; i < vtxHandle_->size(); ++i )
        {
            const auto& vtx = (*vtxHandle_)[i];      
            const float weight = vtx.trackWeight(muon.track());
            if( weight > max_weight )
            {
                max_weight = weight;
                vtx_index = i;
            }
        }    
        
        //---use highest ranked if muon doesn't belong to any vertex
        const reco::Vertex& vtx = (vtx_index == -1 ? (*vtxHandle_)[0] : (*vtxHandle_)[vtx_index]);
        const auto tracks_z  = vertices_to_tracks_z.equal_range(vtx_index == -1 ? 0 : vtx_index);
        const auto tracks_zt3 = vertices_to_tracks_zt3.equal_range(vtx_index == -1 ? 0 : vtx_index);
        const auto tracks_zt4 = vertices_to_tracks_zt4.equal_range(vtx_index == -1 ? 0 : vtx_index);
        const auto tracks_zt5 = vertices_to_tracks_zt5.equal_range(vtx_index == -1 ? 0 : vtx_index);
        
        //---fill muon and vtx information
        outTree_.pt = muon.pt();
        outTree_.eta = muon.eta();
        outTree_.phi = muon.phi();
        outTree_.px = muon.px();
        outTree_.py = muon.py();
        outTree_.pz = muon.pz();
        outTree_.muonZ = muon.track()->dz(vtx.position()) + vtx.z();
        outTree_.isLooseMuon = muon::isLooseMuon(muon);
        outTree_.isMediumMuon = muon::isMediumMuon(muon);
        outTree_.isTightMuon = muon::isTightMuon(muon, (*vtxHandle_)[vtx_index==-1 ? 0 : vtx_index]);
        outTree_.nVtxs = vtxHandle_->size();        
        outTree_.vtxIndex = vtx_index;
        outTree_.vtxX = vtx.x();
        outTree_.vtxY = vtx.y();
        outTree_.vtxZ = vtx.z();
        outTree_.vtxT = vtx.t();
        outTree_.vtxIsFake = vtx.isFake();
        outTree_.vtxNdof = vtx.ndof();
        outTree_.vtxChi2 = vtx.chi2();    

        //---compute the varius isolations for all requested cone sizes
        outTree_.nTracksVtx = std::distance(tracks_z.first, tracks_z.second);
        outTree_.chIsoZCut->resize(isoConeSizes_.size());
        outTree_.chIsoZTCut_3sigma->resize(isoConeSizes_.size());
        outTree_.chIsoZTCut_4sigma->resize(isoConeSizes_.size());
        outTree_.chIsoZTCut_5sigma->resize(isoConeSizes_.size());                
        for(int iDR=0; iDR<(int)isoConeSizes_.size(); ++iDR)
        {
            float DR = isoConeSizes_[iDR];
            outTree_.chIsoDR->push_back(DR);
            for( auto it = tracks_z.first; it != tracks_z.second; ++it ) {
                auto ref = it->second.castTo<reco::TrackRef>();
                if( ref == muon.track() ) continue;
                if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                outTree_.chIsoZCut->at(iDR) += ref->pt();
            }
    
            for( auto it = tracks_zt3.first; it != tracks_zt3.second; ++it ) {
                auto ref = it->second.castTo<reco::TrackRef>();
                if( ref == muon.track() ) continue;
                if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                outTree_.chIsoZTCut_3sigma->at(iDR) += ref->pt();
            }

            for( auto it = tracks_zt4.first; it != tracks_zt4.second; ++it ) {
                auto ref = it->second.castTo<reco::TrackRef>();
                if( ref == muon.track() ) continue;
                if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                outTree_.chIsoZTCut_4sigma->at(iDR) += ref->pt();
            }
    
            for( auto it = tracks_zt5.first; it != tracks_zt5.second; ++it ) {
                auto ref = it->second.castTo<reco::TrackRef>();
                if( ref == muon.track() ) continue;
                if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                outTree_.chIsoZTCut_5sigma->at(iDR) += ref->pt();
            }
        }

        //---Gen matching
        outTree_.genMatched = false;
        outTree_.genMatchedPrompt = false;
        outTree_.genMatchedJet = false;
        outTree_.genPt = -99.;
        outTree_.genEta = -99.;
        outTree_.genPhi = -99.;
    
        double mindr = std::numeric_limits<double>::max();
        for (const reco::GenParticle &p : *genPartHandle_)
        {
            if (p.status() != 1) continue;
            if (std::abs(p.pdgId()) != 13) continue;
      
            double dr = reco::deltaR(muon,p);
            if( dr<0.2 && dr<mindr )
            {
                mindr = dr;	
                outTree_.genMatched = true;
                outTree_.genMatchedPrompt = p.isPromptFinalState();
                outTree_.genPt = p.pt();
                outTree_.genEta = p.eta();
                outTree_.genPhi = p.phi();
            }      
        }

        mindr = std::numeric_limits<double>::max();
        for( const auto& jet : *genJetHandle_ ) {
            if( jet.pt() < 15.0 || jet.hadEnergy()/jet.energy() < 0.3)
                continue;
            double dr = reco::deltaR(muon,jet);
            if( dr < 0.3 && dr < mindr )
            {
                outTree_.genMatchedJet = true;
                outTree_.genPt = jet.pt();
                outTree_.genEta = jet.eta();
                outTree_.genPhi = jet.phi();                
            }
        }

        //---Fill tree
        outTree_.GetTTreePtr()->Fill();
    }
}
  
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FTLMuonIsolation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FTLMuonIsolation);

#endif
