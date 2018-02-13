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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
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
  
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;
    
    explicit FTLMuonIsolation(const edm::ParameterSet&);
    ~FTLMuonIsolation() {};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override {};
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override {};

    //---member data
    edm::EDGetTokenT<genXYZ>                             genXYZToken_;
    edm::Handle<genXYZ>                                  genXYZHandle_;
    edm::EDGetTokenT<float>                              genT0Token_;    
    edm::Handle<float>                                   genT0Handle_;
    edm::EDGetTokenT<vector<SimVertex> >                 genVtxToken_;
    edm::Handle<vector<SimVertex> >                      genVtxHandle_;    
    edm::EDGetTokenT<reco::VertexCollection>             vtxToken_;
    edm::Handle<reco::VertexCollection>                  vtxHandle_;    
    edm::EDGetTokenT<reco::MuonCollection>               muonsToken_;
    edm::Handle<reco::MuonCollection>                    muonsHandle_;    
    edm::EDGetTokenT<edm::View<reco::PFCandidate> >      tracksToken_;
    edm::Handle<edm::View<reco::PFCandidate> >           tracksHandle_;    
    edm::EDGetTokenT<edm::View<reco::GenParticle> >      genPartToken_;
    edm::Handle<edm::View<reco::GenParticle> >           genPartHandle_;    
    edm::EDGetTokenT<vector<reco::GenJet> >              genJetToken_;
    edm::Handle<vector<reco::GenJet> >                   genJetHandle_;    

    //---
    vector<double> targetResolutions_;
    
    //---I/O
    int iEvent_;
    edm::Service<TFileService> fs;
    map<double, FTLMuonIsoTree> outTrees_;

    //---options    
    bool           useMCTruthPV_;
    bool           isTimingSample_;
    bool           saveTracksInfo_;
    float          dzCut_;
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
    genXYZToken_(consumes<genXYZ>(pSet.getUntrackedParameter<edm::InputTag>("genXYZTag"))),    
    genT0Token_(consumes<float>(pSet.getUntrackedParameter<edm::InputTag>("genT0Tag"))),        
    genVtxToken_(consumes<vector<SimVertex> >(pSet.getUntrackedParameter<edm::InputTag>("genVtxTag"))),    
    vtxToken_(consumes<std::vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtxTag"))),    
    muonsToken_(consumes<vector<reco::Muon> >(pSet.getUntrackedParameter<edm::InputTag>("muonsTag"))),
    tracksToken_(consumes<edm::View<reco::PFCandidate> >(pSet.getUntrackedParameter<edm::InputTag>("tracksTag"))),
    genPartToken_(consumes<edm::View<reco::GenParticle> >(pSet.getUntrackedParameter<edm::InputTag>("genPartTag"))),
    genJetToken_(consumes<std::vector<reco::GenJet> >(pSet.getUntrackedParameter<edm::InputTag>("genJetsTag"))),
    targetResolutions_(pSet.getUntrackedParameter<vector<double> >("targetResolutions"))
{
    iEvent_ = 0;
    useMCTruthPV_ = pSet.getUntrackedParameter<bool>("useMCTruthPV");
    isTimingSample_ = pSet.getUntrackedParameter<bool>("isTimingSample");
    saveTracksInfo_ = pSet.getUntrackedParameter<bool>("saveTracksInfo");
    dzCut_ = pSet.getUntrackedParameter<double>("dzCut");
    isoConeSizes_ = pSet.getUntrackedParameter<vector<double> >("isoConeSizes");
    if(isTimingSample_)
    {
        for(auto& res : targetResolutions_)
            outTrees_[res] = FTLMuonIsoTree((pSet.getUntrackedParameter<string>("treeName")+"_"+to_string(int(res*1000))).c_str(),
                                            "Muon tree for FTL studies");
    }
    else
    {
        targetResolutions_.clear();
        targetResolutions_.push_back(0.03);
        outTrees_[0.03] = FTLMuonIsoTree((pSet.getUntrackedParameter<string>("treeName")+"_notiming").c_str(),
                                         "Muon tree for FTL studies");
    }
}

//
// member functions
//

// ------------ method called for each event  ------------
void
FTLMuonIsolation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //---get input collections
    iEvent.getByToken(genXYZToken_, genXYZHandle_);
    iEvent.getByToken(genT0Token_, genT0Handle_);
    iEvent.getByToken(muonsToken_, muonsHandle_);
    iEvent.getByToken(tracksToken_, tracksHandle_);
    if(useMCTruthPV_)
        iEvent.getByToken(genVtxToken_, genVtxHandle_);    
    iEvent.getByToken(vtxToken_, vtxHandle_);
    iEvent.getByToken(genPartToken_, genPartHandle_);
    iEvent.getByToken(genJetToken_, genJetHandle_);

    //---skip bad events (checks on muons and MC-truth vtx)
    unsigned int nmuons = 0;
    for (const auto &muon : *muonsHandle_)
        if (muon.pt() > 5.)
            ++nmuons;

    ++iEvent_;
    
    //---get truth PV
    SimVertex genPV;
    if(genVtxHandle_.isValid())
        genPV = genVtxHandle_.product()->at(0);
    else
    {
        auto xyz = genXYZHandle_.product();
        auto t = *genT0Handle_.product();
        auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
        genPV = SimVertex(v, t);
    }

    //---make a map of vertices to track refs within cuts
    std::unordered_multimap<unsigned, reco::TrackRef >
        vertices_cleaned_tracks_zt3,
        vertices_to_tracks_z,        
        vertices_to_tracks_z1,
        vertices_to_tracks_z3,
        vertices_to_tracks_z5,
        vertices_to_tracks_z7,
        vertices_to_tracks_z10,        
        vertices_to_tracks_zt3,
        vertices_to_tracks_zt4,
        vertices_to_tracks_zt5,
        vertices_to_tracks_zt7,
        vertices_to_tracks_zt10;

    for(auto& iRes : targetResolutions_)
    {
        vertices_to_tracks_z.clear();        
        vertices_to_tracks_z1.clear();
        vertices_to_tracks_z3.clear();
        vertices_to_tracks_z5.clear();
        vertices_to_tracks_z7.clear();
        vertices_to_tracks_z10.clear();        
        vertices_to_tracks_zt3.clear();
        vertices_to_tracks_zt4.clear();
        vertices_to_tracks_zt5.clear();
        vertices_to_tracks_zt7.clear();
        vertices_to_tracks_zt10.clear();

        vector<float> tracks_pt, tracks_eta, tracks_phi, tracks_dz, tracks_dxy, tracks_t;
        vector<bool> tracks_keepz,tracks_keept;
        
        float extra_smearing = std::sqrt(iRes*iRes - 0.03*0.03);
        for(unsigned i = 0; i < tracksHandle_->size(); ++i)
        {
            auto ref = tracksHandle_->refAt(i);

            //---skip neutrals
            if(ref->charge() == 0)
                continue;
            reco::TrackRef track = ref->trackRef();
            if(track.isNull() || !track->quality(reco::TrackBase::highPurity))
                continue;
            
            if(isTimingSample_)
            {
                float time = ref->time();
                float timeReso = ref->timeError() != -1. ? ref->timeError() : 1.;
                time *= gRandom->Gaus(1., extra_smearing);
                timeReso = std::sqrt(timeReso*timeReso + extra_smearing*extra_smearing);
                for(int ivtx = 0; ivtx < (int)vtxHandle_->size(); ++ivtx)
                {
                    const auto& thevtx = (*vtxHandle_)[ivtx];
                    const float dz = std::abs(track->dz(thevtx.position()));
                    const float dxy = std::abs(track->dxy(thevtx.position()));
                    const float dt = std::abs(time - thevtx.t());
                    // const float dz2 = std::pow(ref->dz(thevtx.position()), 2);
                    // const float dt2 = std::pow(time - thevtx.t(), 2);
                    const bool useTime = (thevtx.t() != 0.);

                    const float base_cut = std::sqrt(thevtx.tError()*thevtx.tError()
                                                     + timeReso*timeReso);
            
                    const float time_cut3 = 3.f*base_cut;
                    const float time_cut4 = 4.f*base_cut;
                    const float time_cut5 = 5.f*base_cut;
                    const float time_cut7 = 7.f*base_cut;
                    const float time_cut10 = 20.f*base_cut;

                    const bool keepz = ( dz < dzCut_ && dxy < 0.02 );
                    const bool keept3 = (!useTime || std::isnan(dt) || dt < time_cut3);
                    const bool keept4 = (!useTime || std::isnan(dt) || dt < time_cut4);
                    const bool keept5 = (!useTime || std::isnan(dt) || dt < time_cut5);
                    const bool keept7 = (!useTime || std::isnan(dt) || dt < time_cut7);
                    const bool keept10 = (!useTime || std::isnan(dt) || dt < time_cut10);
            
                    if( keepz ) {
                        if( keept10 ) {
                            vertices_to_tracks_zt10.emplace(ivtx, track);
                        }
                        if( keept7 ) {
                            vertices_to_tracks_zt7.emplace(ivtx, track);
                        }
                        if( keept5 ) {
                            vertices_to_tracks_zt5.emplace(ivtx, track);
                        }
                        if( keept4 ) {
                            vertices_to_tracks_zt4.emplace(ivtx, track);
                        }
                        if( keept3 ) {
                            vertices_to_tracks_zt3.emplace(ivtx, track);
                        }
                        else
                            vertices_cleaned_tracks_zt3.emplace(ivtx, track);
                    }

                    if(saveTracksInfo_)
                    {
                        tracks_pt.push_back(track->pt());
                        tracks_eta.push_back(track->eta());
                        tracks_phi.push_back(track->phi());
                        tracks_dz.push_back(dz);
                        tracks_dxy.push_back(dxy);                            
                        tracks_t.push_back(time);
                        tracks_keepz.push_back(keepz);
                        tracks_keept.push_back(keept3);
                    }
                }
            }

            //---no timing
            for(int ivtx = 0; ivtx < (int)vtxHandle_->size(); ++ivtx)
            {
                const auto& thevtx = (*vtxHandle_)[ivtx];
                const float dz = std::abs(track->dz(thevtx.position()));
                const float dxy = std::abs(track->dxy(thevtx.position()));                
                const float rel_dz = std::abs(dz/std::sqrt(track->dzError()*track->dzError()+
                                                           thevtx.zError()*thevtx.zError()));
                const bool keepz = ( dz < dzCut_ && dxy < 0.02  );
                const bool keepz1 = ( rel_dz < 1 );
                const bool keepz3 = ( rel_dz < 3 );
                const bool keepz5 = ( rel_dz < 5 );
                const bool keepz7 = ( rel_dz < 7 );
                const bool keepz10 = ( rel_dz < 10 );

                if(keepz)
                    vertices_to_tracks_z.emplace(ivtx, track);
                if(!isTimingSample_)
                {
                    if(keepz1)
                        vertices_to_tracks_z1.emplace(ivtx, track);
                    if(keepz3)
                        vertices_to_tracks_z3.emplace(ivtx, track);
                    if(keepz5)
                        vertices_to_tracks_z5.emplace(ivtx, track);
                    if(keepz7)
                        vertices_to_tracks_z7.emplace(ivtx, track);
                    if(keepz10)
                        vertices_to_tracks_z10.emplace(ivtx, track);
                }
            }
        }

        for(unsigned iMuon=0; iMuon<muonsHandle_->size(); ++iMuon)
        {
            auto muon = muonsHandle_->at(iMuon);
            
            //---reset output
            outTrees_[iRes].Reset();

            //---fill global info            
            outTrees_[iRes].event = iEvent.id().event();
            outTrees_[iRes].lumi = iEvent.id().luminosityBlock();
            outTrees_[iRes].run = iEvent.id().run();    
        
            //---basic check
            if(muon.track().isNull() || muon.pt() < 10)
                continue;
    
            int vtx_index = -1;
            
            //---find the vertex this muon is best associated to..
            if(useMCTruthPV_)
            {
                double min_dz = std::numeric_limits<double>::max();
                double min_dzdt = std::numeric_limits<double>::max();
                for( unsigned i = 0; i < vtxHandle_->size(); ++i )
                {
                    const auto& vtx = (*vtxHandle_)[i];
                    const float dz = std::abs(vtx.z() - genPV.position().z());
                    if( dz < min_dz )
                    {
                        min_dz = dz;
                        vtx_index = i;
                    }
                    // if( dz < 0.1 )
                    // {
                    //     if(isTimingSample_)
                    //     {
                    //         const double dzdt = pow((vtx.z() - genPV.position().z())/vtx.zError(), 2) +
                    //             pow((vtx.t()-genPV.position().t())/vtx.tError(), 2);                            
                    //         if( dzdt < min_dzdt )
                    //         {
                    //             min_dzdt = dzdt;
                    //             vtx_index = i;
                    //         }
                    //     }
                    //     else if( dz < min_dz )
                    //     {
                    //             min_dz = dz;
                    //             vtx_index = i;
                    //     }
                    // }
                }
            }
                        
            //---use highest ranked if muon doesn't belong to any vertex
            const reco::Vertex& vtx = (vtx_index == -1 ? (*vtxHandle_)[0] : (*vtxHandle_)[vtx_index]);
            int n_close_vtxs=0;
            for( unsigned i = 0; i < vtxHandle_->size(); ++i )
            {
                if(i == (vtx_index == -1 ? 0 : std::abs(vtx_index)))
                    continue;
                const auto& ivtx = (*vtxHandle_)[i];
                const float dz = std::abs(vtx.z() - ivtx.z());
                if( dz < dzCut_ )
                    ++n_close_vtxs;
            }
            
            const auto tracks_removed = vertices_cleaned_tracks_zt3.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_z  = vertices_to_tracks_z.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_z1  = vertices_to_tracks_z1.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_z3  = vertices_to_tracks_z3.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_z5  = vertices_to_tracks_z5.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_z7  = vertices_to_tracks_z7.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_z10  = vertices_to_tracks_z10.equal_range(vtx_index == -1 ? 0 : vtx_index);                        
            const auto tracks_zt3 = vertices_to_tracks_zt3.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_zt4 = vertices_to_tracks_zt4.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_zt5 = vertices_to_tracks_zt5.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_zt7 = vertices_to_tracks_zt7.equal_range(vtx_index == -1 ? 0 : vtx_index);
            const auto tracks_zt10 = vertices_to_tracks_zt10.equal_range(vtx_index == -1 ? 0 : vtx_index);

            //---event counter
            outTrees_[iRes].iEvent = iEvent_;

            //---fill muon and vtx information
            outTrees_[iRes].iMuon = iMuon;
            outTrees_[iRes].pt = muon.pt();
            outTrees_[iRes].eta = muon.eta();
            outTrees_[iRes].phi = muon.phi();
            outTrees_[iRes].pz = muon.pz();
            outTrees_[iRes].dz = muon.track()->dz(vtx.position());
            outTrees_[iRes].dxy = muon.track()->dxy(vtx.position());
            outTrees_[iRes].isLooseMuon = muon::isLooseMuon(muon);
            outTrees_[iRes].isMediumMuon = muon::isMediumMuon(muon);
            outTrees_[iRes].isTightMuon = muon::isTightMuon(muon, (*vtxHandle_)[vtx_index==-1 ? 0 : vtx_index]);
            outTrees_[iRes].nVtxs = vtxHandle_->size();            
            outTrees_[iRes].vtxX = vtx.x();
            outTrees_[iRes].vtxY = vtx.y();
            outTrees_[iRes].vtxZ = vtx.z();
            outTrees_[iRes].vtxT = vtx.t();            
            outTrees_[iRes].genVtxZ = genPV.position().z();
            outTrees_[iRes].genVtxT = genPV.position().t();
            outTrees_[iRes].vtxIsFake = vtx.isFake();
            outTrees_[iRes].vtxIdx = vtx_index;            
            outTrees_[iRes].vtxNdof = vtx.ndof();
            outTrees_[iRes].vtxChi2 = vtx.chi2();
            outTrees_[iRes].nCloseVtxs = n_close_vtxs;

            //---per tracks info
            if(saveTracksInfo_)
            {
                std::swap(*outTrees_[iRes].tracksPt, tracks_pt);
                std::swap(*outTrees_[iRes].tracksEta, tracks_eta);
                std::swap(*outTrees_[iRes].tracksPhi, tracks_phi);                
                std::swap(*outTrees_[iRes].tracksDZ, tracks_dz);
                std::swap(*outTrees_[iRes].tracksDXY, tracks_dxy);
                std::swap(*outTrees_[iRes].tracksT, tracks_t);
                std::swap(*outTrees_[iRes].tracksKeepZ, tracks_keepz);
                std::swap(*outTrees_[iRes].tracksKeepT, tracks_keept);
            }
            
            //---compute the varius isolations for all requested cone sizes
            outTrees_[iRes].chIsoZCut->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZ1Cut->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZ3Cut->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZ5Cut->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZ7Cut->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZ10Cut->resize(isoConeSizes_.size(), 0);            
            outTrees_[iRes].chIsoZTCut_3sigma->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZTCut_4sigma->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZTCut_5sigma->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZTCut_7sigma->resize(isoConeSizes_.size(), 0);
            outTrees_[iRes].chIsoZTCut_10sigma->resize(isoConeSizes_.size(), 0);            
            // outTrees_[iRes].chIsoZCut->resize(isoConeSizes_.size(), isTimingSample_ ? -1 : 0);
            // outTrees_[iRes].chIsoZTCut_3sigma->resize(isoConeSizes_.size(), isTimingSample_ ? 0 : -1);
            // outTrees_[iRes].chIsoZTCut_4sigma->resize(isoConeSizes_.size(), isTimingSample_ ? 0 : -1);
            // outTrees_[iRes].chIsoZTCut_5sigma->resize(isoConeSizes_.size(), isTimingSample_ ? 0 : -1);
            // outTrees_[iRes].chIsoZTCut_7sigma->resize(isoConeSizes_.size(), isTimingSample_ ? 0 : -1);
            // outTrees_[iRes].chIsoZTCut_10sigma->resize(isoConeSizes_.size(), isTimingSample_ ? 0 : -1);            

            for(int iDR=0; iDR<(int)isoConeSizes_.size(); ++iDR)
            {
                float DR = isoConeSizes_[iDR];
                outTrees_[iRes].chIsoDR->push_back(DR);
                //---fixed dz only
                for( auto it = tracks_z.first; it != tracks_z.second; ++it )
                {
                    auto ref = it->second;
                    if( ref == muon.track() ) continue;
                    if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                    outTrees_[iRes].chIsoZCut->at(iDR) += ref->pt();
                    if(iDR==0)
                        outTrees_[iRes].tracksInConePt->push_back(ref->pt());
                }
                if(!isTimingSample_)
                {
                    //---dz 1 sigma
                    for( auto it = tracks_z1.first; it != tracks_z1.second; ++it )
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZ1Cut->at(iDR) += ref->pt();
                    }                    
                    //---dz 3 sigma
                    for( auto it = tracks_z3.first; it != tracks_z3.second; ++it )
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZ3Cut->at(iDR) += ref->pt();
                    }                                        
                    //---dz 5 sigma
                    for( auto it = tracks_z5.first; it != tracks_z5.second; ++it )
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZ5Cut->at(iDR) += ref->pt();
                    }                    
                    //---dz 7 sigma
                    for( auto it = tracks_z7.first; it != tracks_z7.second; ++it )
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZ7Cut->at(iDR) += ref->pt();
                    }                    
                    //---dz 10 sigma
                    for( auto it = tracks_z10.first; it != tracks_z10.second; ++it )
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZ10Cut->at(iDR) += ref->pt();
                    }                    
                }                
                else if(n_close_vtxs > 0)
                {
                    //--- dz + dt 3 sigma
                    for( auto it = tracks_zt3.first; it != tracks_zt3.second; ++it )
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZTCut_3sigma->at(iDR) += ref->pt();
                    }
                    //---rejected tracks
                    if(iDR==0)
                    {
                        for( auto it = tracks_removed.first; it != tracks_removed.second; ++it )
                        {
                            auto ref = it->second;
                            if( ref == muon.track() ) continue;
                            if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                            outTrees_[iRes].tracksRemovedPt->push_back(ref->pt());
                        }
                    }

                    //--- dz + dt 4 sigma
                    for( auto it = tracks_zt4.first; it != tracks_zt4.second; ++it ) 
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZTCut_4sigma->at(iDR) += ref->pt();
                    }
                    
                    //--- dz + dt 5 sigma
                    for( auto it = tracks_zt5.first; it != tracks_zt5.second; ++it ) 
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZTCut_5sigma->at(iDR) += ref->pt();
                    }

                    //--- dz + dt 7 sigma
                    for( auto it = tracks_zt7.first; it != tracks_zt7.second; ++it ) 
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZTCut_7sigma->at(iDR) += ref->pt();
                    }
                    
                    //--- dz + dt 10 sigma
                    for( auto it = tracks_zt10.first; it != tracks_zt10.second; ++it ) 
                    {
                        auto ref = it->second;
                        if( ref == muon.track() ) continue;
                        if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                        outTrees_[iRes].chIsoZTCut_10sigma->at(iDR) += ref->pt();
                    }
                }
                else
                {
                    outTrees_[iRes].chIsoZTCut_3sigma->at(iDR) = outTrees_[iRes].chIsoZCut->at(iDR);
                    outTrees_[iRes].chIsoZTCut_4sigma->at(iDR) = outTrees_[iRes].chIsoZCut->at(iDR);
                    outTrees_[iRes].chIsoZTCut_5sigma->at(iDR) = outTrees_[iRes].chIsoZCut->at(iDR);
                    outTrees_[iRes].chIsoZTCut_7sigma->at(iDR) = outTrees_[iRes].chIsoZCut->at(iDR);
                    outTrees_[iRes].chIsoZTCut_10sigma->at(iDR) = outTrees_[iRes].chIsoZCut->at(iDR);
                }
            }
            
            //---Gen matching            
            outTrees_[iRes].genMatched = false;
            outTrees_[iRes].genMatchedPrompt = false;
            outTrees_[iRes].genMatchedJet = false;
            outTrees_[iRes].genPt = -99.;
            outTrees_[iRes].genEta = -99.;
            outTrees_[iRes].genPhi = -99.;
            outTrees_[iRes].genIso = 0.;

            double mindr = std::numeric_limits<double>::max();
            double isocone=0.3*0.3;
            for (const auto& p : *genPartHandle_)
            {
                if (p.status() != 1 || !p.isLastCopy() ) continue;
                if (std::abs(p.pdgId()) == 13)
                {      
                    double dr = reco::deltaR2(muon,p);
                    if( dr<0.2*0.2 && dr<mindr )
                    {
                        mindr = dr;	
                        outTrees_[iRes].genMatched = true;
                        outTrees_[iRes].genMatchedPrompt = p.isPromptFinalState();
                        outTrees_[iRes].genPt = p.pt();
                        outTrees_[iRes].genEta = p.eta();
                        outTrees_[iRes].genPhi = p.phi();
                    }
                }
                else if(reco::deltaR2(muon,p)<isocone && p.charge()!=0 && p.pdgId()>40)
                    outTrees_[iRes].genIso += p.pt();
                    
            }

            mindr = std::numeric_limits<double>::max();
            for( const auto& jet : *genJetHandle_ ) {
                if( jet.pt() < 15.0 || jet.hadEnergy()/jet.energy() < 0.3)
                    continue;
                double dr = reco::deltaR2(muon,jet);
                if( dr < isocone && dr < mindr )
                {
                    outTrees_[iRes].genMatchedJet = true;
                    outTrees_[iRes].genPt = jet.pt();
                    outTrees_[iRes].genEta = jet.eta();
                    outTrees_[iRes].genPhi = jet.phi();                
                }
            }

            //---Fill tree
            outTrees_[iRes].GetTTreePtr()->Fill();
        }
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
