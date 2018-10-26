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

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

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
    edm::EDGetTokenT<genXYZ>                            genXYZToken_;
    edm::Handle<genXYZ>                                 genXYZHandle_;
    edm::EDGetTokenT<float>                             genT0Token_;    
    edm::Handle<float>                                  genT0Handle_;
    edm::EDGetTokenT<vector<SimVertex> >                simVtxToken_;
    edm::Handle<vector<SimVertex> >                     simVtxHandle_;    
    edm::EDGetTokenT<reco::VertexCollection>            vtx3DToken_;
    edm::Handle<reco::VertexCollection>                 vtx3DHandle_;    
    edm::EDGetTokenT<reco::VertexCollection>            vtx4DToken_;
    edm::Handle<reco::VertexCollection>                 vtx4DHandle_;    
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

    //---
    vector<double> targetResolutions_;
    double dzCut_;

    
    //---I/O
    int iEvent_;
    edm::Service<TFileService> fs;
    map<double, FTLMuonIsoTree> outTrees_;

    //---options
    bool           useMCTruthPV_;
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
    simVtxToken_(consumes<vector<SimVertex> >(pSet.getUntrackedParameter<edm::InputTag>("genVtxTag"))),    
    vtx3DToken_(consumes<std::vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtxTag3D"))),    
    vtx4DToken_(consumes<std::vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtxTag4D"))),    
    muonsToken_(consumes<reco::MuonCollection>(pSet.getUntrackedParameter<edm::InputTag>("muonsTag"))),
    tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("tracksTag"))),
    timeToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("timeTag"))),
    timeResToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("timeResTag"))),
    genPartToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genPartTag"))),
    genJetToken_(consumes<std::vector<reco::GenJet> >(pSet.getUntrackedParameter<edm::InputTag>("genJetsTag"))),
    targetResolutions_(pSet.getUntrackedParameter<vector<double> >("targetResolutions")),
    dzCut_(pSet.getUntrackedParameter<double>("dzCut"))
{
    iEvent_ = 0;
    for(auto& res : targetResolutions_)
        outTrees_[res] = FTLMuonIsoTree((pSet.getUntrackedParameter<string>("treeName")+"_"+to_string(int(res*1000))).c_str(),
                                         "Muon tree for FTL studies");
    useMCTruthPV_ = pSet.getUntrackedParameter<bool>("useMCTruthPV");
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
    iEvent.getByToken(genXYZToken_, genXYZHandle_);
    iEvent.getByToken(genT0Token_, genT0Handle_);
    iEvent.getByToken(muonsToken_, muonsHandle_);
    iEvent.getByToken(tracksToken_, tracksHandle_);
    iEvent.getByToken(timeToken_, timeHandle_);
    iEvent.getByToken(timeResToken_, timeResHandle_);
    iEvent.getByToken(simVtxToken_, simVtxHandle_);    
    iEvent.getByToken(vtx4DToken_, vtx4DHandle_);
    iEvent.getByToken(vtx3DToken_, vtx3DHandle_);
    iEvent.getByToken(genPartToken_, genPartHandle_);
    iEvent.getByToken(genJetToken_, genJetHandle_);

    edm::ESHandle<TransientTrackBuilder> theTTBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theTTBuilder);

    //---skip bad events (checks on muons and MC-truth vtx)
    unsigned int nmuons = 0;
    for (const auto &muon : *muonsHandle_)
        if (muon.pt() > 5.)
            ++nmuons;
  
    // if(!nmuons || vtx4DHandle_->size()<2 || vtx3DHandle_->size()<2 ||
    //    simVtxHandle_.product()->size() == 0 || simVtxHandle_.product()->at(0).vertexId() != 0
    //     )
    //     return;
    ++iEvent_;
    
    //---get truth PV
    SimVertex simPV;
    if(simVtxHandle_.isValid())
        simPV = simVtxHandle_.product()->at(0);
    else
    {
        auto xyz = genXYZHandle_.product();
        auto t = *genT0Handle_.product();
        auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
        simPV = SimVertex(v, t);
    }

    for(auto& iRes : targetResolutions_)
    {
        // cout << "RESOLUTION : " << iRes << endl;
        //---make a map of vertices to track refs within cuts
        std::unordered_multimap<unsigned,reco::TrackBaseRef>
            vertices_to_tracks_z,
            vertices_to_tracks_zt3,
            vertices_to_tracks_zt4,
            vertices_to_tracks_zt5,
            vertices_to_tracks_zt7,
            vertices_to_tracks_zt10;
        std::map<reco::TrackRef, float> track_times;
        std::map<reco::TrackRef, float> track_timeResos;

        // float extra_smearing = std::sqrt(iRes*iRes - 0.02*0.02);
        for(unsigned i = 0; i < tracksHandle_->size(); ++i)
        {
            auto ref = tracksHandle_->refAt(i);
            float time = (*timeHandle_)[ref];
            float timeReso = (*timeResHandle_)[ref] != 0.f ? (*timeResHandle_)[ref] : 0.170f;
            float extra_smearing = std::sqrt(iRes*iRes - timeReso*timeReso);
            // cout << "Track " << i << ":\n";
            // cout << "              iRes: " << iRes << endl;
            // cout << "    extra_smearing: " << extra_smearing << endl;
            // cout << "     orig timeReso: " << timeReso << endl;
            // cout << "         orig time: " << time << endl;
            time += gRandom->Gaus(0., extra_smearing);
            timeReso = std::sqrt(timeReso*timeReso + extra_smearing*extra_smearing);
            // cout << "          new time: " << time << endl;
            // cout << "      new timeReso: " << timeReso << endl;
            track_times[ref.castTo<reco::TrackRef>()] = time;
            track_timeResos[ref.castTo<reco::TrackRef>()] = timeReso;
            for(int ivtx = 0; ivtx < (int)vtx4DHandle_->size(); ++ivtx)
            {
                const auto& thevtx = (*vtx4DHandle_)[ivtx];
                const float dz = std::abs(ref->dz(thevtx.position()));
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
                const float time_cut10 = 10.f*base_cut;

                const bool keepz = ( dz < dzCut_ );
                const bool keept3 = (!useTime || std::isnan(dt) || dt < time_cut3);
                const bool keept4 = (!useTime || std::isnan(dt) || dt < time_cut4);
                const bool keept5 = (!useTime || std::isnan(dt) || dt < time_cut5);
                const bool keept7 = (!useTime || std::isnan(dt) || dt < time_cut7);
                const bool keept10 = (!useTime || std::isnan(dt) || dt < time_cut10);
                // const bool keept3 = (!useTime || std::isnan(dt) || dt2/time_cut3 + dz2/(0.2*0.2) < 1);
                // const bool keept4 = (!useTime || std::isnan(dt) || dt2/time_cut4 + dz2/(0.2*0.2) < 1);
                // const bool keept5 = (!useTime || std::isnan(dt) || dt2/time_cut5 + dz2/(0.2*0.2) < 1);
                // const bool keept7 = (!useTime || std::isnan(dt) || dt2/time_cut7 + dz2/(0.2*0.2) < 1);
                // const bool keept10 = (!useTime || std::isnan(dt) || dt2/time_cut10 + dz2/(0.2*0.2) < 1);            
            
                if( ref->quality(reco::TrackBase::highPurity) && keepz ) {
                    //vertices_to_tracks_z.emplace(ivtx, ref);        
                    if( keept10 ) {
                        vertices_to_tracks_zt10.emplace(ivtx, ref);
                    }
                    if( keept7 ) {
                        vertices_to_tracks_zt7.emplace(ivtx, ref);
                    }
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
            
            for(int ivtx = 0; ivtx < (int)vtx3DHandle_->size(); ++ivtx)
            {
                const auto& thevtx = (*vtx3DHandle_)[ivtx];
                const float dz = std::abs(ref->dz(thevtx.position()));
                const bool keepz = ( dz < dzCut_ );
                if( ref->quality(reco::TrackBase::highPurity) && keepz ) {
                    vertices_to_tracks_z.emplace(ivtx, ref);        
                }
            }
        }
  
        for(auto &muon : *muonsHandle_)
        {
            //---reset output
            outTrees_[iRes].Reset();

            //---fill global info            
            outTrees_[iRes].event = iEvent.id().event();
            outTrees_[iRes].lumi = iEvent.id().luminosityBlock();
            outTrees_[iRes].run = iEvent.id().run();    
        
            //---basic check
            if(muon.track().isNull() || muon.pt() < 10)
                continue;

            auto muontrackref = muon.get<reco::TrackRef>();
            reco::TransientTrack ttmuon = theTTBuilder->build(muontrackref);

            int vtx3D_index = -1;
            float sip3D=-1;

            //---find the 3D vertex this muon is best associated to..
            if(useMCTruthPV_)
            {
                float min_dz = 999.;
                for( unsigned i = 0; i < vtx3DHandle_->size(); ++i )
                {
                    const auto& vtx = (*vtx3DHandle_)[i];      
                    const float dz = fabs(vtx.z() - simPV.position().z());
                    if( dz < min_dz )
                    {
                        min_dz = dz;
                        vtx3D_index = i;
                    }
                }                
            }
            else
            {
                float max_weight3D = 0.f;
                for( unsigned i = 0; i < vtx3DHandle_->size(); ++i )
                {
                    const auto& vtx = (*vtx3DHandle_)[i];      
                    const float weight = vtx.trackWeight(muon.track());
                    if( weight > max_weight3D )
                    {
                        max_weight3D = weight;
                        vtx3D_index = i;
                    }
                }
                if (vtx3D_index==-1){
                  float minSIP=std::numeric_limits<float>::max();
                  for (unsigned i = 0; i < vtx3DHandle_->size(); ++i){
                    const auto& vtx = (*vtx3DHandle_)[i];
                    std::pair<bool, Measurement1D> IP_Measurement = IPTools::signedImpactParameter3D(
                      ttmuon,
                      GlobalVector(muon.px(), muon.py(), muon.pz()),
                      vtx
                    );
                    float IP_Vtx = IP_Measurement.second.value();
                    float dIP_Vtx = IP_Measurement.second.error();

                    const float valSIP = (dIP_Vtx!=0. ? fabs(IP_Vtx)/dIP_Vtx : 0);
                    if (minSIP>valSIP){
                      minSIP = valSIP;
                      vtx3D_index = i;
                    }
                  }
                }
            }            
                
            //---use highest ranked if muon doesn't belong to any vertex
            const reco::Vertex& vtx3D = (vtx3D_index == -1 ? (*vtx3DHandle_)[0] : (*vtx3DHandle_)[vtx3D_index]);
            const auto tracks_z  = vertices_to_tracks_z.equal_range(vtx3D_index == -1 ? 0 : vtx3D_index);
            {
              std::pair<bool, Measurement1D> IP_Measurement = IPTools::signedImpactParameter3D(
                ttmuon,
                GlobalVector(muon.px(), muon.py(), muon.pz()),
                vtx3D
              );
              float IP_Vtx = IP_Measurement.second.value();
              float dIP_Vtx = IP_Measurement.second.error();
              sip3D = (dIP_Vtx!=0. ? fabs(IP_Vtx)/dIP_Vtx : 0);
            }
            
            int vtx4D_index = -1;
            float sip4D=-1;

            //---find the 4D vertex this muon is best associated to..
            if(useMCTruthPV_)
            {
                double min_dzdt = std::numeric_limits<double>::max();
                for( unsigned i = 0; i < vtx4DHandle_->size(); ++i )
                {
                    const auto& vtx = (*vtx4DHandle_)[i];
                    const float dz = std::abs(vtx.z() - simPV.position().z());
                    const double dzdt = pow((vtx.z() - simPV.position().z())/vtx.zError(), 2) +
                        pow((vtx.t()-simPV.position().t())/vtx.tError(), 2);
                    if( dz < 0.1 && dzdt < min_dzdt )
                    {
                        min_dzdt = dzdt;
                        vtx4D_index = i;
                    }
                }
            }
            else
            {
                float max_weight4D = 0.f;
                for( unsigned i = 0; i < vtx4DHandle_->size(); ++i )
                {
                    const auto& vtx = (*vtx4DHandle_)[i];      
                    const float weight = vtx.trackWeight(muon.track());
                    if( weight > max_weight4D )
                    {
                        max_weight4D = weight;
                        vtx4D_index = i;
                    }
                }    
                if (vtx4D_index==-1){
                  float minSIP=std::numeric_limits<float>::max();
                  for (unsigned i = 0; i < vtx4DHandle_->size(); ++i){
                    const auto& vtx = (*vtx4DHandle_)[i];
                    std::pair<bool, Measurement1D> IP_Measurement = IPTools::signedImpactParameter3D(
                      ttmuon,
                      GlobalVector(muon.px(), muon.py(), muon.pz()),
                      vtx
                    );
                    float IP_Vtx = IP_Measurement.second.value();
                    float dIP_Vtx = IP_Measurement.second.error();
                    if (vtx.t()!=0.){
                      IP_Vtx = sqrt(IP_Vtx*IP_Vtx + pow(vtx.t()-track_times[muontrackref], 2));
                      dIP_Vtx = sqrt(dIP_Vtx*dIP_Vtx + pow(track_timeResos[muontrackref], 2) + pow(vtx.tError(), 2));
                    }

                    const float valSIP = (dIP_Vtx!=0. ? fabs(IP_Vtx)/dIP_Vtx : 0);
                    if (minSIP>valSIP){
                      minSIP = valSIP;
                      vtx4D_index = i;
                    }
                  }
                }
            }
            
            //---use highest ranked if muon doesn't belong to any vertex
            const reco::Vertex& vtx4D = (vtx4D_index == -1 ? (*vtx4DHandle_)[0] : (*vtx4DHandle_)[vtx4D_index]);
            const auto tracks_zt3 = vertices_to_tracks_zt3.equal_range(vtx4D_index == -1 ? 0 : vtx4D_index);
            const auto tracks_zt4 = vertices_to_tracks_zt4.equal_range(vtx4D_index == -1 ? 0 : vtx4D_index);
            const auto tracks_zt5 = vertices_to_tracks_zt5.equal_range(vtx4D_index == -1 ? 0 : vtx4D_index);
            const auto tracks_zt7 = vertices_to_tracks_zt7.equal_range(vtx4D_index == -1 ? 0 : vtx4D_index);
            const auto tracks_zt10 = vertices_to_tracks_zt10.equal_range(vtx4D_index == -1 ? 0 : vtx4D_index);
            {
              std::pair<bool, Measurement1D> IP_Measurement = IPTools::signedImpactParameter3D(
                ttmuon,
                GlobalVector(muon.px(), muon.py(), muon.pz()),
                vtx4D
              );
              float IP_Vtx = IP_Measurement.second.value();
              float dIP_Vtx = IP_Measurement.second.error();
              if (vtx4D.t()!=0.){
                IP_Vtx = sqrt(IP_Vtx*IP_Vtx + pow(vtx4D.t()-track_times[muontrackref], 2));
                dIP_Vtx = sqrt(dIP_Vtx*dIP_Vtx + pow(track_timeResos[muontrackref], 2) + pow(vtx4D.tError(), 2));
              }
              sip4D = (dIP_Vtx!=0. ? fabs(IP_Vtx)/dIP_Vtx : 0);
            }

            //---muon - vtx matching
            auto muonTkRef = muon.track();
            float muonTime = (*timeHandle_)[muonTkRef];
            float muonTimeReso = (*timeResHandle_)[muonTkRef];
            if( useMCTruthPV_ && ( muonTimeReso==0 ||
                                   ( muon.track()->dz(Point(simPV.position().x(),
                                                            simPV.position().y(),
                                                            simPV.position().z())) > 0.1 &&
                                     std::abs(simPV.position().t()-muonTime) > 3*muonTimeReso ) ))
                continue;

            //---event counter
            outTrees_[iRes].iEvent = iEvent_;
            
            //---fill gen vtx info
            outTrees_[iRes].simPVX = simPV.position().x();
            outTrees_[iRes].simPVY = simPV.position().y();
            outTrees_[iRes].simPVZ = simPV.position().z();
            outTrees_[iRes].simPVT = simPV.position().t();

            //---fill muon and vtx information            
            outTrees_[iRes].pt = muon.pt();
            outTrees_[iRes].eta = muon.eta();
            outTrees_[iRes].phi = muon.phi();
            outTrees_[iRes].px = muon.px();
            outTrees_[iRes].py = muon.py();
            outTrees_[iRes].pz = muon.pz();
            outTrees_[iRes].muonZ = muon.track()->dz(vtx4D.position()) + vtx4D.z();
            outTrees_[iRes].isLooseMuon = muon::isLooseMuon(muon);
            outTrees_[iRes].isMediumMuon = muon::isMediumMuon(muon);
            outTrees_[iRes].isTightMuon = muon::isTightMuon(muon, (*vtx4DHandle_)[vtx4D_index==-1 ? 0 : vtx4D_index]);
            outTrees_[iRes].nVtxs = vtx4DHandle_->size();        
            outTrees_[iRes].vtx3DIdx = vtx3D_index;
            outTrees_[iRes].vtx4DIdx = vtx4D_index;
            outTrees_[iRes].vtx3DX = vtx3D.x();
            outTrees_[iRes].vtx3DY = vtx3D.y();
            outTrees_[iRes].vtx3DZ = vtx3D.z();
            outTrees_[iRes].SIP3D = sip3D;
            outTrees_[iRes].vtx3DIsFake = vtx3D.isFake();
            outTrees_[iRes].vtx3DNdof = vtx3D.ndof();
            outTrees_[iRes].vtx3DChi2 = vtx3D.chi2();
            //outTrees_[iRes].vtx3DT = vtx3D.t();
            outTrees_[iRes].vtx4DX = vtx4D.x();
            outTrees_[iRes].vtx4DY = vtx4D.y();
            outTrees_[iRes].vtx4DZ = vtx4D.z();
            outTrees_[iRes].vtx4DT = vtx4D.t();
            outTrees_[iRes].vtx4DTerr = vtx4D.tError();
            outTrees_[iRes].SIP4D = sip4D;
            outTrees_[iRes].vtx4DIsFake = vtx4D.isFake();
            outTrees_[iRes].vtx4DNdof = vtx4D.ndof();
            outTrees_[iRes].vtx4DChi2 = vtx4D.chi2();    

            //---compute the varius isolations for all requested cone sizes
            outTrees_[iRes].nTracksVtx = std::distance(tracks_z.first, tracks_z.second);
            outTrees_[iRes].chIsoZCut->resize(isoConeSizes_.size());
            outTrees_[iRes].chIsoZTCut_3sigma->resize(isoConeSizes_.size());
            outTrees_[iRes].chIsoZTCut_4sigma->resize(isoConeSizes_.size());
            outTrees_[iRes].chIsoZTCut_5sigma->resize(isoConeSizes_.size());                
            outTrees_[iRes].chIsoZTCut_7sigma->resize(isoConeSizes_.size());                
            outTrees_[iRes].chIsoZTCut_10sigma->resize(isoConeSizes_.size());                
            outTrees_[iRes].nTracksZCut->resize(isoConeSizes_.size());
            outTrees_[iRes].nTracksZTCut_3sigma->resize(isoConeSizes_.size());
            outTrees_[iRes].nTracksZTCut_4sigma->resize(isoConeSizes_.size());
            outTrees_[iRes].nTracksZTCut_5sigma->resize(isoConeSizes_.size());                
            outTrees_[iRes].nTracksZTCut_7sigma->resize(isoConeSizes_.size());                
            outTrees_[iRes].nTracksZTCut_10sigma->resize(isoConeSizes_.size());                
            
            outTrees_[iRes].tracksZCut_iso03_pt->clear();
            outTrees_[iRes].tracksZCut_iso03_eta->clear();
            outTrees_[iRes].tracksZCut_iso03_phi->clear();
            outTrees_[iRes].tracksZCut_iso03_dR->clear();
            outTrees_[iRes].tracksZCut_iso03_t->clear();
            outTrees_[iRes].tracksZCut_iso03_dz->clear();
            outTrees_[iRes].tracksZTCut_3sigma_iso03_pt->clear();
            outTrees_[iRes].tracksZTCut_3sigma_iso03_eta->clear();
            outTrees_[iRes].tracksZTCut_3sigma_iso03_phi->clear();
            outTrees_[iRes].tracksZTCut_3sigma_iso03_dR->clear();
            outTrees_[iRes].tracksZTCut_3sigma_iso03_t->clear();
            outTrees_[iRes].tracksZTCut_3sigma_iso03_dz->clear();

            for(int iDR=0; iDR<(int)isoConeSizes_.size(); ++iDR)
            {
                float DR = isoConeSizes_[iDR];
                outTrees_[iRes].chIsoDR->push_back(DR);

                outTrees_[iRes].chIsoZCut->at(iDR) = 0.0;
                outTrees_[iRes].chIsoZTCut_3sigma->at(iDR) = 0.0;
                outTrees_[iRes].chIsoZTCut_4sigma->at(iDR) = 0.0;
                outTrees_[iRes].chIsoZTCut_5sigma->at(iDR) = 0.0;
                outTrees_[iRes].chIsoZTCut_7sigma->at(iDR) = 0.0;
                outTrees_[iRes].chIsoZTCut_10sigma->at(iDR) = 0.0;
                outTrees_[iRes].nTracksZCut->at(iDR) = 0;
                outTrees_[iRes].nTracksZTCut_3sigma->at(iDR) = 0;
                outTrees_[iRes].nTracksZTCut_4sigma->at(iDR) = 0;
                outTrees_[iRes].nTracksZTCut_5sigma->at(iDR) = 0;
                outTrees_[iRes].nTracksZTCut_7sigma->at(iDR) = 0;
                outTrees_[iRes].nTracksZTCut_10sigma->at(iDR) = 0;

                //--- dz only
                for( auto it = tracks_z.first; it != tracks_z.second; ++it ) {
                    auto ref = it->second.castTo<reco::TrackRef>();
                    float this_dr = reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi());
                    if(ref->pt() > 0.90){
                        outTrees_[iRes].nTracksZCut->at(iDR) += 1;
                        if(iDR==0){
                            outTrees_[iRes].tracksZCut_iso03_pt  -> push_back(ref->pt());
                            outTrees_[iRes].tracksZCut_iso03_eta -> push_back(ref->eta());
                            outTrees_[iRes].tracksZCut_iso03_phi -> push_back(ref->phi());
                            outTrees_[iRes].tracksZCut_iso03_dR  -> push_back(this_dr);
                            outTrees_[iRes].tracksZCut_iso03_t   -> push_back(track_times[ref]);
                            outTrees_[iRes].tracksZCut_iso03_dz  -> push_back(std::abs(ref->dz(vtx3D.position())));
                        }
                    }
                    if( ref == muon.track() ) continue;
                    if( this_dr >= DR*DR ) continue;
                    outTrees_[iRes].chIsoZCut->at(iDR) += ref->pt();
                }

                //--- dz + dt 3 sigma
                for( auto it = tracks_zt3.first; it != tracks_zt3.second; ++it ) {
                    auto ref = it->second.castTo<reco::TrackRef>();
                    float this_dr = reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi());
                    if(ref->pt() > 0.90) {
                        outTrees_[iRes].nTracksZTCut_3sigma->at(iDR) += 1;
                        if(iDR==0){
                            outTrees_[iRes].tracksZTCut_3sigma_iso03_pt  -> push_back(ref->pt());
                            outTrees_[iRes].tracksZTCut_3sigma_iso03_eta -> push_back(ref->eta());
                            outTrees_[iRes].tracksZTCut_3sigma_iso03_phi -> push_back(ref->phi());
                            outTrees_[iRes].tracksZTCut_3sigma_iso03_dR  -> push_back(this_dr);
                            outTrees_[iRes].tracksZTCut_3sigma_iso03_t   -> push_back(track_times[ref]);
                            outTrees_[iRes].tracksZTCut_3sigma_iso03_dz  -> push_back(std::abs(ref->dz(vtx4D.position())));
                        }
                    }
                    if( ref == muon.track() ) continue;
                    if( this_dr >= DR*DR ) continue;
                    outTrees_[iRes].chIsoZTCut_3sigma->at(iDR) += ref->pt();
                }

                //--- dz + dt 4 sigma
                for( auto it = tracks_zt4.first; it != tracks_zt4.second; ++it ) {
                    auto ref = it->second.castTo<reco::TrackRef>();
                    if(ref->pt() > 0.90) outTrees_[iRes].nTracksZTCut_4sigma->at(iDR) += 1;
                    if( ref == muon.track() ) continue;
                    if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                    outTrees_[iRes].chIsoZTCut_4sigma->at(iDR) += ref->pt();
                }

                //--- dz + dt 5 sigma
                for( auto it = tracks_zt5.first; it != tracks_zt5.second; ++it ) {
                    auto ref = it->second.castTo<reco::TrackRef>();
                    if(ref->pt() > 0.90) outTrees_[iRes].nTracksZTCut_5sigma->at(iDR) += 1;
                    if( ref == muon.track() ) continue;
                    if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                    outTrees_[iRes].chIsoZTCut_5sigma->at(iDR) += ref->pt();
                }

                //--- dz + dt 7 sigma
                for( auto it = tracks_zt7.first; it != tracks_zt7.second; ++it ) {
                    auto ref = it->second.castTo<reco::TrackRef>();
                    if(ref->pt() > 0.90) outTrees_[iRes].nTracksZTCut_7sigma->at(iDR) += 1;
                    if( ref == muon.track() ) continue;
                    if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                    outTrees_[iRes].chIsoZTCut_7sigma->at(iDR) += ref->pt();
                }

                //--- dz + dt 10 sigma
                for( auto it = tracks_zt10.first; it != tracks_zt10.second; ++it ) {
                    auto ref = it->second.castTo<reco::TrackRef>();
                    if(ref->pt() > 0.90) outTrees_[iRes].nTracksZTCut_10sigma->at(iDR) += 1;
                    if( ref == muon.track() ) continue;
                    if( reco::deltaR2(ref->eta(), ref->phi(), muon.eta(), muon.phi()) >= DR*DR ) continue;
                    outTrees_[iRes].chIsoZTCut_10sigma->at(iDR) += ref->pt();
                }
            }

            //---Gen matching
            outTrees_[iRes].genMatched = false;
            outTrees_[iRes].genMatchedPrompt = false;
            outTrees_[iRes].genMatchedJet = false;
            outTrees_[iRes].genPt = -99.;
            outTrees_[iRes].genEta = -99.;
            outTrees_[iRes].genPhi = -99.;
    
            double mindr = std::numeric_limits<double>::max();

            reco::Vertex::Point matchedPartVtx;
            int globalMatchPartId=-99;
            for (const reco::GenParticle &p : *genPartHandle_){
              if (p.status() != 1) continue;

              double dr = reco::deltaR(muon, p);
              if (dr<mindr){
                mindr = dr;
                matchedPartVtx = p.vertex();
                globalMatchPartId = p.pdgId();
              }
            }
            outTrees_[iRes].globalMatchPartId = globalMatchPartId;
            outTrees_[iRes].globalMatchPartVtxX = matchedPartVtx.x();
            outTrees_[iRes].globalMatchPartVtxY = matchedPartVtx.y();
            outTrees_[iRes].globalMatchPartVtxZ = matchedPartVtx.z();


            mindr = std::numeric_limits<double>::max();
            for (const reco::GenParticle &p : *genPartHandle_)
            {
                if (p.status() != 1) continue;
                if (std::abs(p.pdgId()) != 13) continue;
      
                double dr = reco::deltaR(muon,p);
                if( dr<0.2 && dr<mindr )
                {
                    mindr = dr;	
                    outTrees_[iRes].genMatched = true;
                    outTrees_[iRes].genMatchedPrompt = p.isPromptFinalState();
                    outTrees_[iRes].genPt = p.pt();
                    outTrees_[iRes].genEta = p.eta();
                    outTrees_[iRes].genPhi = p.phi();
                }      
            }

            outTrees_[iRes].genJetE = -99.;
            outTrees_[iRes].genJetPt = -99.;
            outTrees_[iRes].genJetEta = -99.;
            outTrees_[iRes].genJetPhi = -99.;
            mindr = std::numeric_limits<double>::max();
            for( const auto& jet : *genJetHandle_ ) {
                if( jet.pt() < 15.0 || jet.hadEnergy()/jet.energy() < 0.3)
                    continue;
                double dr = reco::deltaR(muon,jet);
                if( dr < 0.3 && dr < mindr )
                {
                    mindr = dr;
                    outTrees_[iRes].genMatchedJet = true;
                    outTrees_[iRes].genJetE = jet.energy();
                    outTrees_[iRes].genJetPt = jet.pt();
                    outTrees_[iRes].genJetEta = jet.eta();
                    outTrees_[iRes].genJetPhi = jet.phi();                
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
