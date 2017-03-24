#ifndef _FTL_DUMP_TRACKS_
#define _FTL_DUMP_TRACKS_

#include "TMath.h"

#include "FWCore/Utilities/interface/BranchType.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/ForwardDetId/interface/FastTimeDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/FastTimeGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "PrecisionTiming/FTLAnalysis/interface/FTLTracksTree.h"

using namespace std;

class FTLDumpTracks : public edm::EDAnalyzer
{
public:
    explicit FTLDumpTracks(const edm::ParameterSet& pSet);
    ~FTLDumpTracks() {};

    //---utils

    //---methods
    virtual void beginJob() override {};
    virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
    virtual void endJob() override {};


private:
    //---inputs
    edm::Handle<reco::GenParticleCollection>          genParticlesHandle_;
    edm::EDGetTokenT<reco::GenParticleCollection>     genParticlesToken_;
    // edm::Handle<edm::SimTrackContainer>            simTkHandle_;
    // edm::EDGetTokenT<edm::SimTrackContainer>       simTkToken_;
    // edm::Handle<edm::SimVertexContainer>           simVtxHandle_;
    // edm::EDGetTokenT<edm::SimVertexContainer>      simVtxToken_;
    edm::Handle<FTLUncalibratedRecHitCollection>      ftlRecHitsHandle_;
    edm::EDGetTokenT<FTLUncalibratedRecHitCollection> ftlRecHitsToken_;    
    edm::EDGetTokenT<edm::View<reco::Track> >         tracksToken_;
    edm::Handle<edm::View<reco::Track> >              tracksHandle_;    
            
    //---outputs
    FTLTracksTree outTree_;
    edm::Service<TFileService> fs_;
};

FTLDumpTracks::FTLDumpTracks(const edm::ParameterSet& pSet):
    genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
    // simTkToken_(consumes<edm::SimTrackContainer>(pSet.getUntrackedParameter<edm::InputTag>("simTkTag"))),
    // simVtxToken_(consumes<edm::SimVertexContainer>(pSet.getUntrackedParameter<edm::InputTag>("simVtxTag"))),
    ftlRecHitsToken_(consumes<FTLUncalibratedRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("ftlRecHitsTag"))),
    tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("tracksTag")))
{
    outTree_ = FTLTracksTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "Tracks tree for FTL studies");
}

void FTLDumpTracks::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    outTree_.Reset();

    //---load gen particles
    event.getByToken(genParticlesToken_, genParticlesHandle_);
    auto genParticles = *genParticlesHandle_.product();

    //---load the tracks
    event.getByToken(tracksToken_, tracksHandle_);
    auto tracks = *tracksHandle_.product();

    //---load the magnetic field, build the track propagator and the FTL cylinder
    edm::ESHandle<MagneticField> magFieldHandle; 
    setup.get<IdealMagneticFieldRecord>().get(magFieldHandle);
    auto magField = magFieldHandle.product();    
    //---FTL cylinder definition
    Surface::PositionType pos(0, 0, 0);
    Surface::RotationType rot;
    auto ftlCylinderEntry = Cylinder::build(118.6, pos, rot);
    auto ftlCylinderExit = Cylinder::build(119.0, pos, rot);
    
    //---load the mc-truth tracker collections
    // event.getByToken(simTkToken_, simTkHandle_);
    // event.getByToken(simVtxToken_, simVtxHandle_);

    //---load the FTL collection if present in the EventContent (avoid crash with standard geometry)
    event.getByToken(ftlRecHitsToken_, ftlRecHitsHandle_);    
    auto ftlRecHits = *ftlRecHitsHandle_.product();

    // //---get the ecal geometry
    // const CaloSubdetectorGeometry* ecalBarrelGeometry = NULL;
    // if(readFTLrecHits_)
    // {
    //     edm::ESHandle<CaloGeometry> caloGeoHandle;
    //     setup.get<CaloGeometryRecord>().get(caloGeoHandle);
    //     ecalBarrelGeometry = caloGeoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    // }

    //---get the FTL geometry
    const FastTimeGeometry* ftlGeometry = NULL;
    edm::ESHandle<FastTimeGeometry> ftlGeoHandle;
    setup.get<IdealGeometryRecord>().get("FastTimeBarrel", ftlGeoHandle);
    ftlGeometry = &(*ftlGeoHandle);
    
    int idx=0;
    outTree_.nTracks = tracks.size();
    for(auto& track : tracks)
    {
        //---track standard info
        outTree_.idx->push_back(idx);
        outTree_.charge->push_back(track.charge());
        outTree_.p->push_back(track.p());        
        outTree_.pt->push_back(track.pt());
        outTree_.eta->push_back(track.eta());
        outTree_.phi->push_back(track.phi());
        outTree_.out_pt->push_back(track.outerPt());
        outTree_.tk_chi2->push_back(track.chi2());

        //---find ftl associated hits
        outTree_.ftlHitsTrkIdx->resize(idx+1);
        outTree_.ftlHitsEnergy->resize(idx+1);
        outTree_.ftlHitsTime->resize(idx+1);
        outTree_.ftlHitsEta->resize(idx+1);
        outTree_.ftlHitsPhi->resize(idx+1);
        outTree_.ftlHitsZ->resize(idx+1);        
        outTree_.ftlHitsEnergySum->resize(idx+1);
        outTree_.ftlNHits->resize(idx+1);
        outTree_.ftlHits3x3Sum->resize(idx+1);
        outTree_.ftlNHits3x3->resize(idx+1);
        outTree_.ftlSeedIdx->resize(idx+1, -1);
        outTree_.ftlSeedEnergy->resize(idx+1, -1);
        outTree_.ftlSeedTime->resize(idx+1, -1);        
        outTree_.ftlClusNHits->resize(idx+1, 0);
        outTree_.ftlClusEnergy->resize(idx+1, 0);
        outTree_.ftlClusTime->resize(idx+1, 0);        
        
        //---propagate track to timing layer cylinder
        //---trajectory definition
        GlobalPoint startingPoint(track.vx(), track.vy(), track.vz());
        GlobalVector startingMomentum(track.px(), track.py(), track.pz());
        FreeTrajectoryState trajectory(startingPoint, startingMomentum, track.charge(), magField);
        //---propagation
        SteppingHelixPropagator tkPropagator(magField);
        auto propagatedTrackFTLEntry = tkPropagator.propagateWithPath(trajectory, *ftlCylinderEntry);
        auto propagatedTrackFTLExit = tkPropagator.propagateWithPath(trajectory, *ftlCylinderExit);        

        float trkZOutFTL=0, trkPhiOutFTL=0;
        if(propagatedTrackFTLEntry.first.isValid())
        {
            outTree_.trkZAtFTL->push_back(propagatedTrackFTLEntry.first.globalPosition().z());
            outTree_.trkEtaAtFTL->push_back(propagatedTrackFTLEntry.first.globalMomentum().eta());
            outTree_.trkPhiAtFTL->push_back(propagatedTrackFTLEntry.first.globalMomentum().phi());
        }
        if(propagatedTrackFTLExit.first.isValid())
        {
            trkZOutFTL = propagatedTrackFTLExit.first.globalPosition().z();
            //trkEtaOutFTL = propagatedTrackFTLExit.first.globalMomentum().eta();
            trkPhiOutFTL = propagatedTrackFTLExit.first.globalMomentum().phi();
        }
        
        outTree_.ftlTotHits = ftlRecHits.size();

        for(auto ftl_hit : ftlRecHits)
        {
            // //---very loose energy cut
            // if(ftl_hit.amplitude() < 0.1)
            //     continue;
                
            FastTimeDetId id = ftl_hit.id();
            // float ftl_hit_eta = -log(tan(atan(1189/(id.iz()*10-0.5))/2))*id.zside();
            // float ftl_hit_phi = id.iphi()<=360 ? id.iphi()/360.*TMath::Pi() : (id.iphi()-720)/360.*TMath::Pi();
            float ftl_hit_eta = ftlGeometry->getPosition(id).eta();
            float ftl_hit_phi = ftlGeometry->getPosition(id).phi();
            float ftl_hit_z = ftlGeometry->getPosition(id).z();

            if(deltaR(ftl_hit_eta, ftl_hit_phi, track.eta(), track.phi()) < 0.1)
            {
                outTree_.ftlHitsTrkIdx->at(idx).push_back(idx);                
                outTree_.ftlHitsEnergy->at(idx).push_back(ftl_hit.amplitude());
                outTree_.ftlHitsTime->at(idx).push_back(ftl_hit.time());
                outTree_.ftlHitsEta->at(idx).push_back(ftl_hit_eta);
                outTree_.ftlHitsPhi->at(idx).push_back(ftl_hit_phi);
                outTree_.ftlHitsZ->at(idx).push_back(ftl_hit_z);
                outTree_.ftlNHits->at(idx)++;
            }
            if(propagatedTrackFTLEntry.first.isValid() &&
               fabs(deltaPhi(ftl_hit_phi, outTree_.trkPhiAtFTL->back())) < M_PI/720.*100. &&
               fabs(ftl_hit_z-outTree_.trkZAtFTL->back()) < 0.5)
            {
                outTree_.ftlSeedIdx->at(idx) = outTree_.ftlNHits->at(idx)-1;                
                outTree_.ftlSeedEnergy->at(idx) = ftl_hit.amplitude();
                outTree_.ftlSeedTime->at(idx) = ftl_hit.time();
            }
            if((propagatedTrackFTLEntry.first.isValid() &&
                fabs(deltaPhi(ftl_hit_phi, outTree_.trkPhiAtFTL->back())) < M_PI/720.*100. &&
                fabs(ftl_hit_z-outTree_.trkZAtFTL->back()) <= 0.5) ||
               (propagatedTrackFTLExit.first.isValid() &&
                fabs(deltaPhi(ftl_hit_phi, trkPhiOutFTL)) <= M_PI/720.*100. &&
                fabs(ftl_hit_z-trkZOutFTL) <= 0.5))
            {
                outTree_.ftlClusNHits->at(idx)++;
                outTree_.ftlClusEnergy->at(idx) += ftl_hit.amplitude();
                outTree_.ftlClusTime->at(idx) += ftl_hit.time()*ftl_hit.amplitude();
            }

        }
        outTree_.ftlClusTime->at(idx) /= outTree_.ftlClusEnergy->at(idx);
        outTree_.ftlClusEnergy->at(idx) /= outTree_.ftlClusNHits->at(idx);
        
        ++idx;
    }

    outTree_.GetTTreePtr()->Fill();
}

DEFINE_FWK_MODULE(FTLDumpTracks);

#endif
