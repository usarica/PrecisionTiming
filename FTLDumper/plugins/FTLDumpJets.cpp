#ifndef _FTL_DUMP_JETS_
#define _FTL_DUMP_JETS_

#include "TMath.h"
#include "TH2.h"

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

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ForwardDetId/interface/FastTimeDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "Geometry/Records/interface/FastTimeGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/FastTimeGeometry.h"

#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"

#include "PrecisionTiming/FTLDumper/interface/FTLJetsTree.h"

using namespace std;

class FTLDumpJets : public edm::EDAnalyzer
{
public:
    explicit FTLDumpJets(const edm::ParameterSet& pSet);
    ~FTLDumpJets() {};

    //---utils

    //---methods
    virtual void beginJob() override {};
    virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
    virtual void endJob() override {};


private:
    //---inputs
    edm::Handle<edm::SimTrackContainer> simTkHandle_;
    edm::EDGetTokenT<edm::SimTrackContainer> simTkToken_;
    edm::Handle<edm::SimVertexContainer> simVtxHandle_;
    edm::EDGetTokenT<edm::SimVertexContainer> simVtxToken_;
    edm::Handle<FTLRecHitCollection> ftlRecHitsHandle_;
    edm::EDGetTokenT<FTLRecHitCollection> ftlRecHitsToken_;    
    edm::Handle<vector<reco::PFJet> > jetsHandle_;
    edm::EDGetTokenT<vector<reco::PFJet> > jetsToken_;    
    edm::Handle<std::vector<reco::Photon> > photonsHandle_;
    edm::EDGetTokenT<std::vector<reco::Photon> > photonsToken_;    
    
    //---options
    float mcTruthPhoEtThr_;
    bool readFTLrecHits_;
    
    //---workers
    PhotonMCTruthFinder photonMCTruthFinder_;
    
    //---outputs
    FTLJetsTree outTree_;
    TH2F* hitsMapHisto_;
    edm::Service<TFileService> fs_;
};

FTLDumpJets::FTLDumpJets(const edm::ParameterSet& pSet):
    simTkToken_(consumes<edm::SimTrackContainer>(pSet.getUntrackedParameter<edm::InputTag>("simTkTag"))),
    simVtxToken_(consumes<edm::SimVertexContainer>(pSet.getUntrackedParameter<edm::InputTag>("simVtxTag"))),
    ftlRecHitsToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("ftlRecHitsTag"))),
    jetsToken_(consumes<vector<reco::PFJet> >(pSet.getUntrackedParameter<edm::InputTag>("jetsTag"))),
    photonsToken_(consumes<std::vector<reco::Photon> >(pSet.getUntrackedParameter<edm::InputTag>("photonsTag"))),    
    mcTruthPhoEtThr_(pSet.getUntrackedParameter<double>("mcTruthPhoEtThr")),
    readFTLrecHits_(pSet.getUntrackedParameter<bool>("readFTLRecHits")),
    photonMCTruthFinder_()    
{
    outTree_ = FTLJetsTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "Jets tree for FTL studies");
    hitsMapHisto_ = fs_->make<TH2F>("hitsMapHisto", "", 720, -360.25, 360.25, 101, -50.5, 50.5);
}

void FTLDumpJets::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    outTree_.Reset();

    //---load the jets
    event.getByToken(jetsToken_, jetsHandle_);
    auto jets = *jetsHandle_.product();

    //---load the photons
    event.getByToken(photonsToken_, photonsHandle_);
    auto photons = *photonsHandle_.product();    

    //---load the mc-truth tracker collections
    event.getByToken(simTkToken_, simTkHandle_);
    event.getByToken(simVtxToken_, simVtxHandle_);

    //---load the FTL collection if present in the EventContent (avoid crash with standard geometry)
    auto ftlRecHits = FTLRecHitCollection();
    if(readFTLrecHits_)
        event.getByToken(ftlRecHitsToken_, ftlRecHitsHandle_);
    if(ftlRecHitsHandle_.isValid())
        ftlRecHits = *ftlRecHitsHandle_.product();

    int idx=0;
    for(auto& jet : jets)
    {
        //---skim
        if(fabs(jet.eta())>1.5 || jet.pt()<30 || jet.getTrackRefs().size()<5)
            continue;
        
        //---jet standard info
        outTree_.idx->push_back(idx);
        outTree_.pt->push_back(jet.pt());
        outTree_.eta->push_back(jet.eta());
        outTree_.phi->push_back(jet.phi());
        outTree_.energy->push_back(jet.energy());
        float minDR=1000;
        for(auto& pho : photons)
        {
            if(pho.pt()<10)
                break;
            if(deltaR(pho.eta(), pho.phi(), jet.eta(), jet.phi()) < minDR)
                minDR = deltaR(pho.eta(), pho.phi(), jet.eta(), jet.phi());
        }
        outTree_.phoEfrac->push_back(minDR);

        if(minDR < 0.4)
            continue;
        
        //---tracks
        for(auto& track : jet.getTrackRefs())
        {
            if(deltaR(jet.eta(), track->eta(), jet.phi(), track->phi())<0.4)
            {
                float delta_phi = deltaPhi(track->outerPhi(), jet.phi())/TMath::Pi()*360;
                float iz_jet = (jet.eta()/fabs(jet.eta())) * (1189/tan(2*atan(exp(fabs(jet.eta()))))-0.5)/10;
                float iz_track = (track->outerEta()/fabs(track->outerEta())) * (1189/tan(2*atan(exp(fabs(track->outerEta()))))-0.5)/10;
                float delta_iz = signbit(jet.eta()) == signbit(track->outerEta()) ?
                    iz_track - iz_jet : iz_track - iz_jet -1;
                hitsMapHisto_->Fill(delta_phi, delta_iz, 1);
            }
        }
        ++idx;
    }

    outTree_.GetTTreePtr()->Fill();
}

DEFINE_FWK_MODULE(FTLDumpJets);

#endif
