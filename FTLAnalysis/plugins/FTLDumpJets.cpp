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
#include "DataFormats/JetReco/interface/GenJet.h"
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

#include "PrecisionTiming/FTLAnalysis/interface/FTLJetsTree.h"

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
    edm::Handle<vector<reco::GenJet> > genJetsHandle_;
    edm::EDGetTokenT<vector<reco::GenJet> > genJetsToken_;    
    edm::Handle<std::vector<reco::Photon> > photonsHandle_;
    edm::EDGetTokenT<std::vector<reco::Photon> > photonsToken_;    
    
    //---options
    bool readFTLrecHits_;
    
    //---outputs
    FTLJetsTree outTree_;
    vector<TH2F*> hitsMapHistos_={5, NULL};
    edm::Service<TFileService> fs_;    
};

FTLDumpJets::FTLDumpJets(const edm::ParameterSet& pSet):
    simTkToken_(consumes<edm::SimTrackContainer>(pSet.getUntrackedParameter<edm::InputTag>("simTkTag"))),
    simVtxToken_(consumes<edm::SimVertexContainer>(pSet.getUntrackedParameter<edm::InputTag>("simVtxTag"))),
    ftlRecHitsToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("ftlRecHitsTag"))),
    jetsToken_(consumes<vector<reco::PFJet> >(pSet.getUntrackedParameter<edm::InputTag>("jetsTag"))),
    genJetsToken_(consumes<vector<reco::GenJet> >(pSet.getUntrackedParameter<edm::InputTag>("genJetsTag"))),    
    photonsToken_(consumes<std::vector<reco::Photon> >(pSet.getUntrackedParameter<edm::InputTag>("photonsTag"))),    
    readFTLrecHits_(pSet.getUntrackedParameter<bool>("readFTLRecHits"))
{
    outTree_ = FTLJetsTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "Jets tree for FTL studies");
    //hitsMapHistos_.resize(5, new TH2F());
    for(unsigned int i=0; i<hitsMapHistos_.size(); ++i)        
        hitsMapHistos_[i] = fs_->make<TH2F>((string("hitsMapHisto_")+to_string(i)).c_str(), "",
            720, -360.25, 360.25, 101, -50.5, 50.5);
}

void FTLDumpJets::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    outTree_.Reset();
    
    //---load the jets
    event.getByToken(jetsToken_, jetsHandle_);
    auto jets = *jetsHandle_.product();

    //---load the gen jets
    // event.getByToken(genJetsToken_, genJetsHandle_);
    // auto genJets = *genJetsHandle_.product();
    
    //---load the photons
    // event.getByToken(photonsToken_, photonsHandle_);
    // auto photons = *photonsHandle_.product();    

    //---load the FTL collection if present in the EventContent (avoid crash with standard geometry)
    // auto ftlRecHits = FTLRecHitCollection();
    // if(readFTLrecHits_)
    //     event.getByToken(ftlRecHitsToken_, ftlRecHitsHandle_);
    // if(ftlRecHitsHandle_.isValid())
    //     ftlRecHits = *ftlRecHitsHandle_.product();

    int idx=0;

    // auto genJet1 = genJets[0];
    // auto genJet2 = genJets[1];    
    for(auto& jet : jets)
    {        
        //---skim
        if(fabs(jet.eta())>1.5 || jet.pt()<30 || jet.pt()>=1000)
        {
            ++idx;
            continue;
        }
        if(idx>1)
            break;

        auto tmpHitsHisto = (TH2F*)hitsMapHistos_[0]->Clone("tmp");
        tmpHitsHisto->Reset();
        
        //---jet standard info
        outTree_.idx->push_back(idx);
        outTree_.pt->push_back(jet.pt());
        outTree_.eta->push_back(jet.eta());
        outTree_.phi->push_back(jet.phi());
        outTree_.energy->push_back(jet.energy());
        // float minDR=1000;
        // for(auto& pho : photons)
        // {
        //     if(pho.pt()<10)
        //         break;
        //     if(deltaR(pho.eta(), pho.phi(), jet.eta(), jet.phi()) < minDR)
        //         minDR = deltaR(pho.eta(), pho.phi(), jet.eta(), jet.phi());
        // }
        outTree_.phoEfrac->push_back(jet.getTrackRefs().size());

        // if(minDR < 0.4)
        //     continue;
        
        //---tracks
        map<int, int> filledBins;
        map<int, double> timeDelta;
        for(auto& track : jet.getTrackRefs())
        {
            float delta_phi = deltaPhi(track->outerPhi(), jet.phi())/TMath::Pi()*360;
            float iz_jet = (jet.eta()/fabs(jet.eta())) * (1189/tan(2*atan(exp(fabs(jet.eta()))))-0.5)/10;
            float iz_track = (track->outerEta()/fabs(track->outerEta())) * (1189/tan(2*atan(exp(fabs(track->outerEta()))))-0.5)/10;
            float delta_iz = signbit(jet.eta()) == signbit(track->outerEta()) ?
                iz_track - iz_jet : iz_track - iz_jet -1;
            int cell = tmpHitsHisto->Fill(delta_phi, delta_iz, 1);
            filledBins[cell]++;
            // if(timeDelta.find(cell) != timeDelta.end())
            //     timeDelta[cell] -= 
        }

        ++idx;
        for(auto& ibin : filledBins)
        {
            if(ibin.second>1)
                hitsMapHistos_[int(jet.pt()/200)]->SetBinContent(ibin.first,
                                                                 hitsMapHistos_[int(jet.pt()/200)]->GetBinContent(ibin.first)+1);
        }
        tmpHitsHisto->Delete();
    }

    outTree_.GetTTreePtr()->Fill();
}

DEFINE_FWK_MODULE(FTLDumpJets);

#endif
