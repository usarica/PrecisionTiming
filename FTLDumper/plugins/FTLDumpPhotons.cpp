#ifndef _FTL_DUMP_PHOTONS_
#define _FTL_DUMP_PHOTONS_

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
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"

#include "PrecisionTiming/FTLDumper/interface/FTLPhotonsTree.h"

using namespace std;

class FTLDumpPhotons : public edm::EDAnalyzer
{
public:
    explicit FTLDumpPhotons(const edm::ParameterSet& pSet);
    ~FTLDumpPhotons() {};

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
    edm::Handle<pat::PhotonCollection> photonsHandle_;
    edm::EDGetTokenT<pat::PhotonCollection> photonsToken_;

    //---options
    float mcTruthPhoEtThr_;
    
    //---workers
    PhotonMCTruthFinder photonMCTruthFinder_;
    
    //---outputs
    FTLPhotonsTree outTree_;
    edm::Service<TFileService> fs_;
};

FTLDumpPhotons::FTLDumpPhotons(const edm::ParameterSet& pSet):
    simTkToken_(consumes<edm::SimTrackContainer>(pSet.getUntrackedParameter<edm::InputTag>("simTkTag"))),
    simVtxToken_(consumes<edm::SimVertexContainer>(pSet.getUntrackedParameter<edm::InputTag>("simVtxTag"))),
    photonsToken_(consumes<pat::PhotonCollection>(pSet.getUntrackedParameter<edm::InputTag>("photonsTag"))),
    mcTruthPhoEtThr_(pSet.getUntrackedParameter<double>("mcTruthPhoEtThr")),
    photonMCTruthFinder_()    
{
    outTree_ = FTLPhotonsTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "Photons tree for FTL studies");
}


void FTLDumpPhotons::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    outTree_.Reset();
    
    event.getByToken(photonsToken_, photonsHandle_);
    event.getByToken(simTkToken_, simTkHandle_);
    event.getByToken(simVtxToken_, simVtxHandle_);
    auto photons = *photonsHandle_.product();

    int idx=0;
    for(auto& pho : photons)
    {
        //---photon standard info
        outTree_.idx->push_back(idx);
        outTree_.pt->push_back(pho.pt());
        outTree_.eta->push_back(pho.eta());
        outTree_.phi->push_back(pho.phi());
        outTree_.energy->push_back(pho.energy());
        outTree_.sc_energy->push_back(pho.superCluster()->energy());
        outTree_.r9->push_back(pho.r9());
        outTree_.sIeIe->push_back(pho.sigmaIetaIeta());
        
        //---MC truth conversion info
        auto mcTruthPhotons = photonMCTruthFinder_.find(*simTkHandle_.product(), *simVtxHandle_.product());
        PhotonMCTruth mcTruthPho;
        bool mcTruthExist=false;
        float minDR=10;
        for(auto& mcpho : mcTruthPhotons)
        {
            if(mcpho.fourMomentum().et() > mcTruthPhoEtThr_  &&
               deltaR(mcpho.fourMomentum().eta(), mcpho.fourMomentum().phi(), pho.eta(), pho.phi()) < minDR)
            {
                minDR = deltaR(mcpho.fourMomentum().eta(), mcpho.fourMomentum().phi(), pho.eta(), pho.phi());
                mcTruthPho = mcpho;
                mcTruthExist = true;
            }
        }        
        if(mcTruthExist)
            outTree_.convRadius->push_back(mcTruthPho.vertex().vect().rho());
        else
            outTree_.convRadius->push_back(-1);;
            
        ++idx;
    }
    if(photons.size() >= 2)
        outTree_.mass = (photons[0].p4() + photons[1].p4()).mass();
    else
        outTree_.mass = -1;

    outTree_.GetTTreePtr()->Fill();
}
    
DEFINE_FWK_MODULE(FTLDumpPhotons);

#endif
