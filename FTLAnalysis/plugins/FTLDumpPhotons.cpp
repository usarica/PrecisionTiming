#ifndef _FTL_DUMP_PHOTONS_
#define _FTL_DUMP_PHOTONS_

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

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/ForwardDetId/interface/FastTimeDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/FastTimeGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"

#include "PrecisionTiming/FTLAnalysis/interface/FTLPhotonsTree.h"

using namespace std;

template<class PhotonCollectionT>
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
    edm::Handle<FTLRecHitCollection> ftlRecHitsHandle_;
    edm::EDGetTokenT<FTLRecHitCollection> ftlRecHitsToken_;    
    edm::Handle<PhotonCollectionT> photonsHandle_;
    edm::EDGetTokenT<PhotonCollectionT> photonsToken_;    
    
    //---options
    float mcTruthPhoEtThr_;
    bool readFTLrecHits_;
    
    //---workers
    PhotonMCTruthFinder photonMCTruthFinder_;
    
    //---outputs
    FTLPhotonsTree outTree_;
    edm::Service<TFileService> fs_;
};

template<class PhotonCollectionT>
FTLDumpPhotons<PhotonCollectionT>::FTLDumpPhotons(const edm::ParameterSet& pSet):
    simTkToken_(consumes<edm::SimTrackContainer>(pSet.getUntrackedParameter<edm::InputTag>("simTkTag"))),
    simVtxToken_(consumes<edm::SimVertexContainer>(pSet.getUntrackedParameter<edm::InputTag>("simVtxTag"))),
    ftlRecHitsToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("ftlRecHitsTag"))),
    photonsToken_(consumes<PhotonCollectionT>(pSet.getUntrackedParameter<edm::InputTag>("photonsTag"))),
    mcTruthPhoEtThr_(pSet.getUntrackedParameter<double>("mcTruthPhoEtThr")),
    readFTLrecHits_(pSet.getUntrackedParameter<bool>("readFTLRecHits")),
    photonMCTruthFinder_()    
{
    outTree_ = FTLPhotonsTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "Photons tree for FTL studies");
}

template<class PhotonCollectionT>
void FTLDumpPhotons<PhotonCollectionT>::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    outTree_.Reset();

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

    //---get the ecal geometry
    edm::ESHandle<CaloGeometry> caloGeoHandle;
    setup.get<CaloGeometryRecord>().get(caloGeoHandle);
    const CaloSubdetectorGeometry* ecalBarrelGeometry =
        caloGeoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

    //---get the FTL geometry
    edm::ESHandle<FastTimeGeometry> ftlGeoHandle;
    setup.get<IdealGeometryRecord>().get("FastTimeBarrel", ftlGeoHandle);
    const FastTimeGeometry& ftlGeometry = *ftlGeoHandle;
    
    //---fill global info
    outTree_.event = event.id().event();
    outTree_.lumi = event.id().luminosityBlock();
    outTree_.run = event.id().run();    
    
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
        {
            outTree_.convRadius->push_back(mcTruthPho.vertex().vect().rho());
            outTree_.convPhi->push_back(mcTruthPho.vertex().vect().phi());
            outTree_.convZ->push_back(mcTruthPho.vertex().vect().z());            
            outTree_.genEnergy->push_back(mcTruthPho.fourMomentum().e());
            outTree_.genEta->push_back(mcTruthPho.fourMomentum().eta());
            outTree_.genPhi->push_back(mcTruthPho.fourMomentum().phi());
            outTree_.genPt->push_back(mcTruthPho.fourMomentum().et());
        }
        else
        {
            outTree_.convRadius->push_back(-1);;
            outTree_.genEnergy->push_back(-1);
            outTree_.genEta->push_back(-1);
            outTree_.genPhi->push_back(-1);
            outTree_.genPt->push_back(-1);
        }
        
        //---find ftl associated hits
        outTree_.ftlHitsPhoIdx->resize(idx+1);
        outTree_.ftlHitsEnergy->resize(idx+1);
        outTree_.ftlHitsTime->resize(idx+1);
        outTree_.ftlHitsEta->resize(idx+1);
        outTree_.ftlHitsPhi->resize(idx+1);
        outTree_.ftlHitsZ->resize(idx+1);        
        outTree_.ftlHitsEnergySum->resize(idx+1);
        outTree_.ftlNHits->resize(idx+1);
        outTree_.ftlHits3x3Sum->resize(idx+1);
        outTree_.ftlNHits3x3->resize(idx+1);

        outTree_.ftlTotHits = ftlRecHits.size();
        for(auto ftl_hit : ftlRecHits)
        {
            FastTimeDetId id = ftl_hit.id();
            // float ftl_hit_eta = -log(tan(atan(1189/(id.iz()*10-0.5))/2))*id.zside();
            // float ftl_hit_phi = id.iphi()<=360 ? id.iphi()/360.*TMath::Pi() : (id.iphi()-720)/360.*TMath::Pi();
            float ftl_hit_eta = ftlGeometry.getPosition(id).eta();
            float ftl_hit_phi = ftlGeometry.getPosition(id).phi();
            float ftl_hit_z = ftlGeometry.getPosition(id).z();
            for(auto& ecal_hit : pho.superCluster()->hitsAndFractions())
            {
                //---match FTL hits to ECAL hits
                //   (selection is loose: 0.0175 is the size of a barrel ECAL crystal in eta/phi)
                const CaloCellGeometry *cellGeometry = ecalBarrelGeometry->getGeometry(EBDetId(ecal_hit.first));
                float ecal_hit_eta = cellGeometry->getPosition().eta();
                float ecal_hit_phi = cellGeometry->getPosition().phi();
                if(fabs(ecal_hit_eta-ftl_hit_eta) < cellGeometry->etaSpan() &&
                   fabs(deltaPhi(ecal_hit_phi, ftl_hit_phi)) < cellGeometry->phiSpan())
                {
                    outTree_.ftlHitsPhoIdx->at(idx).push_back(idx);
                    outTree_.ftlHitsEnergy->at(idx).push_back(ftl_hit.energy());
                    outTree_.ftlHitsTime->at(idx).push_back(ftl_hit.time());
                    outTree_.ftlHitsEta->at(idx).push_back(ftl_hit_eta);
                    outTree_.ftlHitsPhi->at(idx).push_back(ftl_hit_phi);
                    outTree_.ftlHitsZ->at(idx).push_back(ftl_hit_z);                    
                    outTree_.ftlHitsEnergySum->at(idx)+= ftl_hit.energy();
                    outTree_.ftlNHits->at(idx)++;

                    float seedEta = EBDetId(pho.superCluster()->seed()->seed()).approxEta();
                    float seedPhi = EBDetId(pho.superCluster()->seed()->seed()).iphi();
                    seedPhi = seedPhi<=180 ? seedPhi/180.*TMath::Pi() : (seedPhi-360)/180.*TMath::Pi();
                    if(fabs(seedEta-ftl_hit_eta) < 0.0175*1.5 &&
                       fabs(deltaPhi(seedPhi, ftl_hit_phi)) < 0.0175*1.5)
                    {
                        outTree_.ftlHits3x3Sum->at(idx)+= ftl_hit.energy();
                        outTree_.ftlNHits3x3->at(idx)++;
                    }

                    //---if ftl hit matches ECAL hit break to avoid double counting
                    break;
                }
            }
        }
        
        ++idx;
    }
    if(photons.size() >= 2)
        outTree_.mass = (photons[0].p4() + photons[1].p4()).mass();
    else
        outTree_.mass = -1;

    outTree_.GetTTreePtr()->Fill();
}

typedef FTLDumpPhotons<std::vector<reco::Photon> > FTLDumpPhotonsRECO;
typedef FTLDumpPhotons<pat::PhotonCollection> FTLDumpPhotonsPAT;

DEFINE_FWK_MODULE(FTLDumpPhotonsRECO);
DEFINE_FWK_MODULE(FTLDumpPhotonsPAT);

#endif
