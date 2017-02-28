#ifndef _FTL_DUMP_ELECTRONS_
#define _FTL_DUMP_ELECTRONS_

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

#include "DataFormats/PatCandidates/interface/Electron.h"
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

#include "PrecisionTiming/FTLDumper/interface/FTLElectronsTree.h"

using namespace std;

template<class ElectronCollectionT>
class FTLDumpElectrons : public edm::EDAnalyzer
{
public:
    explicit FTLDumpElectrons(const edm::ParameterSet& pSet);
    ~FTLDumpElectrons() {};

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
    edm::Handle<ElectronCollectionT> electronsHandle_;
    edm::EDGetTokenT<ElectronCollectionT> electronsToken_;    
    
    //---options
    float mcTruthEleEtThr_;
    bool readFTLrecHits_;
    
    //---workers
    ElectronMCTruthFinder electronMCTruthFinder_;
    
    //---outputs
    FTLElectronsTree outTree_;
    edm::Service<TFileService> fs_;
};

template<class ElectronCollectionT>
FTLDumpElectrons<ElectronCollectionT>::FTLDumpElectrons(const edm::ParameterSet& pSet):
    simTkToken_(consumes<edm::SimTrackContainer>(pSet.getUntrackedParameter<edm::InputTag>("simTkTag"))),
    simVtxToken_(consumes<edm::SimVertexContainer>(pSet.getUntrackedParameter<edm::InputTag>("simVtxTag"))),
    ftlRecHitsToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("ftlRecHitsTag"))),
    electronsToken_(consumes<ElectronCollectionT>(pSet.getUntrackedParameter<edm::InputTag>("electronsTag"))),
    mcTruthEleEtThr_(pSet.getUntrackedParameter<double>("mcTruthEleEtThr")),
    readFTLrecHits_(pSet.getUntrackedParameter<bool>("readFTLRecHits")),
    electronMCTruthFinder_()    
{
    outTree_ = FTLElectronsTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "Electrons tree for FTL studies");
}

template<class ElectronCollectionT>
void FTLDumpElectrons<ElectronCollectionT>::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    outTree_.Reset();

    //---load the electrons
    event.getByToken(electronsToken_, electronsHandle_);
    auto electrons = *electronsHandle_.product();

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
    for(auto& ele : electrons)
    {
        //---electron standard info
        outTree_.idx->push_back(idx);
        outTree_.pt->push_back(ele.pt());
        outTree_.eta->push_back(ele.eta());
        outTree_.phi->push_back(ele.phi());
        outTree_.energy->push_back(ele.energy());
        outTree_.sc_energy->push_back(ele.superCluster()->energy());
        outTree_.sc_eop->push_back(ele.eSuperClusterOverP());        
        outTree_.r9->push_back(ele.r9());
        outTree_.sIeIe->push_back(ele.sigmaIetaIeta());
        outTree_.nBrem->push_back(ele.numberOfBrems());        
        outTree_.fBrem->push_back(ele.fbrem());
        
        //---MC truth conversion info (catch seg fault from EGammaTools)
        vector<ElectronMCTruth> mcTruthElectrons;
        //mcTruthElectrons = electronMCTruthFinder_.find(*simTkHandle_.product(), *simVtxHandle_.product());
        ElectronMCTruth mcTruthEle;
        bool mcTruthExist=false;
        float minDR=10;
        for(auto& mcele : mcTruthElectrons)
        {
            if(mcele.fourMomentum().et() > mcTruthEleEtThr_  &&
               deltaR(mcele.fourMomentum().eta(), mcele.fourMomentum().phi(), ele.eta(), ele.phi()) < minDR)
            {
                minDR = deltaR(mcele.fourMomentum().eta(), mcele.fourMomentum().phi(), ele.eta(), ele.phi());
                mcTruthEle = mcele;
                mcTruthExist = true;
            }
        }        
        if(mcTruthExist)
        {
            float tot_brem=0;
            for(auto& brem : mcTruthEle.eloss())
                tot_brem+=brem;
            outTree_.firstBremRadius->push_back(tot_brem > 0 ? mcTruthEle.bremVertices()[0].rho() : -1);
            outTree_.genBrem->push_back(tot_brem);
            outTree_.genNBrem->push_back(tot_brem > 0 ? mcTruthEle.eloss().size() : -1);            
            outTree_.genEnergy->push_back(mcTruthEle.fourMomentum().e());
            outTree_.genEta->push_back(mcTruthEle.fourMomentum().eta());
            outTree_.genPhi->push_back(mcTruthEle.fourMomentum().phi());
            outTree_.genPt->push_back(mcTruthEle.fourMomentum().et());
        }
        else
        {
            outTree_.firstBremRadius->push_back(-1);
            outTree_.genBrem->push_back(-1);
            outTree_.genNBrem->push_back(-1);            
            outTree_.genEnergy->push_back(-1);
            outTree_.genEta->push_back(-1);
            outTree_.genPhi->push_back(-1);
            outTree_.genPt->push_back(-1);
        }
        
        //---find ftl associated hits
        for(auto ftl_hit : ftlRecHits)
        {
            outTree_.ftlHitsEleIdx->resize(idx+1);
            outTree_.ftlHitsEnergy->resize(idx+1);
            outTree_.ftlHitsTime->resize(idx+1);
            outTree_.ftlHitsEta->resize(idx+1);
            outTree_.ftlHitsPhi->resize(idx+1);
            outTree_.ftlHitsEnergySum->resize(idx+1);
            outTree_.ftlNHits->resize(idx+1);
            outTree_.ftlHits3x3Sum->resize(idx+1);
            outTree_.ftlNHits3x3->resize(idx+1);
            FastTimeDetId id = ftl_hit.id();
            float hit_eta = -log(tan(atan(1189/(id.iz()*10-0.5))/2))*id.zside();
            float hit_phi = id.iphi()<=360 ? id.iphi()/360.*TMath::Pi() : (id.iphi()-720)/360.*TMath::Pi();
            if(fabs(ele.eta()-hit_eta) < ele.superCluster()->etaWidth() &&
               fabs(deltaPhi(ele.phi(), hit_phi)) < ele.superCluster()->phiWidth())
            {
                outTree_.ftlHitsEleIdx->at(idx).push_back(idx);
                outTree_.ftlHitsEnergy->at(idx).push_back(ftl_hit.energy());
                outTree_.ftlHitsTime->at(idx).push_back(ftl_hit.time());
                outTree_.ftlHitsEta->at(idx).push_back(hit_eta);
                outTree_.ftlHitsPhi->at(idx).push_back(hit_phi);
                outTree_.ftlHitsEnergySum->at(idx)+= ftl_hit.energy();
                outTree_.ftlNHits->at(idx)++;
            }
            float seedEta = EBDetId(ele.superCluster()->seed()->seed()).approxEta();
            float seedPhi = EBDetId(ele.superCluster()->seed()->seed()).iphi();
            seedPhi = seedPhi<=180 ? seedPhi/180.*TMath::Pi() : (seedPhi-360)/180.*TMath::Pi();
            if(fabs(seedEta-hit_eta) < 0.0175*1.5 &&
               fabs(deltaPhi(seedPhi, hit_phi)) < 0.0175*1.5)
            {
                outTree_.ftlHits3x3Sum->at(idx)+= ftl_hit.energy();
                outTree_.ftlNHits3x3->at(idx)++;
            }
        }
        
        ++idx;
    }
    if(electrons.size() >= 2)
        outTree_.mass = (electrons[0].p4() + electrons[1].p4()).mass();
    else
        outTree_.mass = -1;

    outTree_.GetTTreePtr()->Fill();
}

typedef FTLDumpElectrons<std::vector<reco::GsfElectron> > FTLDumpElectronsRECO;
typedef FTLDumpElectrons<pat::ElectronCollection> FTLDumpElectronsPAT;

DEFINE_FWK_MODULE(FTLDumpElectronsRECO);
DEFINE_FWK_MODULE(FTLDumpElectronsPAT);

#endif
