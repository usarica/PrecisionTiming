#ifndef _FTL_DUMPER_
#define _FTL_DUMPER_

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

#include "PrecisionTiming/FTLDumper/interface/FTLElectronsTree.h"

using namespace std;

class FTLDumper : public edm::EDAnalyzer
{
public:
    explicit FTLDumper(const edm::ParameterSet& pSet);
    ~FTLDumper() {};

    //---utils

    //---methods
    virtual void beginJob() override {};
    virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
    virtual void endJob() override {};


private:
    //---inputs
    edm::Handle<pat::ElectronCollection> electronsHandle_;
    edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;

    //---outputs
    FTLElectronsTree outTree_;
    edm::Service<TFileService> fs_;
};

FTLDumper::FTLDumper(const edm::ParameterSet& pSet):
    electronsToken_(consumes<pat::ElectronCollection>(pSet.getUntrackedParameter<edm::InputTag>("electronsTag")))
{
    outTree_ = FTLElectronsTree("t_ele", "Electrons tree for FTL studies");
}


void FTLDumper::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    outTree_.Reset();
    
    event.getByToken(electronsToken_, electronsHandle_);
    auto electrons = *electronsHandle_.product();

    int idx=0;
    for(auto& ele : electrons)
    {
        outTree_.idx->push_back(idx);
        outTree_.pt->push_back(ele.pt());
        outTree_.eta->push_back(ele.eta());
        outTree_.phi->push_back(ele.phi());
        outTree_.energy->push_back(ele.energy());
        outTree_.sc_energy->push_back(ele.superCluster()->energy());
        outTree_.r9->push_back(ele.r9());
        outTree_.sIeIe->push_back(ele.sigmaIetaIeta());
    }
    if(electrons.size() >= 2)
        outTree_.mass = (electrons[0].p4() + electrons[1].p4()).mass();
    else
        outTree_.mass = -1;

    outTree_.GetTTreePtr()->Fill();
}
    
DEFINE_FWK_MODULE(FTLDumper);

#endif
