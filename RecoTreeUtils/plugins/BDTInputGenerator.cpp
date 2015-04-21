#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/EcalDetId/interface/EKDetId.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"

using namespace std;

struct FTEcalRecHit {
    int    ix;
    int    iy;
    double z;
    double time;
    double energy;
};

class BDTInputGenerator : public edm::EDAnalyzer
{
public:
    explicit BDTInputGenerator(const edm::ParameterSet&);
    ~BDTInputGenerator() {};

    //---utils---
    void BuildRecHitsMatrix(vector<EcalRecHit*> recHits, DetId seed, int sqrt_n);
        
private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    const CaloGeometry* skGeometry_;
    const MagneticField* magField_;
    //---output file---
    edm::Service<TFileService> fs_;
    TFile* outFile_;
    TTree* outTree_;    
    float* times_;
    float* energies_;
    float* pos_x_;
    float* pos_y_;
    float true_time_;
    //---objects interfaces---
    edm::ESHandle<MagneticField> magFieldHandle_;             
    edm::ESHandle<CaloGeometry> geoHandle_;
    edm::Handle<vector<SimVertex> > genSigVtxHandle_;
    edm::Handle<vector<reco::Vertex> > recoVtxHandle_;
    edm::Handle<vector<reco::PFCandidate> > candHandle_;
    edm::Handle<edm::SortedCollection<EcalRecHit, 
                                      edm::StrictWeakOrdering<EcalRecHit > > > recSort_;
    //---FT objects---
    vector<VertexWithFT> recoVtxCollection_;
    vector<PFCandidateWithFT> particlesCollection_;
    map<int, FTEcalRecHit> recHitsMatrix_;
    //---Options---
    int particleType_;
    int matrixSide_;
};

BDTInputGenerator::BDTInputGenerator(const edm::ParameterSet& Config)
{
    particleType_ = Config.getUntrackedParameter<int>("particleType", 4);
    matrixSide_ = Config.getUntrackedParameter<int>("matrixSide", 5);
    times_ = new float[25]();
    energies_ = new float[25]();
    pos_x_ = new float[25]();
    pos_y_ = new float[25]();
    true_time_ = 0;
}

void BDTInputGenerator::beginJob()
{
    outFile_ = &fs_->file();
    outFile_->cd();
    outTree_ = new TTree("ft_bdt_input", "ft_bdt_input");
    for(int iRec=0; iRec<matrixSide_*matrixSide_; ++iRec)
    {
        TString timeBranch = TString::Format("time_xstal_%.2d", iRec);
        TString energyBranch = TString::Format("energy_xstal_%.2d", iRec);
        TString posxBranch = TString::Format("posx_xstal_%.2d", iRec);
        TString posyBranch = TString::Format("posy_xstal_%.2d", iRec);
        outTree_->Branch(timeBranch.Data(), &times_[iRec], (timeBranch+"/F").Data());
        outTree_->Branch(energyBranch.Data(), &energies_[iRec], (energyBranch+"/F").Data());
        outTree_->Branch(posxBranch.Data(), &pos_x_[iRec], (posxBranch+"/F").Data());
        outTree_->Branch(posyBranch.Data(), &pos_y_[iRec], (posyBranch+"/F").Data());
    }
    outTree_->Branch("true_time", &true_time_, "true_time/F");
}

void BDTInputGenerator::endJob()
{
    outFile_->cd();
    outTree_->Write("ft_bdt_input");
}

void BDTInputGenerator::analyze(const edm::Event& Event, const edm::EventSetup& Setup)
{
    //---get the magnetic field---
    Setup.get<IdealMagneticFieldRecord>().get(magFieldHandle_);
    magField_ = magFieldHandle_.product();
    //---get the geometry---
    Setup.get<CaloGeometryRecord>().get(geoHandle_);
    skGeometry_ = geoHandle_.product();
    //---get signal gen vertex time---
    const SimVertex* genSignalVtx=NULL;
    Event.getByLabel("g4SimHits", genSigVtxHandle_);
    if(genSigVtxHandle_.product()->size() == 0 || genSigVtxHandle_.product()->at(0).vertexId() != 0)
        return;
    genSignalVtx = &genSigVtxHandle_.product()->at(0);
    //---get reco primary vtxs---
    VertexWithFT* recoSignalVtx=NULL;
    //---get EK detailed time RecHits---
    Event.getByLabel(edm::InputTag("ecalDetailedTimeRecHit", "EcalRecHitsEK", "RECO"),
                      recSort_);
    if(!recSort_.isValid())
        return;
    vector<EcalRecHit>* recVect = (vector<EcalRecHit>*)recSort_.product();
    //---get all particles---
    Event.getByLabel("particleFlow", candHandle_);

    //---loop over all particles---
    for(unsigned int iCand=0; iCand<candHandle_.product()->size(); iCand++)
    {
        PFCandidateWithFT particle(&candHandle_.product()->at(iCand), recVect,
                                   skGeometry_, magField_, genSignalVtx, recoSignalVtx);

        if(particle.particleId() != particleType_ || !particle.hasTime())
            continue;
        BuildRecHitsMatrix(particle.GetRecHits(), particle.GetPFCluster()->seed(), matrixSide_);
        for(auto& mtxEle : recHitsMatrix_)
        {
            times_[mtxEle.first] = mtxEle.second.energy!=-1 ? (mtxEle.second.time) : -1;
            energies_[mtxEle.first] = mtxEle.second.energy!=-1 ? mtxEle.second.energy : -1;
            pos_x_[mtxEle.first] = mtxEle.second.ix;
            pos_y_[mtxEle.first] = mtxEle.second.iy;
        }
        float Dz = fabs(recHitsMatrix_[matrixSide_*matrixSide_/2].z - genSignalVtx->position().z());
        true_time_ = genSignalVtx->position().t()*1E9 + particle.p()*Dz/(fabs(particle.pz())*30);
        outTree_->Fill();
        return;
    }
}

void BDTInputGenerator::BuildRecHitsMatrix(vector<EcalRecHit*> recHits, DetId seed, int sqrt_n)
{
    recHitsMatrix_.clear();    

    int nRH=0;
    int seed_ix = EKDetId(seed).ix();
    int seed_iy = EKDetId(seed).iy();
    int* ixs = new int[sqrt_n*sqrt_n];
    int* iys = new int[sqrt_n*sqrt_n];
    for(int iX=-sqrt_n/2; iX<=sqrt_n/2; ++iX)
    {
        for(int iY=-sqrt_n/2; iY<=sqrt_n/2; ++iY)
        {
            ixs[nRH] = seed_ix+iX;
            iys[nRH] = seed_iy+iY;
            ++nRH;
        }
    }
    nRH=0;
    while(nRH < sqrt_n*sqrt_n)
    {
        for(auto& recHit : recHits)
        {
            if(EKDetId(recHit->id()).ix() == ixs[nRH] && EKDetId(recHit->id()).iy() == iys[nRH])
            {
                const CaloCellGeometry* cell=skGeometry_->getGeometry(recHit->id());
                GlobalPoint recHitPos = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0);
                FTEcalRecHit tmp = {ixs[nRH],
                                    iys[nRH],
                                    recHitPos.z(),
                                    recHit->time()*1E9+recHitPos.mag()/30,
                                    recHit->energy()};
                recHitsMatrix_[nRH] = tmp;
            }
        }
        if(recHitsMatrix_.find(nRH) == recHitsMatrix_.end())
            recHitsMatrix_[nRH] = FTEcalRecHit{ixs[nRH], iys[nRH], 0, 0,-1};
        ++nRH;                        
    }
}    

//define this as a plugin
DEFINE_FWK_MODULE(BDTInputGenerator);
