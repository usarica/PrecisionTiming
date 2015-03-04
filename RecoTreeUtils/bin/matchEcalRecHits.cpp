#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "FastTiming/RecoTreeUtils/interface/PFCandidateWithFT.h"
#include "FastTiming/RecoTreeUtils/interface/FTFile.h"

using namespace std;

//****************************************************************************************

int main(int argc, char* argv[])
{
    gSystem->Load("libFWCoreFWLite.so"); 
    AutoLibraryLoader::enable();

    if(argc < 2)
    {
        cout << "Usage : " << argv[0] << " [parameters.py]" << endl;
        return 0;
    }
    if(!edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process"))
    {
        cout << " ERROR: ParametersSet 'process' is missing in your configuration file"
             << endl;
        return 0;
    }
    //---get the python configuration---
    const edm::ParameterSet& process = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
    const edm::ParameterSet& filesOpt = process.getParameter<edm::ParameterSet>("ioFiles");
    //---io files option---
    vector<string> filesList = filesOpt.getParameter<vector<string> >("inputFiles");    
    TFile* outFile = TFile::Open(filesOpt.getParameter<string>("outputFile").c_str(), "recreate");
    outFile->cd();
    FTParticlesTree outTree;

    int iEvent=0;
    for(unsigned int iFile=0; iFile<filesList.size(); iFile++)
    {
        TFile* inFile = TFile::Open(filesList.at(iFile).c_str());
	std::cout << " >>> " << filesList.at(iFile) << std::endl;
        //---FWLite interfaces---
        fwlite::Event event(inFile);
        fwlite::Handle<vector<SimVertex> > genVtxHandle;
        fwlite::Handle<vector<reco::PFCandidate> > candHandle;
        // fwlite::Handle<vector<reco::PFJet> > jetsHandle;
        // fwlite::Handle<vector<reco::GenJet> > genJetsHandle;
        fwlite::Handle<edm::SortedCollection<EcalRecHit, 
                                             edm::StrictWeakOrdering<EcalRecHit > > > recSort;    
        //---events loop---
        for(event.toBegin(); !event.atEnd(); ++event)
        {
            outTree.event_n = iEvent;
            cout << "\r### EVENT: " << iEvent << flush;
            iEvent++;            
            //---get gen vertex time---
            const SimVertex* primaryVtx=NULL;
            genVtxHandle.getByLabel(event, "g4SimHits");
            if(genVtxHandle.ptr()->size() == 0 || genVtxHandle.ptr()->at(0).vertexId() != 0)
                continue;
            primaryVtx = &genVtxHandle.ptr()->at(0);
            //---fill gen vtx infos 
            outTree.gen_vtx_z = primaryVtx->position().z();
            outTree.gen_vtx_t = primaryVtx->position().t()*1E9;                
            //---get EK detailed time RecHits---
            recSort.getByLabel(event, "ecalDetailedTimeRecHit", "EcalRecHitsEK", "RECO");
            if(!recSort.isValid())
                continue;
            vector<EcalRecHit>* recVect = (vector<EcalRecHit>*)recSort.ptr();
            //---loop over all particles---
            candHandle.getByLabel(event, "particleFlow");
            for(unsigned int iCand=0; iCand<candHandle.ptr()->size(); iCand++)
            {                
                PFCandidateWithFT particle(&candHandle.ptr()->at(iCand), recVect,
                                           primaryVtx);
                if(particle.particleId() > 4 || !particle.GetPFCluster())
                    continue;
                outTree.particle_n = iCand;
                outTree.particle_type = particle.particleId();
                outTree.particle_E = particle.energy();
                outTree.particle_pt = particle.pt();
                outTree.particle_eta = particle.eta();
                outTree.particle_phi = particle.phi();
                outTree.maxE_time = particle.GetRecHitTimeMaxE().first;
                outTree.maxE_energy = particle.GetRecHitTimeMaxE().second;                
                outTree.all_time.clear();
                outTree.all_energy.clear();
                outTree.track_length = particle.GetTrackLength();
                outTree.trackCluster_dr = particle.GetDrTrackCluster();                
                vector<pair<float, float> > TandE = particle.GetRecHitsTimeE();
                if(TandE.size() == 0)
                    continue;
                for(unsigned int iRec=0; iRec<TandE.size(); iRec++)
                {
                    outTree.all_time.push_back(TandE.at(iRec).first);
                    outTree.all_energy.push_back(TandE.at(iRec).second);
                }
                outTree.Fill();
            }
        }
    }
    outFile->cd();
    outTree.Write("fast_timing");
    outFile->Close();
    return 0;
}
