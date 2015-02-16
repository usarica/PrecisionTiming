#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "SimDataFormats/Vertex/interface/SimVertex.h"
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

#include "FastTiming/RecoTreeUtils/interface/ParticleWithFT.hpp"

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
    TTree* outTree = new TTree();
    int event_n=0, particle_n=0, particle_type=0;
    float maxE_time=0, maxE_energy=0;
    vector<float> all_time, all_energy;
    outTree->Branch("event", &event_n, "event/I");
    outTree->Branch("particle_n", &particle_n, "particle_n/I");
    outTree->Branch("particle_type", &particle_type, "particle_type/I");
    outTree->Branch("maxE_time", &maxE_time, "maxE_time/F");
    outTree->Branch("maxE_energy", &maxE_energy, "maxE_energy/F");
    outTree->Branch("all_time", "std::vector<float>", &all_time);
    outTree->Branch("all_energy", "std::vector<float>", &all_energy);

    int iEvent=0;
    for(unsigned int iFile=0; iFile<filesList.size(); iFile++)
    {
        TFile* inFile = TFile::Open(filesList.at(iFile).c_str());
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
            event_n = iEvent;
            cout << "\r### EVENT: " << iEvent << flush;
            iEvent++;            
            //---get gen vertex time---
            float primaryVtxTime=-1;
            genVtxHandle.getByLabel(event, "g4SimHits");
            if(genVtxHandle.ptr()->size() > 0 && genVtxHandle.ptr()->at(0).vertexId() == 0)
                primaryVtxTime = genVtxHandle.ptr()->at(0).position().t();
            //---get EK detailed time RecHits---
            recSort.getByLabel(event, "ecalDetailedTimeRecHit", "EcalRecHitsEK", "RECO");
            if(!recSort.isValid())
                continue;
            vector<EcalRecHit>* recVect = (vector<EcalRecHit>*)recSort.ptr();
            //---loop over all particles---
            candHandle.getByLabel(event, "particleFlow");
            for(unsigned int iCand=0; iCand<candHandle.ptr()->size(); iCand++)
            {                
                ParticleWithFT particle(&candHandle.ptr()->at(iCand), recVect, primaryVtxTime);                
                particle_n = iCand;
                particle_type = particle.Type();
                if(particle_type > 4)
                    continue;
                maxE_time = particle.GetRecHitTimeMaxE().first;
                maxE_energy = particle.GetRecHitTimeMaxE().second;
                all_time.clear();
                all_energy.clear();
                vector<pair<float, float> > TandE = particle.GetRecHitsTimeE();
                for(unsigned int iRec=0; iRec<TandE.size(); iRec++)
                {
                    all_time.push_back(TandE.at(iRec).first);
                    all_energy.push_back(TandE.at(iRec).second);
                }
                outTree->Fill();
            }
        }
    }
    outFile->cd();
    outTree->Write("fast_timing");
    outFile->Close();
    return 0;
}
