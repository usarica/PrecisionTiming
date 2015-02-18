#ifndef FTTREE_H
#define FTTREE_H

#include <string>
#include <vector>

#include "TTree.h"

using namespace std;

class FTTree
{
public: 

    //---ctors---
    FTTree();
    FTTree(TTree* tree);
    //---dtor---
    ~FTTree() {};
    //---wrappers
    inline void Fill() {tree_->Fill();};
    inline void Write(const char* name) {tree_->Write(name);};
    inline void Write(string name) {tree_->Write(name.c_str());};
    
    //---branches variables---
    int event_n;
    int particle_n;
    int particle_type;
    float particle_E;
    float particle_pt;
    float particle_eta;
    float particle_phi;
    float maxE_time;
    float maxE_energy;
    float track_phiIn;
    float track_phiOut;
    int track_charge;
    float trackCluster_dr;
    vector<float> all_time;
    vector<float> all_energy;


private:

    TTree* tree_;
    
};

FTTree::FTTree()
{
    tree_ = new TTree();
    //---init
    event_n=0;
    particle_n=0;
    particle_type=0;
    particle_E=0;
    particle_pt=0;
    particle_eta=0;
    particle_phi=0;
    maxE_time=0;
    maxE_energy=0;
    track_phiIn = 0.;
    track_phiOut = 0.;
    track_charge = 0;
    trackCluster_dr = 0.;

    //---create branches
    tree_->Branch("event", &event_n, "event/I");
    tree_->Branch("particle_n", &particle_n, "particle_n/I");
    tree_->Branch("particle_type", &particle_type, "particle_type/I");
    tree_->Branch("particle_E", &particle_E, "particle_E/F");
    tree_->Branch("particle_pt", &particle_pt, "particle_pt/F");
    tree_->Branch("particle_eta", &particle_eta, "particle_eta/F");
    tree_->Branch("particle_phi", &particle_phi, "particle_phi/F");
    tree_->Branch("maxE_time", &maxE_time, "maxE_time/F");
    tree_->Branch("maxE_energy", &maxE_energy, "maxE_energy/F");
    tree_->Branch("track_phiIn", &track_phiIn, "track_phiIn/F");
    tree_->Branch("track_phiOut", &track_phiOut, "track_phiOut/F");
    tree_->Branch("track_charge", &track_charge, "track_charge/I");
    tree_->Branch("trackCluster_dr", &trackCluster_dr, "trackCluster_dr/F");
    tree_->Branch("all_time", "std::vector<float>", &all_time);
    tree_->Branch("all_energy", "std::vector<float>", &all_energy);
}

FTTree::FTTree(TTree* tree)
{
    tree_ = tree;
    event_n=0;
    particle_n=0;
    particle_type=0;
    maxE_time=0;
    maxE_energy=0;
    track_phiIn = 0.;
    track_phiOut = 0.;
    track_charge = 0;
    trackCluster_dr = 0.;
}

#endif 
