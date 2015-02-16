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
    float maxE_time;
    float maxE_energy;
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
    maxE_time=0;
    maxE_energy=0;

    //---create branches
    tree_->Branch("event", &event_n, "event/I");
    tree_->Branch("particle_n", &particle_n, "particle_n/I");
    tree_->Branch("particle_type", &particle_type, "particle_type/I");
    tree_->Branch("maxE_time", &maxE_time, "maxE_time/F");
    tree_->Branch("maxE_energy", &maxE_energy, "maxE_energy/F");
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
}

#endif 
