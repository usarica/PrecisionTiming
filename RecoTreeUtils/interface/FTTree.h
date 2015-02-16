#ifndef FTTREE_H
#define FTTREE_H

#include "TTree.h"

class FTTree:

public: 

outTree->Branch("event", &event_n, "event/I");
    outTree->Branch("particle_n", &particle_n, "particle_n/I");
    outTree->Branch("particle_type", &particle_type, "particle_type/I");
    outTree->Branch("maxE_time", &maxE_time, "maxE_time/F");
    outTree->Branch("maxE_energy", &maxE_energy, "maxE_energy/F");
    outTree->Branch("all_time", "std::vector<float>", &all_time);
    outTree->Branch("all_energy", "std::vector<float>", &all_energy);

#endif FTTREE_H
