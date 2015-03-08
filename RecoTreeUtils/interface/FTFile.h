#ifndef FTFILE_H
#define FTFILE_H

#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using namespace std;

class FTParticlesTree
{
public: 

    //---ctors---
    FTParticlesTree();
    //---dtor---
    ~FTParticlesTree() {};
    //---wrappers
    inline void Fill() {tree_->Fill();};
    inline void Write(const char* name) {tree_->Write(name);};
    inline void Write(string name) {tree_->Write(name.c_str());};
    
    //---branches variables---
    int event_n;
    float gen_vtx_z;
    float gen_vtx_t;
    int particle_n;
    int particle_type;
    float particle_E;
    float particle_pt;
    float particle_eta;
    float particle_phi;
    float maxE_time;
    float maxE_energy;
    int reco_vtx_index;
    float reco_vtx_z;
    float reco_vtx_t;
    float track_length;
    float track_length_helix;
    float track_dz;
    float track_dxy;
    float trackCluster_dr;
    vector<float> all_time;
    vector<float> all_energy;

private:

    TTree* tree_;
    
};

FTParticlesTree::FTParticlesTree()
{
    tree_ = new TTree();
    //---init
    event_n=0;
    gen_vtx_z=0;
    gen_vtx_t=0;
    particle_n=0;
    particle_type=0;
    particle_E=0;
    particle_pt=0;
    particle_eta=0;
    particle_phi=0;
    maxE_time=0;
    maxE_energy=0;
    reco_vtx_index=0;
    reco_vtx_z=0;
    reco_vtx_t=0;
    track_length=0;
    track_length_helix=0;
    track_dz=0;
    track_dxy=0;
    trackCluster_dr=0;

    //---create branches
    tree_->Branch("event", &event_n, "event/I");
    tree_->Branch("gen_vtx_z", &gen_vtx_z, "gen_vtx_z/F");
    tree_->Branch("gen_vtx_t", &gen_vtx_t, "gen_vtx_t/F");
    tree_->Branch("particle_n", &particle_n, "particle_n/I");
    tree_->Branch("particle_type", &particle_type, "particle_type/I");
    tree_->Branch("particle_E", &particle_E, "particle_E/F");
    tree_->Branch("particle_pt", &particle_pt, "particle_pt/F");
    tree_->Branch("particle_eta", &particle_eta, "particle_eta/F");
    tree_->Branch("particle_phi", &particle_phi, "particle_phi/F");
    tree_->Branch("maxE_time", &maxE_time, "maxE_time/F");
    tree_->Branch("maxE_energy", &maxE_energy, "maxE_energy/F");
    tree_->Branch("reco_vtx_index", &reco_vtx_index, "reco_vtx_index/I");
    tree_->Branch("reco_vtx_z", &reco_vtx_z, "reco_vtx_z/F");
    tree_->Branch("reco_vtx_t", &reco_vtx_t, "reco_vtx_t/F");
    tree_->Branch("track_length", &track_length, "track_length/F");
    tree_->Branch("track_length_helix", &track_length_helix, "track_length_helix/F");
    tree_->Branch("track_dz", &track_dz, "track_dz/F");
    tree_->Branch("track_dxy", &track_dxy, "track_dxy/F");
    tree_->Branch("trackCluster_dr", &trackCluster_dr, "trackCluster_dr/F");
    tree_->Branch("all_time", "std::vector<float>", &all_time);
    tree_->Branch("all_energy", "std::vector<float>", &all_energy);
}

class FTVerticesTree
{
public: 

    //---ctors---
    FTVerticesTree();
    //---dtor---
    ~FTVerticesTree() {};
    //---wrappers
    inline void Fill() {tree_->Fill();};
    inline void Write(const char* name) {tree_->Write(name);};
    inline void Write(string name) {tree_->Write(name.c_str());};
    
    //---branches variables---
    int event_n;
    float gen_vtx_z;
    float gen_vtx_t;
    int reco_vtx_index;
    int reco_vtx_n_part;
    float reco_vtx_ndof;
    float reco_vtx_chi2;
    float reco_vtx_sumpt2;
    float reco_vtx_seed_pt;
    float reco_vtx_seed_t;
    float reco_vtx_z;
    float reco_vtx_t;

private:

    TTree* tree_;
    
};

FTVerticesTree::FTVerticesTree()
{
    tree_ = new TTree();
    //---init
    event_n=0;
    gen_vtx_z=0;
    gen_vtx_t=0;
    reco_vtx_index=0;
    reco_vtx_n_part=0;
    reco_vtx_ndof=0;
    reco_vtx_chi2=0;
    reco_vtx_sumpt2=0;
    reco_vtx_seed_pt=0;
    reco_vtx_seed_t=0;
    reco_vtx_z=0;
    reco_vtx_t=0;

    //---create branches
    tree_->Branch("event", &event_n, "event/I");
    tree_->Branch("gen_vtx_z", &gen_vtx_z, "gen_vtx_z/F");
    tree_->Branch("gen_vtx_t", &gen_vtx_t, "gen_vtx_t/F");
    tree_->Branch("reco_vtx_index", &reco_vtx_index, "reco_vtx_index/I");
    tree_->Branch("reco_vtx_n_part", &reco_vtx_n_part, "reco_vtx_n_part/I");
    tree_->Branch("reco_vtx_ndof", &reco_vtx_ndof, "reco_vtx_ndof/F");
    tree_->Branch("reco_vtx_chi2", &reco_vtx_chi2, "reco_vtx_chi2/F");
    tree_->Branch("reco_vtx_sumpt2", &reco_vtx_sumpt2, "reco_vtx_sumpt2/F");
    tree_->Branch("reco_vtx_seed_pt", &reco_vtx_seed_pt, "reco_vtx_seed_pt/F");
    tree_->Branch("reco_vtx_seed_t", &reco_vtx_seed_t, "reco_vtx_seed_t/F");
    tree_->Branch("reco_vtx_z", &reco_vtx_z, "reco_vtx_z/F");
    tree_->Branch("reco_vtx_t", &reco_vtx_t, "reco_vtx_t/F");
}

class FTFile
{
public:    
    
    FTFile();
    FTFile(TFile* file);

    inline void Close() {file_->Close();};
    inline void cd() {file_->cd();};

    FTParticlesTree particlesTree;
    FTVerticesTree  verticesTree;
    
private:
    
    TFile*          file_;
};

FTFile::FTFile()
{}

FTFile::FTFile(TFile* file)
{
    file_ = file;
    file_->cd();
}

#endif 
