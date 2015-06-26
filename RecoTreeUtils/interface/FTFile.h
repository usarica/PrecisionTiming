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
    float particle_p;
    float particle_pz;
    float particle_pt;
    float particle_eta;
    float particle_phi;
    float cluster_E;
    float cluster_eta;
    float cluster_phi;
    float maxE_time;
    float maxE_energy;
    int reco_vtx_index;
    float reco_vtx_z;
    float reco_vtx_t;
    float sumEt_nocut;
    float sumEt_t_cut;
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
    particle_p=0;
    particle_pz=0;
    particle_pt=0;
    particle_eta=0;
    particle_phi=0;
    cluster_E=0;
    cluster_eta=0;
    cluster_phi=0;
    maxE_time=0;
    maxE_energy=0;
    reco_vtx_index=0;
    reco_vtx_z=0;
    reco_vtx_t=0;
    sumEt_nocut=0;
    sumEt_t_cut=0;
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
    tree_->Branch("particle_p", &particle_p, "particle_p/F");
    tree_->Branch("particle_pz", &particle_pz, "particle_pz/F");
    tree_->Branch("particle_pt", &particle_pt, "particle_pt/F");
    tree_->Branch("particle_eta", &particle_eta, "particle_eta/F");
    tree_->Branch("particle_phi", &particle_phi, "particle_phi/F");
    tree_->Branch("cluster_E", &cluster_E, "cluster_E/F");
    tree_->Branch("cluster_eta", &cluster_eta, "cluster_eta/F");
    tree_->Branch("cluster_phi", &cluster_phi, "cluster_phi/F");
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
    int gen_vtx_id;
    float gen_vtx_x;
    float gen_vtx_y;
    float gen_vtx_z;
    float gen_vtx_t;
    int reco_vtx_index;
    int reco_vtx_n_part;
    int reco_vtx_n_cha;
    int reco_vtx_n_neu;
    int reco_vtx_n_part_EE;
    int reco_vtx_n_cha_EE;
    int reco_vtx_n_neu_EE;
    float reco_vtx_ndof;
    float reco_vtx_chi2;
    float reco_vtx_sumpt2;
    float reco_vtx_sumpt2_EE;
    float reco_vtx_seed_pt;
    float reco_vtx_seed_eta;
    float reco_vtx_seed_t;
    float reco_vtx_z;
    float reco_vtx_t;
    float reco_vtx_cha_t;
    float reco_vtx_neu_t;
    float reco_vtx_t_EE;
    float reco_vtx_cha_t_EE;
    float reco_vtx_neu_t_EE;

private:

    TTree* tree_;
    
};

FTVerticesTree::FTVerticesTree()
{
    tree_ = new TTree();
    //---init
    event_n=0;    
    gen_vtx_id=-1;
    gen_vtx_x=0;
    gen_vtx_y=0;
    gen_vtx_z=0;
    gen_vtx_t=0;
    reco_vtx_index=0;
    reco_vtx_n_part=0;
    reco_vtx_n_cha=0;
    reco_vtx_n_neu=0;
    reco_vtx_n_part_EE=0;
    reco_vtx_n_cha_EE=0;
    reco_vtx_n_neu_EE=0;
    reco_vtx_ndof=0;
    reco_vtx_chi2=0;
    reco_vtx_sumpt2=0;
    reco_vtx_sumpt2_EE=0;
    reco_vtx_seed_pt=0;
    reco_vtx_seed_eta=0;
    reco_vtx_seed_t=0;
    reco_vtx_z=0;
    reco_vtx_t=0;
    reco_vtx_cha_t=0;
    reco_vtx_neu_t=0;
    reco_vtx_t_EE=0;
    reco_vtx_cha_t_EE=0;
    reco_vtx_neu_t_EE=0;

    //---create branches
    tree_->Branch("event", &event_n, "event/I");
    tree_->Branch("gen_vtx_id", &gen_vtx_id, "gen_vtx_id/I");
    tree_->Branch("gen_vtx_x", &gen_vtx_x, "gen_vtx_x/F");
    tree_->Branch("gen_vtx_y", &gen_vtx_y, "gen_vtx_y/F");
    tree_->Branch("gen_vtx_z", &gen_vtx_z, "gen_vtx_z/F");
    tree_->Branch("gen_vtx_t", &gen_vtx_t, "gen_vtx_t/F");
    tree_->Branch("reco_vtx_index", &reco_vtx_index, "reco_vtx_index/I");
    tree_->Branch("reco_vtx_n_part", &reco_vtx_n_part, "reco_vtx_n_part/I");
    tree_->Branch("reco_vtx_n_cha", &reco_vtx_n_cha, "reco_vtx_n_cha/I");
    tree_->Branch("reco_vtx_n_neu", &reco_vtx_n_neu, "reco_vtx_n_neu/I");
    tree_->Branch("reco_vtx_n_part_EE", &reco_vtx_n_part_EE, "reco_vtx_n_part_EE/I");
    tree_->Branch("reco_vtx_n_cha_EE", &reco_vtx_n_cha_EE, "reco_vtx_n_cha_EE/I");
    tree_->Branch("reco_vtx_n_neu_EE", &reco_vtx_n_neu_EE, "reco_vtx_n_neu_EE/I");
    tree_->Branch("reco_vtx_ndof", &reco_vtx_ndof, "reco_vtx_ndof/F");
    tree_->Branch("reco_vtx_chi2", &reco_vtx_chi2, "reco_vtx_chi2/F");
    tree_->Branch("reco_vtx_sumpt2", &reco_vtx_sumpt2, "reco_vtx_sumpt2/F");
    tree_->Branch("reco_vtx_sumpt2_EE", &reco_vtx_sumpt2_EE, "reco_vtx_sumpt2_EE/F");
    tree_->Branch("reco_vtx_seed_pt", &reco_vtx_seed_pt, "reco_vtx_seed_pt/F");
    tree_->Branch("reco_vtx_seed_eta", &reco_vtx_seed_eta, "reco_vtx_seed_eta/F");
    tree_->Branch("reco_vtx_seed_t", &reco_vtx_seed_t, "reco_vtx_seed_t/F");
    tree_->Branch("reco_vtx_z", &reco_vtx_z, "reco_vtx_z/F");
    tree_->Branch("reco_vtx_t", &reco_vtx_t, "reco_vtx_t/F");
    tree_->Branch("reco_vtx_cha_t", &reco_vtx_cha_t, "reco_vtx_cha_t/F");
    tree_->Branch("reco_vtx_neu_t", &reco_vtx_neu_t, "reco_vtx_neu_t/F");
    tree_->Branch("reco_vtx_t_EE", &reco_vtx_t, "reco_vtx_t/F");
    tree_->Branch("reco_vtx_cha_t_EE", &reco_vtx_cha_t_EE, "reco_vtx_cha_t_EE/F");
    tree_->Branch("reco_vtx_neu_t_EE", &reco_vtx_neu_t_EE, "reco_vtx_neu_t_EE/F");
}

class FTGlobalTree
{
public: 

    //---ctors---
    FTGlobalTree();
    //---dtor---
    ~FTGlobalTree() {};
    //---wrappers
    inline void Fill() {tree_->Fill();};    
    inline void Write(const char* name) {tree_->Write(name);};
    inline void Write(string name) {tree_->Write(name.c_str());};
    void        Reset();
    
    //---branches variables---
    float sumEt_nocut;
    float sumEt_t_cut[200];
    float sumEt_gen;
    int   nEEplus[200];
    int   nEEminus[200];
    int   n_vtx;
    int   vtx_id[200];
    
private:

    TTree* tree_;
    
};

FTGlobalTree::FTGlobalTree()
{
    tree_ = new TTree();

    //---init
    Reset();
    
    //---create branches
    tree_->Branch("sumEt_nocut", &sumEt_nocut, "sumEt_nocut/F");
    tree_->Branch("sumEt_t_cut", &sumEt_t_cut, "sumEt_t_cut[200]/F");
    tree_->Branch("sumEt_gen", &sumEt_gen, "sumEt_gen/F");
    tree_->Branch("nEEplus", &nEEplus, "nEEplus[200]/I");
    tree_->Branch("nEEminus", &nEEminus, "nEEminus[200]/I");
    tree_->Branch("n_vtx", &n_vtx, "n_vtx/I");
    tree_->Branch("vtx_id", &vtx_id, "vtx_id[200]/I");
}

void FTGlobalTree::Reset()
{
    //---reset
    sumEt_nocut=0;
    sumEt_gen=0;
    for(unsigned int i=0; i<200; ++i)
    {
        sumEt_t_cut[i]=0;
        nEEplus[i]=0;
        nEEminus[i]=0;
        vtx_id[i]=0;
    }
}

class FTJetsTree
{
public: 

    //---ctors---
    FTJetsTree();
    //---dtor---
    ~FTJetsTree() {};
    //---wrappers
    inline void Fill() {tree_->Fill();};
    inline void Write(const char* name) {tree_->Write(name);};
    inline void Write(string name) {tree_->Write(name.c_str());};
    
    //---branches variables---
    int gen_n;
    float gen_j1_pt;
    float gen_j2_pt;
    float gen_j1_eta;
    float gen_j2_eta;
    float gen_j1_E;
    float gen_j2_E;
    float gen_jj_m;
    int chs_n;
    float chs_j1_pt;
    float chs_j2_pt;
    float chs_j1_eta;
    float chs_j2_eta;
    float chs_j1_E;
    float chs_j2_E;
    float chs_j1_dR;
    float chs_j2_dR;
    float chs_jj_m;
    int tcut_n;
    float tcut_j1_pt;
    float tcut_j2_pt;
    float tcut_j1_eta;
    float tcut_j2_eta;

    float tcut_j1_pVtx_SeedEta;
    float tcut_j1_pVtx_NPart;
    float tcut_j1_pVtx_NPartEE;
    float tcut_j2_pVtx_SeedEta;
    float tcut_j2_pVtx_NPart;
    float tcut_j2_pVtx_NPartEE;

    float tcut_j1_E;
    float tcut_j2_E;
    float tcut_j1_dR;
    float tcut_j2_dR;
    float tcut_jj_m;
    
private:

    TTree* tree_;
    
};

FTJetsTree::FTJetsTree()
{
    tree_ = new TTree();
    
    //---init
    gen_n=0;
    gen_j1_pt=0;
    gen_j2_pt=0;
    gen_j1_eta=0;
    gen_j2_eta=0;
    gen_j1_E=0;
    gen_j2_E=0;
    gen_jj_m=0;
    chs_n=0;
    chs_j1_pt=0;
    chs_j2_pt=0;
    chs_j1_eta=0;
    chs_j2_eta=0;
    chs_j1_E=0;
    chs_j2_E=0;
    chs_j1_dR=0;
    chs_j2_dR=0;
    chs_jj_m=0;
    tcut_n=0;
    tcut_j1_pt=0;
    tcut_j2_pt=0;
    tcut_j1_eta=0;
    tcut_j2_eta=0;

    tcut_j1_pVtx_SeedEta = 0;
    tcut_j1_pVtx_NPart = 0;
    tcut_j1_pVtx_NPartEE = 0;
    tcut_j2_pVtx_SeedEta = 0;
    tcut_j2_pVtx_NPart = 0;
    tcut_j2_pVtx_NPartEE = 0;

    tcut_j1_E=0;
    tcut_j2_E=0;
    tcut_j1_dR=0;
    tcut_j2_dR=0;
    tcut_jj_m=0;

    //---create branches
    tree_->Branch("gen_n", &gen_n, "gen_n/I");
    tree_->Branch("gen_j1_pt", &gen_j1_pt, "gen_j1_pt/F");
    tree_->Branch("gen_j2_pt", &gen_j2_pt, "gen_j2_pt/F");    
    tree_->Branch("gen_j1_eta", &gen_j1_eta, "gen_j1_eta/F");
    tree_->Branch("gen_j2_eta", &gen_j2_eta, "gen_j2_eta/F");
    tree_->Branch("gen_j1_E", &gen_j1_E, "gen_j1_E/F");
    tree_->Branch("gen_j2_E", &gen_j2_E, "gen_j2_E/F");
    tree_->Branch("gen_jj_m", &gen_jj_m, "gen_jj_m/F");
    tree_->Branch("chs_n", &chs_n, "chs_n/I");
    tree_->Branch("chs_j1_pt", &chs_j1_pt, "chs_j1_pt/F");
    tree_->Branch("chs_j2_pt", &chs_j2_pt, "chs_j2_pt/F");    
    tree_->Branch("chs_j1_eta", &chs_j1_eta, "chs_j1_eta/F");
    tree_->Branch("chs_j2_eta", &chs_j2_eta, "chs_j2_eta/F");
    tree_->Branch("chs_j1_E", &chs_j1_E, "chs_j1_E/F");
    tree_->Branch("chs_j2_E", &chs_j2_E, "chs_j2_E/F");
    tree_->Branch("chs_j1_dR", &chs_j1_dR, "chs_j1_dR/F");
    tree_->Branch("chs_j2_dR", &chs_j2_dR, "chs_j2_dR/F");
    tree_->Branch("chs_jj_m", &chs_jj_m, "chs_jj_m/F");
    tree_->Branch("tcut_n", &tcut_n, "tcut_n/I");
    tree_->Branch("tcut_j1_pt", &tcut_j1_pt, "tcut_j1_pt/F");
    tree_->Branch("tcut_j2_pt", &tcut_j2_pt, "tcut_j2_pt/F");    
    tree_->Branch("tcut_j1_eta", &tcut_j1_eta, "tcut_j1_eta/F");
    tree_->Branch("tcut_j2_eta", &tcut_j2_eta, "tcut_j2_eta/F");
    tree_->Branch("tcut_j1_pVtx_SeedEta", &tcut_j1_pVtx_SeedEta, "tcut_j1_pVtx_SeedEta/F");
    tree_->Branch("tcut_j1_pVtx_NPart", &tcut_j1_pVtx_NPart, "tcut_j1_pVtx_NPart/I");
    tree_->Branch("tcut_j1_pVtx_NPartEE", &tcut_j1_pVtx_NPartEE, "tcut_j1_pVtx_NPartEE/I");
    tree_->Branch("tcut_j2_pVtx_SeedEta", &tcut_j2_pVtx_SeedEta, "tcut_j2_pVtx_SeedEta/F");
    tree_->Branch("tcut_j2_pVtx_NPart", &tcut_j2_pVtx_NPart, "tcut_j2_pVtx_NPart/I");
    tree_->Branch("tcut_j2_pVtx_NPartEE", &tcut_j2_pVtx_NPartEE, "tcut_j2_pVtx_NPartEE/I");
    tree_->Branch("tcut_j1_E", &tcut_j1_E, "tcut_j1_E/F");
    tree_->Branch("tcut_j2_E", &tcut_j2_E, "tcut_j2_E/F");
    tree_->Branch("tcut_j1_dR", &tcut_j1_dR, "tcut_j1_dR/F");
    tree_->Branch("tcut_j2_dR", &tcut_j2_dR, "tcut_j2_dR/F");
    tree_->Branch("tcut_jj_m", &tcut_jj_m, "tcut_jj_m/F");                                                                                                                                                                           
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
    FTGlobalTree    globalTree;
    FTJetsTree      jetsTree;
    
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
