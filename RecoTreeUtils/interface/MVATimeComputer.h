#ifndef _MVA_TIME_CUMPUTER_H_
#define _MVA_TIME_CUMPUTER_H_

#include<vector>

#include "TRandom.h"

#include "TMVA/Reader.h"     
#include "TMVA/Tools.h"      
#include "TMVA/MethodCuts.h"

using namespace std;

class FTEcalRecHit
{
public:
    //---ctors---
    FTEcalRecHit() {};
    FTEcalRecHit(unsigned int c_id, int c_ix, int c_iy, float c_z=0, float c_time=0, float c_energy=-1):
        id(c_id), ix(c_ix), iy(c_iy), z(c_z), time(c_time), energy(c_energy) {};

    //---dtor---
    ~FTEcalRecHit() {};

    //---utils---
    bool operator>(const FTEcalRecHit& other) const {return energy>other.energy;};
    bool operator<(const FTEcalRecHit& other) const {return energy<other.energy;};

    unsigned int id;
    int          ix;
    int          iy;
    float       z;
    float       time;
    float       energy;
};

class MVATimeComputer
{
public:
    //---ctors---
    MVATimeComputer();
    MVATimeComputer(string weightfile);

    //---dtor---
    ~MVATimeComputer() {};

    //---getters---
    float GetMVATime(vector<FTEcalRecHit>* recHits, float E, float pt, float pz, float smearing);

private:
    float* times_;
    float* E_ratios_;
    float  pt_;
    float  pz_;
    float  tan_theta_;
    float  smearing_;
    TMVA::Reader* reader;
};

#endif
