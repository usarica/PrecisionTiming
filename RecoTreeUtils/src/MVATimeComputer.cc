#include "FastTiming/RecoTreeUtils/interface/MVATimeComputer.h"

using namespace std;

//****************************************************************************************

//----------ctor with weights file--------------------------------------------------------
MVATimeComputer::MVATimeComputer(string weightfile)
{    
    times_ = new float[9]();
    E_ratios_ = new float[9]();
    pt_=0, pz_=0, tan_theta_=0;    

    //---construct TMVA reader
    reader = new TMVA::Reader("!Color:Silent");
    
    //---add variables
    for(int iRec=0; iRec<9; ++iRec)
    {          
        TString timeBranch = TString::Format("time_xstal_%.2d", iRec);
        TString energyBranch = TString::Format("energy_xstal_%.2d", iRec);
        TString posxBranch = TString::Format("posx_xstal_%.2d-posx_xstal_00", iRec);
        TString posyBranch = TString::Format("posy_xstal_%.2d-posy_xstal_00", iRec);
        reader->AddVariable(timeBranch, &times_[iRec]); 
        reader->AddVariable(energyBranch, &E_ratios_[iRec]);
    }
    reader->AddVariable("tan_theta_:=TMath::Abs(pt/pz)", &tan_theta_);
    reader->AddVariable("pz:=TMath::Abs(pz)", &pz_);

    //---prepare the BDTG 
    reader->BookMVA("BDTG method", weightfile);
}

//----------Compute the time from the BDTG------------------------------------------------
float MVATimeComputer::GetMVATime(vector<FTEcalRecHit>* recHits, float E, float pt, float pz, float smearing)
{
    //---set MVA inputs
    pt_ = pt;
    pz_ = TMath::Abs(pz);
    smearing_ = smearing;
    tan_theta_ = pt_/pz_;

    for(unsigned int iRec=0; iRec<9; ++iRec)
    {
        times_[iRec] = iRec < recHits->size() ?
            gRandom->Gaus(recHits->at(iRec).time, smearing) : -1;
        E_ratios_[iRec] = iRec < recHits->size() ?
            recHits->at(iRec).energy/E : -1;
    }
    //---get the regression output
    return (reader->EvaluateRegression("BDTG method"))[0];
}

