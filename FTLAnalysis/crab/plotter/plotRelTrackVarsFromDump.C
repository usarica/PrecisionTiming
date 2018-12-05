#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <ctime>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <functional>
#include <sys/types.h>
#include <dirent.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TObject.h"
#include "TKey.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TSpline.h"
#include "TLine.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TText.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"


// From StdExtensions.h
namespace std{

  template<> struct hash<TString>{
    typedef TString argument_type;
    typedef size_t result_type;
    result_type operator()(argument_type const& arg) const{ return hash<string>{}(arg.Data()); }
  };

}


using namespace std;


struct VariableSpec{
  TString name;
  TString title;
  unsigned int color;

  VariableSpec() :
    name(""),
    title(""),
    color(0)
  {}
  VariableSpec(
    TString const& name_,
    TString const& title_,
    unsigned int const& color_
  ) :
    name(name_),
    title(title_),
    color(color_)
  {}
  VariableSpec(VariableSpec const& other) :
    name(other.name),
    title(other.title),
    color(other.color)
  {}

};


template<typename T> void appendVector(std::vector<T>& a, std::vector<T> const& b){ a.insert(a.end(), b.cbegin(), b.cend()); }

template<typename T> void addByLowest(std::vector<T>& valArray, T val, bool unique){
  bool inserted = false;
  if (unique){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it==val){
        inserted=true;
        break;
      }
    }
  }
  if (!inserted){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it>=val){
        inserted=true;
        valArray.insert(it, val);
        break;
      }
    }
  }
  if (!inserted) valArray.push_back(val);
}
template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, T val, U index){
  bool inserted = false;
  for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it).first>=val){
      inserted=true;
      if ((*it).second!=index) valArray.insert(it, std::pair<T, U>(val, index));
      break;
    }
  }
  if (!inserted) valArray.push_back(std::pair<T, U>(val, index));
}
template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive, bool inputordered){
  if (consecutive){
    bool inserted = false;
    typename std::vector<std::pair<T, U>>::iterator inbegin = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inend = inArray.end();
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first>=(*inbegin).first){
        inserted=true;
        if ((*it).second!=(*inbegin).second) valArray.insert(it, inbegin, inend);
        break;
      }
    }
    if (!inserted) appendVector<std::pair<T, U>>(valArray, inArray);
  }
  else if (!inputordered){
    for (typename std::vector<std::pair<T, U>>::iterator init = inArray.begin(); init<inArray.end(); init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
        if ((*it).first>=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.push_back(*init);
    }
  }
  else if (inArray.size()>0){
    typename std::vector<std::pair<T, U>>::iterator infirst = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inlast = inArray.end()-1;
    typename std::vector<std::pair<T, U>>::iterator valfirst = valArray.begin();
    typename std::vector<std::pair<T, U>>::iterator vallast = valArray.end()-1;
    while ((*valfirst).first<(*infirst).first) valfirst++;
    while ((*vallast).first>=(*inlast).first) vallast--;
    vallast++;
    inlast++;

    for (typename std::vector<std::pair<T, U>>::iterator init = infirst; init<inlast; init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valfirst; it<vallast; it++){
        if ((*it).first>=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.insert(vallast, *init);
    }
  }
}
template<typename T> void addByHighest(std::vector<T>& valArray, T val, bool unique){
  bool inserted = false;
  if (unique){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it==val){
        inserted=true;
        break;
      }
    }
  }
  if (!inserted){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it<=val){
        inserted=true;
        valArray.insert(it, val);
        break;
      }
    }
  }
  if (!inserted) valArray.push_back(val);
}
template<typename T, typename U> void addByHighest(std::vector<std::pair<T, U>>& valArray, T val, U index){
  bool inserted = false;
  for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it).first<=val){
      inserted=true;
      if ((*it).second!=index) valArray.insert(it, std::pair<T, U>(val, index));
      break;
    }
  }
  if (!inserted) valArray.push_back(std::pair<T, U>(val, index));
}
template<typename T, typename U> void addByHighest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive, bool inputordered){
  if (consecutive){
    bool inserted = false;
    typename std::vector<std::pair<T, U>>::iterator inbegin = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inend = inArray.end();
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first<=(*inbegin).first){
        inserted=true;
        if ((*it).second!=(*inbegin).second) valArray.insert(it, inbegin, inend);
        break;
      }
    }
    if (!inserted) appendVector<std::pair<T, U>>(valArray, inArray);
  }
  else if (!inputordered){
    for (typename std::vector<std::pair<T, U>>::iterator init = inArray.begin(); init<inArray.end(); init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
        if ((*it).first<=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.push_back(*init);
    }
  }
  else if (inArray.size()>0){
    typename std::vector<std::pair<T, U>>::iterator infirst = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inlast = inArray.end()-1;
    typename std::vector<std::pair<T, U>>::iterator valfirst = valArray.begin();
    typename std::vector<std::pair<T, U>>::iterator vallast = valArray.end()-1;
    while ((*valfirst).first>(*infirst).first) valfirst++;
    while ((*vallast).first<=(*inlast).first) vallast--;
    vallast++;
    inlast++;

    for (typename std::vector<std::pair<T, U>>::iterator init = infirst; init<inlast; init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valfirst; it<vallast; it++){
        if ((*it).first<=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.insert(vallast, *init);
    }
  }
}

template<typename T> TGraph* makeGraphFromPair(std::vector<std::pair<T, T>> points, TString name){
  if (points.empty()) return nullptr;
  unsigned int nbins = points.size();
  double* xy[2];
  for (unsigned int ix=0; ix<2; ix++) xy[ix] = new double[nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    xy[0][bin] = points[bin].first;
    xy[1][bin] = points[bin].second;
  }
  TGraph* tg = new TGraph(nbins, xy[0], xy[1]);
  tg->SetName(name);
  tg->SetTitle("");
  for (unsigned int ix=0; ix<2; ix++) delete[] xy[ix];
  return tg;
}

void addPoint(TGraph* tg, double x, double y){
  TString strname = tg->GetName();
  TString strtitle = tg->GetTitle();
  TString strxtitle = tg->GetXaxis()->GetTitle();
  TString strytitle = tg->GetYaxis()->GetTitle();

  vector<double> xarray;
  vector<double> yarray;
  xarray.push_back(x);
  yarray.push_back(y);
  for (int ip=0; ip<tg->GetN(); ip++){
    if (tg->GetX()[ip]!=x){
      xarray.push_back(tg->GetX()[ip]);
      yarray.push_back(tg->GetY()[ip]);
    }
  }
  vector<pair<double, int>> xorder;
  for (unsigned int ip=0; ip<xarray.size(); ip++) addByLowest<double, int>(xorder, xarray.at(ip), ip);

  double* xynew[2];
  for (unsigned int i=0; i<2; i++) xynew[i] = new double[xorder.size()];
  for (unsigned int ip=0; ip<xarray.size(); ip++){
    unsigned int pos = xorder[ip].second;
    xynew[0][ip] = xarray[pos];
    xynew[1][ip] = yarray[pos];
  }

  *tg = TGraph(xorder.size(), xynew[0], xynew[1]);
  tg->SetName(strname);
  tg->SetTitle(strtitle);
  tg->GetXaxis()->SetTitle(strxtitle);
  tg->GetYaxis()->SetTitle(strytitle);
  for (unsigned int i=0; i<2; i++) delete[] xynew[i];
}
void addPoint(TGraphErrors* tg, double x, double ex, double y, double ey){
  TString strname = tg->GetName();
  TString strtitle = tg->GetTitle();
  TString strxtitle = tg->GetXaxis()->GetTitle();
  TString strytitle = tg->GetYaxis()->GetTitle();

  const unsigned int nbins = tg->GetN();
  double* xexyey[4]={
    tg->GetX(),
    tg->GetEX(),
    tg->GetY(),
    tg->GetEY()
  };

  double* xexyey_new[4];
  for (unsigned int ix=0; ix<4; ix++) xexyey_new[ix] = new double[nbins+1];

  unsigned int lowbin=0;
  for (unsigned int iy=0; iy<nbins; iy++){ if (xexyey[0][iy]>=x){ lowbin=iy; break; } }

  unsigned int ctr=0;
  if (nbins==0){
    xexyey_new[0][ctr] = x;
    xexyey_new[1][ctr] = ex;
    xexyey_new[2][ctr] = y;
    xexyey_new[3][ctr] = ey;
    ctr++;
  }
  else{
    for (unsigned int iy=0; iy<nbins; iy++){
      if (iy==lowbin){
        xexyey_new[0][ctr] = x;
        xexyey_new[1][ctr] = ex;
        xexyey_new[2][ctr] = y;
        xexyey_new[3][ctr] = ey;
        ctr++;
      }
      for (unsigned int ix=0; ix<4; ix++) xexyey_new[ix][ctr] = xexyey[ix][iy];
      ctr++;
    }
  }

  *tg = TGraphErrors(ctr, xexyey_new[0], xexyey_new[2], xexyey_new[1], xexyey_new[3]);
  tg->SetName(strname);
  tg->SetTitle(strtitle);
  tg->GetXaxis()->SetTitle(strxtitle);
  tg->GetYaxis()->SetTitle(strytitle);
  for (unsigned int ix=0; ix<4; ix++) delete[] xexyey_new[ix];
}

std::vector<TString> lsdir(TString const& indir){
  std::vector<TString> res;

  struct dirent* ep;
  DIR* dp = opendir(indir.Data());
  if (dp != NULL){
    while ((ep = readdir(dp))) res.push_back(ep->d_name);
    closedir(dp);
  }
  else cerr << "Couldn't open the directory" << endl;

  return res;
}

void addFilesInDirectory(TChain* tree, TString dirname){
  vector<TString> lscontent = lsdir(dirname);
  for (auto const& s:lscontent){
    vector<TString> lsversion = lsdir(dirname+"/"+s);
    for (auto const& sv:lsversion){
      vector<TString> lsfnum = lsdir(dirname+"/"+s+"/"+sv);
      for (auto const& sn:lsfnum){
        if (sn.Contains(".root")) tree->Add(dirname+"/"+s+"/"+sv+"/"+sn);
      }
    }
  }
}

TGraph* createROCFromDistributions(TH1* hA, TH1* hB, TString name){
  if (!hA || !hB) return nullptr;
  assert(hA->GetNbinsX()==hB->GetNbinsX());
  const int nbins=hA->GetNbinsX();
  double integral[2]={ hA->Integral(0, nbins), hB->Integral(0, nbins) };
  vector<pair<float, float>> sumWgtsPerBin; sumWgtsPerBin.assign(nbins, pair<float, float>(0, 0));
  for (int ix=1; ix<=nbins; ix++){
    for (int jx=1; jx<=ix; jx++){
      sumWgtsPerBin.at(ix-1).second += hA->GetBinContent(jx)/integral[0];
      sumWgtsPerBin.at(ix-1).first += hB->GetBinContent(jx)/integral[1];
    }
  }
  TGraph* tg=makeGraphFromPair(sumWgtsPerBin, name);
  addPoint(tg, 0, 0);
  addPoint(tg, 1, 1);
  tg->GetYaxis()->SetRangeUser(0, 1);
  tg->GetYaxis()->SetTitle(TString(hA->GetTitle())+" eff.");
  tg->GetXaxis()->SetRangeUser(0, 1);
  tg->GetXaxis()->SetTitle(TString(hB->GetTitle())+" eff.");
  return tg;
}

float calculateSIP3D(float const& ip, float const& d_ip){ return (d_ip==0.f ? 0. : ip/d_ip); }

float getTrackMarkerSize(float const& refval, float const& refmin, float const& refmax){
  constexpr float vmin = 0.7;
  constexpr float vmax = 1.5;
  const float vrefmin = std::max(refmin, 0.f);
  const float vrefmax = std::min(refmax, 6.f);
  float res = std::max(vmin, std::min(vmax, ((refval-vrefmin)/(vrefmax-vrefmin)*(vmax-vmin) + vmin)));
  if (refval>vrefmax) res = vmax + 0.5;
  return res;
}

void plotRelTrackVarsFromDump(){
  gStyle->SetOptStat(0);

  std::unordered_map<TString, VariableSpec> variableList;
  variableList["dz_3D"]=VariableSpec("dz_3D", "#Deltaz wrt. 3D vtx. (cm)", kRed);
  variableList["dz_4D"]=VariableSpec("dz_4D", "#Deltat wrt. 4D vtx. (cm)", kOrange);
  variableList["dt_4D"]=VariableSpec("dt_4D", "#Deltat wrt. 4D vtx. (cm)", kViolet);

  int eventId;
  float simPVZ, simPVT;
  std::vector<float>* muon_pt=nullptr;
  std::vector<float>* muon_eta=nullptr;
  std::vector<float>* muon_phi=nullptr;
  std::vector<float>* muon_px=nullptr;
  std::vector<float>* muon_py=nullptr;
  std::vector<float>* muon_pz=nullptr;
  std::vector<float>* muon_vx=nullptr;
  std::vector<float>* muon_vy=nullptr;
  std::vector<float>* muon_vz=nullptr;
  std::vector<float>* muon_t=nullptr;
  std::vector<float>* muon_terr=nullptr;
  std::vector<unsigned int>* isLooseMuon=nullptr;
  std::vector<unsigned int>* isMediumMuon=nullptr;
  std::vector<unsigned int>* isTightMuon=nullptr;
  std::vector<unsigned int>* muonGenMatched=nullptr;
  std::vector<unsigned int>* muonGenMatchedPrompt=nullptr;
  std::vector<unsigned int>* muonGenMatchedJet=nullptr;
  std::vector<float>* muonGenPt=nullptr;
  std::vector<float>* muonGenEta=nullptr;
  std::vector<float>* muonGenPhi=nullptr;
  std::vector<float>* muonGenJetE=nullptr;
  std::vector<float>* muonGenJetPt=nullptr;
  std::vector<float>* muonGenJetEta=nullptr;
  std::vector<float>* muonGenJetPhi=nullptr;
  std::vector<int>* muonTrkId=nullptr;
  std::vector<unsigned int>* muonVtx3DId=nullptr;
  std::vector<float>* muonIP3DVtx3D=nullptr;
  std::vector<float>* muondIP3DVtx3D=nullptr;
  std::vector<unsigned int>* muonVtx4DId=nullptr;
  std::vector<float>* muonIP3DVtx4D=nullptr;
  std::vector<float>* muondIP3DVtx4D=nullptr;
  std::vector<float>* track_pt=nullptr;
  std::vector<float>* track_eta=nullptr;
  std::vector<float>* track_phi=nullptr;
  std::vector<float>* track_px=nullptr;
  std::vector<float>* track_py=nullptr;
  std::vector<float>* track_pz=nullptr;
  std::vector<float>* track_vx=nullptr;
  std::vector<float>* track_vy=nullptr;
  std::vector<float>* track_vz=nullptr;
  std::vector<float>* track_t=nullptr;
  std::vector<float>* track_terr=nullptr;
  std::vector<unsigned int>* trackVtx3DId=nullptr;
  std::vector<float>* trackIP3DVtx3D=nullptr;
  std::vector<float>* trackdIP3DVtx3D=nullptr;
  std::vector<unsigned int>* trackVtx3DAssociationRank=nullptr;
  std::vector<unsigned int>* trackVtx4DId=nullptr;
  std::vector<float>* trackIP3DVtx4D=nullptr;
  std::vector<float>* trackdIP3DVtx4D=nullptr;
  std::vector<unsigned int>* trackVtx4DAssociationRank=nullptr;
  std::vector<float>* vtx3D_vx=nullptr;
  std::vector<float>* vtx3D_vy=nullptr;
  std::vector<float>* vtx3D_vz=nullptr;
  std::vector<unsigned int>* vtx3D_ntrks=nullptr;
  std::vector<float>* vtx3D_ndof=nullptr;
  std::vector<float>* vtx3D_chisq=nullptr;
  std::vector<float>* vtx4D_vx=nullptr;
  std::vector<float>* vtx4D_vy=nullptr;
  std::vector<float>* vtx4D_vz=nullptr;
  std::vector<float>* vtx4D_t=nullptr;
  std::vector<float>* vtx4D_terr=nullptr;
  std::vector<unsigned int>* vtx4D_ntrks=nullptr;
  std::vector<float>* vtx4D_ndof=nullptr;
  std::vector<float>* vtx4D_chisq=nullptr;

  vector<TString> strSampleTitle, strSampleLabel;
  std::vector<TChain*> treeList;

  treeList.push_back(new TChain("muon_tree_30"));
  //addFilesInDirectory(&(treeList.back()), "/hadoop/cms/store/user/usarica/MTD/UPS/MuonIsolation/20181024/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/crab_DY_MuonIsolationMTD_200PU_932_HGCparam_nTracks_dz0p1_pt0p9_v7");
  treeList.back()->Add("/afs/cern.ch/work/u/usarica/scratch-0/CMSSW_9_3_2/src/PrecisionTiming/FTLAnalysis/test/testfile.root");
  strSampleTitle.emplace_back("DY");
  strSampleLabel.emplace_back("DY");

  /*
  treeList.emplace_back("muon_tree_30");
  addFilesInDirectory(tree_bkg, "/hadoop/cms/store/user/usarica/MTD/UPS/MuonIsolation/20181024/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/crab_TTbar_MuonIsolationMTD_200PU_932_HGCparam_nTracks_dz0p1_pt0p9_v7");
  strSampleTitle.emplace_back("ttbar");
  strSampleLabel.emplace_back("t#bar{t}");
  */

  TFile* foutput = TFile::Open("trkvtx_2Dplots.root", "recreate");

  unsigned int const nSamples = treeList.size();
  for (unsigned int is=0; is<nSamples; is++){
    TChain* tree = treeList.at(is);

    int const nEntries = tree->GetEntries();
    tree->SetBranchAddress("event", &eventId);
    tree->SetBranchAddress("simPVZ", &simPVZ);
    tree->SetBranchAddress("simPVT", &simPVT);
    tree->SetBranchAddress("muon_pt", &muon_pt);
    tree->SetBranchAddress("muon_eta", &muon_eta);
    tree->SetBranchAddress("muon_phi", &muon_phi);
    tree->SetBranchAddress("muon_px", &muon_px);
    tree->SetBranchAddress("muon_py", &muon_py);
    tree->SetBranchAddress("muon_pz", &muon_pz);
    tree->SetBranchAddress("muon_vx", &muon_vx);
    tree->SetBranchAddress("muon_vy", &muon_vy);
    tree->SetBranchAddress("muon_vz", &muon_vz);
    tree->SetBranchAddress("muon_t", &muon_t);
    tree->SetBranchAddress("muon_terr", &muon_terr);
    tree->SetBranchAddress("isLooseMuon", &isLooseMuon);
    tree->SetBranchAddress("isMediumMuon", &isMediumMuon);
    tree->SetBranchAddress("isTightMuon", &isTightMuon);
    tree->SetBranchAddress("muonGenMatched", &muonGenMatched);
    tree->SetBranchAddress("muonGenMatchedPrompt", &muonGenMatchedPrompt);
    tree->SetBranchAddress("muonGenMatchedJet", &muonGenMatchedJet);
    tree->SetBranchAddress("muonGenPt", &muonGenPt);
    tree->SetBranchAddress("muonGenEta", &muonGenEta);
    tree->SetBranchAddress("muonGenPhi", &muonGenPhi);
    tree->SetBranchAddress("muonGenJetE", &muonGenJetE);
    tree->SetBranchAddress("muonGenJetPt", &muonGenJetPt);
    tree->SetBranchAddress("muonGenJetEta", &muonGenJetEta);
    tree->SetBranchAddress("muonGenJetPhi", &muonGenJetPhi);
    tree->SetBranchAddress("muonTrkId", &muonTrkId);
    tree->SetBranchAddress("muonVtx3DId", &muonVtx3DId);
    tree->SetBranchAddress("muonIP3DVtx3D", &muonIP3DVtx3D);
    tree->SetBranchAddress("muondIP3DVtx3D", &muondIP3DVtx3D);
    tree->SetBranchAddress("muonVtx4DId", &muonVtx4DId);
    tree->SetBranchAddress("muonIP3DVtx4D", &muonIP3DVtx4D);
    tree->SetBranchAddress("muondIP3DVtx4D", &muondIP3DVtx4D);
    tree->SetBranchAddress("track_pt", &track_pt);
    tree->SetBranchAddress("track_eta", &track_eta);
    tree->SetBranchAddress("track_phi", &track_phi);
    tree->SetBranchAddress("track_px", &track_px);
    tree->SetBranchAddress("track_py", &track_py);
    tree->SetBranchAddress("track_pz", &track_pz);
    tree->SetBranchAddress("track_vx", &track_vx);
    tree->SetBranchAddress("track_vy", &track_vy);
    tree->SetBranchAddress("track_vz", &track_vz);
    tree->SetBranchAddress("track_t", &track_t);
    tree->SetBranchAddress("track_terr", &track_terr);
    tree->SetBranchAddress("trackVtx3DId", &trackVtx3DId);
    tree->SetBranchAddress("trackIP3DVtx3D", &trackIP3DVtx3D);
    tree->SetBranchAddress("trackdIP3DVtx3D", &trackdIP3DVtx3D);
    tree->SetBranchAddress("trackVtx3DAssociationRank", &trackVtx3DAssociationRank);
    tree->SetBranchAddress("trackVtx4DId", &trackVtx4DId);
    tree->SetBranchAddress("trackIP3DVtx4D", &trackIP3DVtx4D);
    tree->SetBranchAddress("trackdIP3DVtx4D", &trackdIP3DVtx4D);
    tree->SetBranchAddress("trackVtx4DAssociationRank", &trackVtx4DAssociationRank);
    tree->SetBranchAddress("vtx3D_vx", &vtx3D_vx);
    tree->SetBranchAddress("vtx3D_vy", &vtx3D_vy);
    tree->SetBranchAddress("vtx3D_vz", &vtx3D_vz);
    tree->SetBranchAddress("vtx3D_ntrks", &vtx3D_ntrks);
    tree->SetBranchAddress("vtx3D_ndof", &vtx3D_ndof);
    tree->SetBranchAddress("vtx3D_chisq", &vtx3D_chisq);
    tree->SetBranchAddress("vtx4D_vx", &vtx4D_vx);
    tree->SetBranchAddress("vtx4D_vy", &vtx4D_vy);
    tree->SetBranchAddress("vtx4D_vz", &vtx4D_vz);
    tree->SetBranchAddress("vtx4D_t", &vtx4D_t);
    tree->SetBranchAddress("vtx4D_terr", &vtx4D_terr);
    tree->SetBranchAddress("vtx4D_ntrks", &vtx4D_ntrks);
    tree->SetBranchAddress("vtx4D_ndof", &vtx4D_ndof);
    tree->SetBranchAddress("vtx4D_chisq", &vtx4D_chisq);

    // Process entry with most 3D vertices
    int ev_maxVtx3D = -1; unsigned int nMaxVtx3D=0;
    for (int ev=0; ev<nEntries; ev++){
      tree->GetEntry(ev);
      if (nMaxVtx3D<vtx3D_vz->size()){
        nMaxVtx3D = vtx3D_vz->size();
        ev_maxVtx3D = ev;
      }
    }

    // Process entry with most 3D vertices
    tree->GetEntry(ev_maxVtx3D);
    cout << "Plotting entry " << ev_maxVtx3D << " (event = " << eventId << ")" << endl;

    float minSIP3D=std::numeric_limits<float>::max();
    float maxSIP3D=std::numeric_limits<float>::min();

    // Prepare TGraphs
    TGraphErrors vtx3D_dz3D_vs_zerodt;
    vtx3D_dz3D_vs_zerodt.SetName(Form("vtx3D_dz3D_vs_zerodt_Sample%i", is));
    vtx3D_dz3D_vs_zerodt.SetTitle("");
    for (unsigned int ivtx=0; ivtx<vtx3D_vz->size(); ivtx++) addPoint(&vtx3D_dz3D_vs_zerodt, 0, 0, vtx3D_vz->at(ivtx), 0);
    vtx3D_dz3D_vs_zerodt.SetMarkerColor(kOrange+3);
    vtx3D_dz3D_vs_zerodt.SetLineColor(kOrange+3);
    vtx3D_dz3D_vs_zerodt.SetMarkerSize(1.2);
    vtx3D_dz3D_vs_zerodt.SetMarkerStyle(30);
    foutput->WriteTObject(&vtx3D_dz3D_vs_zerodt);

    TGraphErrors vtx4D_dz4D_vs_dt4D;
    vtx4D_dz4D_vs_dt4D.SetName(Form("vtx4D_dz4D_vs_dt4D_Sample%i", is));
    vtx4D_dz4D_vs_dt4D.SetTitle("");
    for (unsigned int ivtx=0; ivtx<vtx4D_vz->size(); ivtx++) addPoint(&vtx4D_dz4D_vs_dt4D, vtx4D_t->at(ivtx), vtx4D_terr->at(ivtx), vtx4D_vz->at(ivtx), 0);
    vtx4D_dz4D_vs_dt4D.SetMarkerColor(kGreen+2);
    vtx4D_dz4D_vs_dt4D.SetLineColor(kGreen+2);
    vtx4D_dz4D_vs_dt4D.SetMarkerSize(1.2);
    vtx4D_dz4D_vs_dt4D.SetMarkerStyle(29);
    foutput->WriteTObject(&vtx4D_dz4D_vs_dt4D);

    TGraphErrors trk3D_dz3D_vs_zerodt_trkAssoc0;
    trk3D_dz3D_vs_zerodt_trkAssoc0.SetName(Form("trk3D_dz3D_vs_zerodt_trkAssoc0_Sample%i", is));
    trk3D_dz3D_vs_zerodt_trkAssoc0.SetTitle("");
    TGraphErrors trk3D_dz3D_vs_zerodt_trkAssoc1;
    trk3D_dz3D_vs_zerodt_trkAssoc1.SetName(Form("trk3D_dz3D_vs_zerodt_trkAssoc1_Sample%i", is));
    trk3D_dz3D_vs_zerodt_trkAssoc1.SetTitle("");
    TGraphErrors trk4D_dz4D_vs_dt4D_trkAssoc0;
    trk4D_dz4D_vs_dt4D_trkAssoc0.SetName(Form("trk4D_dz4D_vs_dt4D_trkAssoc0_Sample%i", is));
    trk4D_dz4D_vs_dt4D_trkAssoc0.SetTitle("");
    TGraphErrors trk4D_dz4D_vs_dt4D_trkAssoc1;
    trk4D_dz4D_vs_dt4D_trkAssoc1.SetName(Form("trk4D_dz4D_vs_dt4D_trkAssoc1_Sample%i", is));
    trk4D_dz4D_vs_dt4D_trkAssoc1.SetTitle("");

    TGraphErrors trk4D_dz4D_vs_zerodt_trkAssoc0;
    trk4D_dz4D_vs_zerodt_trkAssoc0.SetName(Form("trk4D_dz4D_vs_zerodt_trkAssoc0_Sample%i", is));
    trk4D_dz4D_vs_zerodt_trkAssoc0.SetTitle("");
    TGraphErrors trk4D_dz4D_vs_zerodt_trkAssoc2;
    trk4D_dz4D_vs_zerodt_trkAssoc2.SetName(Form("trk4D_dz4D_vs_zerodt_trkAssoc2_Sample%i", is));
    trk4D_dz4D_vs_zerodt_trkAssoc2.SetTitle("");
    for (unsigned int itrk=0; itrk<track_vz->size(); itrk++){
      if (track_terr->at(itrk)<=0.){
        switch (trackVtx3DAssociationRank->at(itrk)){
        case 0:
          addPoint(&trk3D_dz3D_vs_zerodt_trkAssoc0, 0, 0, track_vz->at(itrk), 0);
          break;
        case 1:
          addPoint(&trk3D_dz3D_vs_zerodt_trkAssoc1, 0, 0, track_vz->at(itrk), 0);
          break;
        }
        switch (trackVtx4DAssociationRank->at(itrk)){
        case 0:
          addPoint(&trk4D_dz4D_vs_zerodt_trkAssoc0, 0, 0, track_vz->at(itrk), 0);
          break;
        case 2:
          addPoint(&trk4D_dz4D_vs_zerodt_trkAssoc2, 0, 0, track_vz->at(itrk), 0);
          break;
        }
      }
      else{
        switch (trackVtx4DAssociationRank->at(itrk)){
        case 0:
          addPoint(&trk4D_dz4D_vs_dt4D_trkAssoc0, track_t->at(itrk), track_terr->at(itrk), track_vz->at(itrk), 0);
          break;
        case 1:
          addPoint(&trk4D_dz4D_vs_dt4D_trkAssoc1, track_t->at(itrk), track_terr->at(itrk), track_vz->at(itrk), 0);
          break;
        }
      }
      float valSIP3D = fabs(calculateSIP3D(trackIP3DVtx4D->at(itrk), trackdIP3DVtx4D->at(itrk)));
      minSIP3D=std::min(minSIP3D, valSIP3D);
      maxSIP3D=std::max(maxSIP3D, valSIP3D);
    }
    foutput->WriteTObject(&trk3D_dz3D_vs_zerodt_trkAssoc0);
    foutput->WriteTObject(&trk3D_dz3D_vs_zerodt_trkAssoc1);
    foutput->WriteTObject(&trk4D_dz4D_vs_dt4D_trkAssoc0);
    foutput->WriteTObject(&trk4D_dz4D_vs_dt4D_trkAssoc1);
    foutput->WriteTObject(&trk4D_dz4D_vs_zerodt_trkAssoc0);
    foutput->WriteTObject(&trk4D_dz4D_vs_zerodt_trkAssoc2);

    std::vector<std::pair<unsigned int, float>> muon_pt_ordering;
    for (unsigned int imuon=0; imuon<muon_pt->size(); imuon++) addByHighest(muon_pt_ordering, imuon, muon_pt->at(imuon));

    // Plot vertex - track associations
    for (unsigned int ivtx=0; ivtx<vtx4D_vz->size(); ivtx++){
      TString canvasname=Form("c_trk_dz_vs_dt_Vtx4D%i_Sample%i", ivtx, is);
      TCanvas canvas(canvasname, "", 8, 30, 800, 800);
      canvas.cd();
      gStyle->SetOptStat(0);
      canvas.SetFillColor(0);
      canvas.SetBorderMode(0);
      canvas.SetBorderSize(2);
      canvas.SetTickx(1);
      canvas.SetTicky(1);
      canvas.SetLeftMargin(0.17);
      canvas.SetRightMargin(0.05);
      canvas.SetTopMargin(0.07);
      canvas.SetBottomMargin(0.13);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);

      std::vector</*TLine*/TMarker> trk_vtx_lines_Assoc0; trk_vtx_lines_Assoc0.reserve(trk4D_dz4D_vs_dt4D_trkAssoc0.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_Assoc1; trk_vtx_lines_Assoc1.reserve(trk4D_dz4D_vs_dt4D_trkAssoc1.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_notime_Assoc0; trk_vtx_lines_notime_Assoc0.reserve(trk4D_dz4D_vs_zerodt_trkAssoc0.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_notime_Assoc2; trk_vtx_lines_notime_Assoc2.reserve(trk4D_dz4D_vs_zerodt_trkAssoc2.GetN());

      TMarker model_trk_vtx_lines_Assoc0(0, 0, 20); model_trk_vtx_lines_Assoc0.SetMarkerSize(1); model_trk_vtx_lines_Assoc0.SetMarkerColor(kBlue);
      TMarker model_trk_vtx_lines_Assoc1(0, 0, 24); model_trk_vtx_lines_Assoc1.SetMarkerSize(1); model_trk_vtx_lines_Assoc1.SetMarkerColor(kCyan+2);
      TMarker model_trk_vtx_lines_notime_Assoc0(0, 0, 21); model_trk_vtx_lines_notime_Assoc0.SetMarkerSize(1); model_trk_vtx_lines_notime_Assoc0.SetMarkerColor(kRed);
      TMarker model_trk_vtx_lines_notime_Assoc2(0, 0, 25); model_trk_vtx_lines_notime_Assoc2.SetMarkerSize(1); model_trk_vtx_lines_notime_Assoc2.SetMarkerColor(kViolet);

      TH2F hdummy("hdummy", Form("Event %i, 4D vertex %i", eventId, ivtx), 100, -0.8, 0.8, 100, -35, 35);
      hdummy.SetXTitle("t (ns)");
      hdummy.SetYTitle("z (cm)");
      hdummy.Draw();

      float minSIP3DPerVtx=std::numeric_limits<float>::max();
      float maxSIP3DPerVtx=std::numeric_limits<float>::min();

      for (unsigned int itrk=0; itrk<track_vz->size(); itrk++){
        if (trackVtx4DId->at(itrk)!=ivtx) continue;

        float valSIP3D = fabs(calculateSIP3D(trackIP3DVtx4D->at(itrk), trackdIP3DVtx4D->at(itrk)));
        minSIP3DPerVtx = std::min(minSIP3DPerVtx, valSIP3D);
        maxSIP3DPerVtx = std::max(maxSIP3DPerVtx, valSIP3D);

        TMarker* theMarker=nullptr;
        if (track_terr->at(itrk)>0.){
          switch (trackVtx4DAssociationRank->at(itrk)){
          case 0:
            //trk_vtx_lines_Assoc0.emplace_back(track_t->at(itrk), track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_Assoc0.emplace_back(track_t->at(itrk), track_vz->at(itrk), 20);
            theMarker=&(trk_vtx_lines_Assoc0.back());
            break;
          case 1:
            //trk_vtx_lines_Assoc1.emplace_back(track_t->at(itrk), track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_Assoc1.emplace_back(track_t->at(itrk), track_vz->at(itrk), 24);
            theMarker=&(trk_vtx_lines_Assoc1.back());
            break;
          }
        }
        else{
          switch (trackVtx4DAssociationRank->at(itrk)){
          case 0:
            //trk_vtx_lines_notime_Assoc0.emplace_back(0, track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_notime_Assoc0.emplace_back(0, track_vz->at(itrk), 21);
            theMarker=&(trk_vtx_lines_notime_Assoc0.back());
            break;
          case 2:
            //trk_vtx_lines_notime_Assoc2.emplace_back(0, track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_notime_Assoc2.emplace_back(0, track_vz->at(itrk), 25);
            theMarker=&(trk_vtx_lines_notime_Assoc2.back());
            break;
          }
        }
        if (theMarker){
          theMarker->SetMarkerSize(getTrackMarkerSize(valSIP3D, minSIP3D, maxSIP3D));
        }
      }
      TMarker vtxMarker(vtx4D_t->at(ivtx), vtx4D_vz->at(ivtx), 33);
      //vtxMarker.SetMarkerStyle(20);
      vtxMarker.SetMarkerSize(3);
      vtxMarker.SetMarkerColor(kBlack);
      vtxMarker.Draw("same");
      TMarker pvMarker(simPVT, simPVZ, 29);
      //pvMarker.SetMarkerStyle(29);
      pvMarker.SetMarkerSize(3);
      pvMarker.SetMarkerColor(kGreen+2);
      pvMarker.Draw("same");
      //vtx4D_dz4D_vs_dt4D.Draw("e1psame");
      for (auto& line:trk_vtx_lines_notime_Assoc2){
        line.SetMarkerColor(kViolet);
        //line.SetLineStyle(7);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_notime_Assoc0){
        line.SetMarkerColor(kRed);
        //line.SetLineStyle(2);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_Assoc1){
        line.SetMarkerColor(kCyan+2);
        //line.SetLineStyle(7);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_Assoc0){
        line.SetMarkerColor(kBlue);
        //line.SetLineStyle(2);
        line.Draw("same");
      }

      const float legend_minX = 0.50;
      const float legend_maxY = 0.90;
      const float legend_maxX = 0.90;
      float legend_minY = legend_maxY;
      if (!trk_vtx_lines_Assoc0.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_Assoc1.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_notime_Assoc0.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_notime_Assoc2.empty()) legend_minY -= 0.05;
      TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(0.04);
      legend.SetLineColor(0);
      legend.SetLineStyle(0);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(0);

      if (!trk_vtx_lines_Assoc0.empty()) legend.AddEntry(&model_trk_vtx_lines_Assoc0, Form("4D tracks, A (n=%lu)", trk_vtx_lines_Assoc0.size()), "p");
      if (!trk_vtx_lines_Assoc1.empty()) legend.AddEntry(&model_trk_vtx_lines_Assoc1, Form("4D tracks, B (n=%lu)", trk_vtx_lines_Assoc1.size()), "p");
      if (!trk_vtx_lines_notime_Assoc0.empty()) legend.AddEntry(&model_trk_vtx_lines_notime_Assoc0, Form("3D tracks, A (n=%lu)", trk_vtx_lines_notime_Assoc0.size()), "p");
      if (!trk_vtx_lines_notime_Assoc2.empty()) legend.AddEntry(&model_trk_vtx_lines_notime_Assoc2, Form("3D tracks, C (n=%lu)", trk_vtx_lines_notime_Assoc2.size()), "p");
      legend.Draw("same");

      TText* text=nullptr;

      const float paveSIP3D_minX = 0.20;
      const float paveSIP3D_maxX = legend_minX-0.05;
      const float paveSIP3D_maxY = legend_maxY;
      const float paveSIP3D_minY = paveSIP3D_maxY - 0.05;
      TPaveText paveSIP3D(paveSIP3D_minX, paveSIP3D_minY, paveSIP3D_maxX, paveSIP3D_maxY, "brNDC");
      paveSIP3D.SetBorderSize(0);
      paveSIP3D.SetFillStyle(0);
      paveSIP3D.SetTextAlign(12);
      paveSIP3D.SetTextFont(42);
      paveSIP3D.SetTextSize(0.04);
      text = paveSIP3D.AddText(0.02, 0.45, Form("SIP_{3D} = [%.1f, %.1f]", minSIP3DPerVtx, maxSIP3DPerVtx));
      paveSIP3D.Draw("same");

      canvas.RedrawAxis();
      canvas.Modified();
      canvas.Update();
      canvas.SaveAs(canvasname+".png");
      foutput->WriteTObject(&canvas);
      canvas.Close();
    }

    // Plot muon - track dR
    unsigned int jmuon=0;
    for(std::pair<unsigned int, float> const& iorder:muon_pt_ordering){
      unsigned int const& imuon = iorder.first;
      unsigned int const& ivtx = muonVtx4DId->at(imuon);

      TString canvasname=Form("c_muon_dR_dt_Vtx4D%i_muon%i_Sample%i", ivtx, jmuon, is);
      TCanvas canvas(canvasname, "", 8, 30, 800, 800);
      canvas.cd();
      gStyle->SetOptStat(0);
      canvas.SetFillColor(0);
      canvas.SetBorderMode(0);
      canvas.SetBorderSize(2);
      canvas.SetTickx(1);
      canvas.SetTicky(1);
      canvas.SetLeftMargin(0.17);
      canvas.SetRightMargin(0.05);
      canvas.SetTopMargin(0.07);
      canvas.SetBottomMargin(0.13);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);

      std::vector</*TLine*/TMarker> trk_vtx_lines_Assoc0; trk_vtx_lines_Assoc0.reserve(trk4D_dz4D_vs_dt4D_trkAssoc0.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_Assoc1; trk_vtx_lines_Assoc1.reserve(trk4D_dz4D_vs_dt4D_trkAssoc1.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_notime_Assoc0; trk_vtx_lines_notime_Assoc0.reserve(trk4D_dz4D_vs_zerodt_trkAssoc0.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_notime_Assoc2; trk_vtx_lines_notime_Assoc2.reserve(trk4D_dz4D_vs_zerodt_trkAssoc2.GetN());

      TMarker model_trk_vtx_lines_Assoc0(0, 0, 20); model_trk_vtx_lines_Assoc0.SetMarkerSize(1); model_trk_vtx_lines_Assoc0.SetMarkerColor(kBlue);
      TMarker model_trk_vtx_lines_Assoc1(0, 0, 24); model_trk_vtx_lines_Assoc1.SetMarkerSize(1); model_trk_vtx_lines_Assoc1.SetMarkerColor(kCyan+2);
      TMarker model_trk_vtx_lines_notime_Assoc0(0, 0, 21); model_trk_vtx_lines_notime_Assoc0.SetMarkerSize(1); model_trk_vtx_lines_notime_Assoc0.SetMarkerColor(kRed);
      TMarker model_trk_vtx_lines_notime_Assoc2(0, 0, 25); model_trk_vtx_lines_notime_Assoc2.SetMarkerSize(1); model_trk_vtx_lines_notime_Assoc2.SetMarkerColor(kViolet);

      TVector3 pMuon(muon_px->at(imuon), muon_py->at(imuon), muon_pz->at(imuon));

      float minSIP3DPerVtx=std::numeric_limits<float>::max();
      float maxSIP3DPerVtx=std::numeric_limits<float>::min();

      float maxAbsDeltaT=0;
      float maxDeltaR=0;

      for (unsigned int itrk=0; itrk<track_vz->size(); itrk++){
        if (trackVtx4DId->at(itrk)!=ivtx) continue;
        if (itrk==muonTrkId->at(imuon)) continue;
        TVector3 pTrk(track_px->at(itrk), track_py->at(itrk), track_pz->at(itrk));
        float dR = pMuon.DeltaR(pTrk);
        maxDeltaR = std::max(maxDeltaR, dR);

        float valSIP3D = fabs(calculateSIP3D(trackIP3DVtx4D->at(itrk), trackdIP3DVtx4D->at(itrk)));
        minSIP3DPerVtx = std::min(minSIP3DPerVtx, valSIP3D);
        maxSIP3DPerVtx = std::max(maxSIP3DPerVtx, valSIP3D);

        TMarker* theMarker=nullptr;
        if (track_terr->at(itrk)>0.){
          switch (trackVtx4DAssociationRank->at(itrk)){
          case 0:
            //trk_vtx_lines_Assoc0.emplace_back(track_t->at(itrk), track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_Assoc0.emplace_back(track_t->at(itrk) - muon_t->at(imuon), dR, 20);
            maxAbsDeltaT = std::max(maxAbsDeltaT, fabs(track_t->at(itrk) - muon_t->at(imuon)));
            theMarker=&(trk_vtx_lines_Assoc0.back());
            break;
          case 1:
            //trk_vtx_lines_Assoc1.emplace_back(track_t->at(itrk), track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_Assoc1.emplace_back(track_t->at(itrk) - muon_t->at(imuon), dR, 24);
            maxAbsDeltaT = std::max(maxAbsDeltaT, fabs(track_t->at(itrk) - muon_t->at(imuon)));
            theMarker=&(trk_vtx_lines_Assoc1.back());
            break;
          }
        }
        else{
          switch (trackVtx4DAssociationRank->at(itrk)){
          case 0:
            //trk_vtx_lines_notime_Assoc0.emplace_back(0, track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_notime_Assoc0.emplace_back(0, dR, 21);
            theMarker=&(trk_vtx_lines_notime_Assoc0.back());
            break;
          case 2:
            //trk_vtx_lines_notime_Assoc2.emplace_back(0, track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_notime_Assoc2.emplace_back(0, dR, 25);
            theMarker=&(trk_vtx_lines_notime_Assoc2.back());
            break;
          }
        }
        if (theMarker){
          theMarker->SetMarkerSize(getTrackMarkerSize(valSIP3D, minSIP3D, maxSIP3D));
        }
      }

      TH2F hdummy("hdummy", Form("Event %i, muon %i, 4D vertex %i", eventId, jmuon, ivtx), 100, -maxAbsDeltaT*1.1, maxAbsDeltaT*1.1, 100, 0, maxDeltaR*1.5);
      hdummy.SetXTitle("t_{trk}-t_{#mu} (ns)");
      hdummy.SetYTitle("#DeltaR");
      hdummy.Draw();

      for (auto& line:trk_vtx_lines_notime_Assoc2){
        line.SetMarkerColor(kViolet);
        //line.SetLineStyle(7);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_notime_Assoc0){
        line.SetMarkerColor(kRed);
        //line.SetLineStyle(2);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_Assoc1){
        line.SetMarkerColor(kCyan+2);
        //line.SetLineStyle(7);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_Assoc0){
        line.SetMarkerColor(kBlue);
        //line.SetLineStyle(2);
        line.Draw("same");
      }

      const float legend_minX = 0.50;
      const float legend_maxY = 0.90;
      const float legend_maxX = 0.90;
      float legend_minY = legend_maxY;
      if (!trk_vtx_lines_Assoc0.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_Assoc1.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_notime_Assoc0.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_notime_Assoc2.empty()) legend_minY -= 0.05;
      TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(0.04);
      legend.SetLineColor(0);
      legend.SetLineStyle(0);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(0);

      if (!trk_vtx_lines_Assoc0.empty()) legend.AddEntry(&model_trk_vtx_lines_Assoc0, Form("4D tracks, A (n=%lu)", trk_vtx_lines_Assoc0.size()), "p");
      if (!trk_vtx_lines_Assoc1.empty()) legend.AddEntry(&model_trk_vtx_lines_Assoc1, Form("4D tracks, B (n=%lu)", trk_vtx_lines_Assoc1.size()), "p");
      if (!trk_vtx_lines_notime_Assoc0.empty()) legend.AddEntry(&model_trk_vtx_lines_notime_Assoc0, Form("3D tracks, A (n=%lu)", trk_vtx_lines_notime_Assoc0.size()), "p");
      if (!trk_vtx_lines_notime_Assoc2.empty()) legend.AddEntry(&model_trk_vtx_lines_notime_Assoc2, Form("3D tracks, C (n=%lu)", trk_vtx_lines_notime_Assoc2.size()), "p");
      legend.Draw("same");

      TText* text=nullptr;

      const float paveSIP3D_minX = 0.20;
      const float paveSIP3D_maxX = legend_minX-0.05;
      const float paveSIP3D_maxY = legend_maxY;
      const float paveSIP3D_minY = paveSIP3D_maxY - 0.05;
      TPaveText paveSIP3D(paveSIP3D_minX, paveSIP3D_minY, paveSIP3D_maxX, paveSIP3D_maxY, "brNDC");
      paveSIP3D.SetBorderSize(0);
      paveSIP3D.SetFillStyle(0);
      paveSIP3D.SetTextAlign(12);
      paveSIP3D.SetTextFont(42);
      paveSIP3D.SetTextSize(0.04);
      text = paveSIP3D.AddText(0.02, 0.45, Form("SIP_{3D} = [%.1f, %.1f]", minSIP3DPerVtx, maxSIP3DPerVtx));
      paveSIP3D.Draw("same");

      canvas.RedrawAxis();
      canvas.Modified();
      canvas.Update();
      canvas.SaveAs(canvasname+".png");
      foutput->WriteTObject(&canvas);
      canvas.Close();

      jmuon++;
    }

    // Plot track dR vs pt/pt_mu
    jmuon=0;
    for (std::pair<unsigned int, float> const& iorder:muon_pt_ordering){
      unsigned int const& imuon = iorder.first;
      unsigned int const& ivtx = muonVtx4DId->at(imuon);

      TString canvasname=Form("c_dR_ptRatio_Vtx4D%i_muon%i_Sample%i", ivtx, jmuon, is);
      TCanvas canvas(canvasname, "", 8, 30, 800, 800);
      canvas.cd();
      gStyle->SetOptStat(0);
      canvas.SetFillColor(0);
      canvas.SetBorderMode(0);
      canvas.SetBorderSize(2);
      canvas.SetTickx(1);
      canvas.SetTicky(1);
      canvas.SetLeftMargin(0.17);
      canvas.SetRightMargin(0.05);
      canvas.SetTopMargin(0.07);
      canvas.SetBottomMargin(0.13);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);

      std::vector</*TLine*/TMarker> trk_vtx_lines_Assoc0; trk_vtx_lines_Assoc0.reserve(trk4D_dz4D_vs_dt4D_trkAssoc0.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_Assoc1; trk_vtx_lines_Assoc1.reserve(trk4D_dz4D_vs_dt4D_trkAssoc1.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_notime_Assoc0; trk_vtx_lines_notime_Assoc0.reserve(trk4D_dz4D_vs_zerodt_trkAssoc0.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_notime_Assoc2; trk_vtx_lines_notime_Assoc2.reserve(trk4D_dz4D_vs_zerodt_trkAssoc2.GetN());

      TMarker model_trk_vtx_lines_Assoc0(0, 0, 20); model_trk_vtx_lines_Assoc0.SetMarkerSize(1); model_trk_vtx_lines_Assoc0.SetMarkerColor(kBlue);
      TMarker model_trk_vtx_lines_Assoc1(0, 0, 24); model_trk_vtx_lines_Assoc1.SetMarkerSize(1); model_trk_vtx_lines_Assoc1.SetMarkerColor(kCyan+2);
      TMarker model_trk_vtx_lines_notime_Assoc0(0, 0, 21); model_trk_vtx_lines_notime_Assoc0.SetMarkerSize(1); model_trk_vtx_lines_notime_Assoc0.SetMarkerColor(kRed);
      TMarker model_trk_vtx_lines_notime_Assoc2(0, 0, 25); model_trk_vtx_lines_notime_Assoc2.SetMarkerSize(1); model_trk_vtx_lines_notime_Assoc2.SetMarkerColor(kViolet);

      TVector3 pMuon(muon_px->at(imuon), muon_py->at(imuon), muon_pz->at(imuon));

      float minSIP3DPerVtx=std::numeric_limits<float>::max();
      float maxSIP3DPerVtx=std::numeric_limits<float>::min();

      float maxPtRatio=0;
      float maxDeltaR=0;

      for (unsigned int itrk=0; itrk<track_vz->size(); itrk++){
        if (trackVtx4DId->at(itrk)!=ivtx) continue;
        if (itrk==muonTrkId->at(imuon)) continue;
        TVector3 pTrk(track_px->at(itrk), track_py->at(itrk), track_pz->at(itrk));

        float dR = pMuon.DeltaR(pTrk);
        float ptRatio = pTrk.Pt() / pMuon.Pt();
        maxDeltaR = std::max(maxDeltaR, dR);
        maxPtRatio = std::max(maxPtRatio, ptRatio);

        float valSIP3D = fabs(calculateSIP3D(trackIP3DVtx4D->at(itrk), trackdIP3DVtx4D->at(itrk)));
        minSIP3DPerVtx = std::min(minSIP3DPerVtx, valSIP3D);
        maxSIP3DPerVtx = std::max(maxSIP3DPerVtx, valSIP3D);

        TMarker* theMarker=nullptr;
        if (track_terr->at(itrk)>0.){
          switch (trackVtx4DAssociationRank->at(itrk)){
          case 0:
            //trk_vtx_lines_Assoc0.emplace_back(track_t->at(itrk), track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_Assoc0.emplace_back(ptRatio, dR, 20);
            theMarker=&(trk_vtx_lines_Assoc0.back());
            break;
          case 1:
            //trk_vtx_lines_Assoc1.emplace_back(track_t->at(itrk), track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_Assoc1.emplace_back(ptRatio, dR, 24);
            theMarker=&(trk_vtx_lines_Assoc1.back());
            break;
          }
        }
        else{
          switch (trackVtx4DAssociationRank->at(itrk)){
          case 0:
            //trk_vtx_lines_notime_Assoc0.emplace_back(0, track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_notime_Assoc0.emplace_back(ptRatio, dR, 21);
            theMarker=&(trk_vtx_lines_notime_Assoc0.back());
            break;
          case 2:
            //trk_vtx_lines_notime_Assoc2.emplace_back(0, track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
            trk_vtx_lines_notime_Assoc2.emplace_back(ptRatio, dR, 25);
            theMarker=&(trk_vtx_lines_notime_Assoc2.back());
            break;
          }
        }
        if (theMarker){
          theMarker->SetMarkerSize(getTrackMarkerSize(valSIP3D, minSIP3D, maxSIP3D));
        }
      }

      TH2F hdummy("hdummy", Form("Event %i, muon %i, 4D vertex %i", eventId, jmuon, ivtx), 100, 0, maxPtRatio*1.1, 100, 0, maxDeltaR*1.5);
      hdummy.SetXTitle("p_{T}^{trk} / p_{T}^{#mu}");
      hdummy.SetYTitle("#DeltaR");
      hdummy.Draw();

      for (auto& line:trk_vtx_lines_notime_Assoc2){
        line.SetMarkerColor(kViolet);
        //line.SetLineStyle(7);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_notime_Assoc0){
        line.SetMarkerColor(kRed);
        //line.SetLineStyle(2);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_Assoc1){
        line.SetMarkerColor(kCyan+2);
        //line.SetLineStyle(7);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_Assoc0){
        line.SetMarkerColor(kBlue);
        //line.SetLineStyle(2);
        line.Draw("same");
      }

      const float legend_minX = 0.50;
      const float legend_maxY = 0.90;
      const float legend_maxX = 0.90;
      float legend_minY = legend_maxY;
      if (!trk_vtx_lines_Assoc0.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_Assoc1.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_notime_Assoc0.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_notime_Assoc2.empty()) legend_minY -= 0.05;
      TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(0.04);
      legend.SetLineColor(0);
      legend.SetLineStyle(0);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(0);

      if (!trk_vtx_lines_Assoc0.empty()) legend.AddEntry(&model_trk_vtx_lines_Assoc0, Form("4D tracks, A (n=%lu)", trk_vtx_lines_Assoc0.size()), "p");
      if (!trk_vtx_lines_Assoc1.empty()) legend.AddEntry(&model_trk_vtx_lines_Assoc1, Form("4D tracks, B (n=%lu)", trk_vtx_lines_Assoc1.size()), "p");
      if (!trk_vtx_lines_notime_Assoc0.empty()) legend.AddEntry(&model_trk_vtx_lines_notime_Assoc0, Form("3D tracks, A (n=%lu)", trk_vtx_lines_notime_Assoc0.size()), "p");
      if (!trk_vtx_lines_notime_Assoc2.empty()) legend.AddEntry(&model_trk_vtx_lines_notime_Assoc2, Form("3D tracks, C (n=%lu)", trk_vtx_lines_notime_Assoc2.size()), "p");
      legend.Draw("same");

      TText* text=nullptr;

      const float paveSIP3D_minX = 0.20;
      const float paveSIP3D_maxX = legend_minX-0.05;
      const float paveSIP3D_maxY = legend_maxY;
      const float paveSIP3D_minY = paveSIP3D_maxY - 0.05;
      TPaveText paveSIP3D(paveSIP3D_minX, paveSIP3D_minY, paveSIP3D_maxX, paveSIP3D_maxY, "brNDC");
      paveSIP3D.SetBorderSize(0);
      paveSIP3D.SetFillStyle(0);
      paveSIP3D.SetTextAlign(12);
      paveSIP3D.SetTextFont(42);
      paveSIP3D.SetTextSize(0.04);
      text = paveSIP3D.AddText(0.02, 0.45, Form("SIP_{3D} = [%.1f, %.1f]", minSIP3DPerVtx, maxSIP3DPerVtx));
      paveSIP3D.Draw("same");

      canvas.RedrawAxis();
      canvas.Modified();
      canvas.Update();
      canvas.SaveAs(canvasname+".png");
      foutput->WriteTObject(&canvas);
      canvas.Close();

      jmuon++;
    }

    // Plot track dR vs pt/pt_mu (3D vertex)
    jmuon=0;
    for (std::pair<unsigned int, float> const& iorder:muon_pt_ordering){
      unsigned int const& imuon = iorder.first;
      unsigned int const& ivtx = muonVtx3DId->at(imuon);

      TString canvasname=Form("c_dR_ptRatio_Vtx3D%i_muon%i_Sample%i", ivtx, jmuon, is);
      TCanvas canvas(canvasname, "", 8, 30, 800, 800);
      canvas.cd();
      gStyle->SetOptStat(0);
      canvas.SetFillColor(0);
      canvas.SetBorderMode(0);
      canvas.SetBorderSize(2);
      canvas.SetTickx(1);
      canvas.SetTicky(1);
      canvas.SetLeftMargin(0.17);
      canvas.SetRightMargin(0.05);
      canvas.SetTopMargin(0.07);
      canvas.SetBottomMargin(0.13);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);
      canvas.SetFrameFillStyle(0);
      canvas.SetFrameBorderMode(0);

      std::vector</*TLine*/TMarker> trk_vtx_lines_notime_Assoc0; trk_vtx_lines_notime_Assoc0.reserve(trk3D_dz3D_vs_zerodt_trkAssoc0.GetN());
      std::vector</*TLine*/TMarker> trk_vtx_lines_notime_Assoc1; trk_vtx_lines_notime_Assoc1.reserve(trk3D_dz3D_vs_zerodt_trkAssoc1.GetN());

      TMarker model_trk_vtx_lines_notime_Assoc0(0, 0, 21); model_trk_vtx_lines_notime_Assoc0.SetMarkerSize(1); model_trk_vtx_lines_notime_Assoc0.SetMarkerColor(kRed);
      TMarker model_trk_vtx_lines_notime_Assoc1(0, 0, 25); model_trk_vtx_lines_notime_Assoc1.SetMarkerSize(1); model_trk_vtx_lines_notime_Assoc1.SetMarkerColor(kViolet);

      TVector3 pMuon(muon_px->at(imuon), muon_py->at(imuon), muon_pz->at(imuon));

      float minSIP3DPerVtx=std::numeric_limits<float>::max();
      float maxSIP3DPerVtx=std::numeric_limits<float>::min();

      float maxPtRatio=0;
      float maxDeltaR=0;

      for (unsigned int itrk=0; itrk<track_vz->size(); itrk++){
        if (trackVtx3DId->at(itrk)!=ivtx) continue;
        if (itrk==muonTrkId->at(imuon)) continue;
        TVector3 pTrk(track_px->at(itrk), track_py->at(itrk), track_pz->at(itrk));

        float dR = pMuon.DeltaR(pTrk);
        float ptRatio = pTrk.Pt() / pMuon.Pt();
        maxDeltaR = std::max(maxDeltaR, dR);
        maxPtRatio = std::max(maxPtRatio, ptRatio);

        float valSIP3D = fabs(calculateSIP3D(trackIP3DVtx3D->at(itrk), trackdIP3DVtx3D->at(itrk)));
        minSIP3DPerVtx = std::min(minSIP3DPerVtx, valSIP3D);
        maxSIP3DPerVtx = std::max(maxSIP3DPerVtx, valSIP3D);

        TMarker* theMarker=nullptr;
        switch (trackVtx4DAssociationRank->at(itrk)){
        case 0:
          //trk_vtx_lines_notime_Assoc0.emplace_back(0, track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
          trk_vtx_lines_notime_Assoc0.emplace_back(ptRatio, dR, 21);
          theMarker=&(trk_vtx_lines_notime_Assoc0.back());
          break;
        case 1:
          //trk_vtx_lines_notime_Assoc2.emplace_back(0, track_vz->at(itrk), vtx4D_t->at(trackVtx4DId->at(itrk)), vtx4D_vz->at(trackVtx4DId->at(itrk)));
          trk_vtx_lines_notime_Assoc1.emplace_back(ptRatio, dR, 25);
          theMarker=&(trk_vtx_lines_notime_Assoc1.back());
          break;
        }
        if (theMarker){
          theMarker->SetMarkerSize(getTrackMarkerSize(valSIP3D, minSIP3D, maxSIP3D));
        }
      }

      TH2F hdummy("hdummy", Form("Event %i, muon %i, 3D vertex %i", eventId, jmuon, ivtx), 100, 0, maxPtRatio*1.1, 100, 0, maxDeltaR*1.5);
      hdummy.SetXTitle("p_{T}^{trk} / p_{T}^{#mu}");
      hdummy.SetYTitle("#DeltaR");
      hdummy.Draw();

      for (auto& line:trk_vtx_lines_notime_Assoc1){
        line.SetMarkerColor(kViolet);
        //line.SetLineStyle(7);
        line.Draw("same");
      }
      for (auto& line:trk_vtx_lines_notime_Assoc0){
        line.SetMarkerColor(kRed);
        //line.SetLineStyle(2);
        line.Draw("same");
      }

      const float legend_minX = 0.50;
      const float legend_maxY = 0.90;
      const float legend_maxX = 0.90;
      float legend_minY = legend_maxY;
      if (!trk_vtx_lines_notime_Assoc0.empty()) legend_minY -= 0.05;
      if (!trk_vtx_lines_notime_Assoc1.empty()) legend_minY -= 0.05;
      TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(0.04);
      legend.SetLineColor(0);
      legend.SetLineStyle(0);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(0);

      if (!trk_vtx_lines_notime_Assoc0.empty()) legend.AddEntry(&model_trk_vtx_lines_notime_Assoc0, Form("3D tracks, A (n=%lu)", trk_vtx_lines_notime_Assoc0.size()), "p");
      if (!trk_vtx_lines_notime_Assoc1.empty()) legend.AddEntry(&model_trk_vtx_lines_notime_Assoc1, Form("3D tracks, C (n=%lu)", trk_vtx_lines_notime_Assoc1.size()), "p");
      legend.Draw("same");

      TText* text=nullptr;

      const float paveSIP3D_minX = 0.20;
      const float paveSIP3D_maxX = legend_minX-0.05;
      const float paveSIP3D_maxY = legend_maxY;
      const float paveSIP3D_minY = paveSIP3D_maxY - 0.05;
      TPaveText paveSIP3D(paveSIP3D_minX, paveSIP3D_minY, paveSIP3D_maxX, paveSIP3D_maxY, "brNDC");
      paveSIP3D.SetBorderSize(0);
      paveSIP3D.SetFillStyle(0);
      paveSIP3D.SetTextAlign(12);
      paveSIP3D.SetTextFont(42);
      paveSIP3D.SetTextSize(0.04);
      text = paveSIP3D.AddText(0.02, 0.45, Form("SIP_{3D} = [%.1f, %.1f]", minSIP3DPerVtx, maxSIP3DPerVtx));
      paveSIP3D.Draw("same");

      canvas.RedrawAxis();
      canvas.Modified();
      canvas.Update();
      canvas.SaveAs(canvasname+".png");
      foutput->WriteTObject(&canvas);
      canvas.Close();

      jmuon++;
    }

  }

  foutput->Close();
  for (TChain*& tree:treeList) delete tree;
}
