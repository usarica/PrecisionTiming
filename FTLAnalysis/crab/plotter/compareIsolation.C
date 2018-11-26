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
  double integral[2]={ hA->Integral(0, nbins+1), hB->Integral(0, nbins+1) };
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

void divideHistByIntegral(TH1F* h){
  if (!h) return;
  float integral = h->Integral(0, h->GetNbinsX()+1);
  if (integral!=0.) h->Scale(1./integral);
}
void appendOverUnderflow(TH1F* h){
  if (!h) return;

  float overval=0, overerr=0;
  for (int bin=h->GetNbinsX(); bin<=h->GetNbinsX()+1; bin++){
    overval+=h->GetBinContent(bin);
    overerr+=pow(h->GetBinError(bin), 2);
    h->SetBinContent(bin, 0);
    h->SetBinError(bin, 0);
  }
  h->SetBinContent(h->GetNbinsX(), overval);
  h->SetBinError(h->GetNbinsX(), sqrt(overerr));

  float underval=0, undererr=0;
  for (int bin=0; bin<=1; bin++){
    underval+=h->GetBinContent(bin);
    undererr+=pow(h->GetBinError(bin), 2);
    h->SetBinContent(bin, 0);
    h->SetBinError(bin, 0);
  }
  h->SetBinContent(1, underval);
  h->SetBinError(1, sqrt(undererr));
}
void postprocessHist(TH1F* h){
  if (!h) return;
  divideHistByIntegral(h);
  appendOverUnderflow(h);
}
void postprocessHist(TH1F& h){ postprocessHist(&h); }

TH1F getPlottableHistogram(TH1F const& h){
  TH1F res(h); res.SetName(Form("%s_plotable", h.GetName()));
  for (int ix=1; ix<=res.GetNbinsX(); ix++){
    double bincontent = res.GetBinContent(ix);
    double binerror = res.GetBinError(ix);
    double binwidth = res.GetXaxis()->GetBinWidth(ix);
    bincontent /= binwidth;
    binerror /= binwidth;
    res.SetBinContent(ix, bincontent);
    res.SetBinError(ix, binerror);
  }
  return res;
}


void compareIsolation(TString strSelection="GenMatched"){
  gStyle->SetOptStat(0);

  TString strSelectionLower = strSelection; strSelectionLower.ToLower();

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

  std::vector<float>* muon_isosumtrackpt_vtx3D_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_nodzcut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_sipcut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_sipcut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_sipcut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_associationrank_2=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_nodzcut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_nodzcut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_nodzcut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_nodzcut_associationrank_2=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_sipcut_unassociated=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_sipcut_associationrank_0=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_sipcut_associationrank_1=nullptr;
  std::vector<float>* muon_isosumtrackpt_vtx4D_sipcut_associationrank_2=nullptr;


  vector<TString> strSampleTitle, strSampleLabel;
  std::vector<TChain*> treeList;

  treeList.push_back(new TChain("muon_tree_30"));
  addFilesInDirectory(treeList.back(), "/hadoop/cms/store/user/usarica/MTD/UPS/MuonIsolation/20181024/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/crab_DY_MuonIsolationMTD_200PU_932_HGCparam_new_v3");
  strSampleTitle.emplace_back("DY");
  strSampleLabel.emplace_back("DY");

  treeList.push_back(new TChain("muon_tree_30"));
  //addFilesInDirectory(treeList.back(), "/hadoop/cms/store/user/usarica/MTD/UPS/MuonIsolation/20181024/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/crab_TTbar_MuonIsolationMTD_200PU_932_HGCparam_new_v1");
  addFilesInDirectory(treeList.back(), "/hadoop/cms/store/user/usarica/MTD/UPS/MuonIsolation/20181024/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/crab_TTbar_ext_MuonIsolationMTD_200PU_932_HGCparam_new_nonlocal_v3");
  strSampleTitle.emplace_back("ttbar");
  strSampleLabel.emplace_back("t#bar{t}");

  TFile* foutput = TFile::Open(Form("muon_%s_iso_plots.root", strSelectionLower.Data()), "recreate");

  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_associationrank_1;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_nodzcut_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_sipcut_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_sipcut_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_sipcut_associationrank_1;

  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1;

  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0;
  std::vector<TH1F> list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1;

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
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_unassociated", &muon_isosumtrackpt_vtx3D_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_associationrank_0", &muon_isosumtrackpt_vtx3D_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_associationrank_1", &muon_isosumtrackpt_vtx3D_associationrank_1);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_nodzcut_unassociated", &muon_isosumtrackpt_vtx3D_nodzcut_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0", &muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1", &muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_sipcut_unassociated", &muon_isosumtrackpt_vtx3D_sipcut_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_sipcut_associationrank_0", &muon_isosumtrackpt_vtx3D_sipcut_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_sipcut_associationrank_1", &muon_isosumtrackpt_vtx3D_sipcut_associationrank_1);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated", &muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0", &muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1", &muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated", &muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0", &muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1", &muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated", &muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0", &muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1", &muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0);
    tree->SetBranchAddress("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1", &muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1);

    std::vector<float> binarray;
    addByLowest(binarray, 0.f, true);
    addByLowest(binarray, 0.02f, true);
    addByLowest(binarray, 100.f, true);
    for (unsigned int ix=1; ix<=100; ix++) addByLowest(binarray, float(0.05)*float(ix), true);
    for (unsigned int ix=1; ix<=50; ix++) addByLowest(binarray, float(0.1)*float(ix)+5.f, true);
    for (unsigned int ix=1; ix<=90; ix++) addByLowest(binarray, float(1.)*float(ix)+10.f, true);
    for (unsigned int ix=1; ix<=180; ix++) addByLowest(binarray, float(5.)*float(ix)+100.f, true);
    const int n_binarray = binarray.size();

    list_muon_isosumtrackpt_vtx3D_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_nodzcut_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_nodzcut_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_sipcut_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_sipcut_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_sipcut_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_sipcut_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_sipcut_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_sipcut_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    TH1F& h_muon_isosumtrackpt_vtx3D_unassociated=list_muon_isosumtrackpt_vtx3D_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_associationrank_0=list_muon_isosumtrackpt_vtx3D_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_associationrank_1=list_muon_isosumtrackpt_vtx3D_associationrank_1.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_nodzcut_unassociated=list_muon_isosumtrackpt_vtx3D_nodzcut_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0=list_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1=list_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_sipcut_unassociated=list_muon_isosumtrackpt_vtx3D_sipcut_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_sipcut_associationrank_0=list_muon_isosumtrackpt_vtx3D_sipcut_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_sipcut_associationrank_1=list_muon_isosumtrackpt_vtx3D_sipcut_associationrank_1.back();
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated=list_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0=list_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1=list_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated=list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0=list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1=list_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated=list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0=list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1=list_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1.back();
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0_Sample%i", is), "", n_binarray-1, binarray.data());
    list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1.emplace_back(Form("muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1_Sample%i", is), "", n_binarray-1, binarray.data());
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0.back();
    TH1F& h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1=list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1.back();

    bool applySignalSelection=(is==0);
    cout << "Will loop over " << nEntries << " entries in Sample " << is << endl;
    for (int ev=0; ev<nEntries; ev++){
      tree->GetEntry(ev);
      for (unsigned int imuon=0; imuon<muon_pt->size(); imuon++){
        if (muon_terr->at(imuon)<=0.) continue; // Only muons with time!

        if (strSelectionLower.Contains("genmatched") &&
            !(
            (applySignalSelection && muonGenMatchedPrompt->at(imuon))
              ||
              (!applySignalSelection && !muonGenMatchedPrompt->at(imuon) && muonGenMatchedJet->at(imuon))
              )
            ) continue;
        if (strSelectionLower.Contains("loose") && !isLooseMuon->at(imuon)) continue;
        if (strSelectionLower.Contains("endcap") && !(fabs(muon_eta->at(imuon))<2.4 && fabs(muon_eta->at(imuon))>=1.5)) continue;
        if (strSelectionLower.Contains("barrel") && !(fabs(muon_eta->at(imuon))<1.5)) continue;
        if (strSelectionLower.Contains("highpt") && !(muon_pt->at(imuon)>20.)) continue;

        if (
          (strSelectionLower.Contains("lowsip") ? (muondIP3DVtx3D->at(imuon)!=0. ? std::abs(muonIP3DVtx3D->at(imuon)/muondIP3DVtx3D->at(imuon))<4. : true) : true)
          //&&
          //()
          ){
          h_muon_isosumtrackpt_vtx3D_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_unassociated.Fill((muon_isosumtrackpt_vtx3D_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_unassociated->at(imuon))/muon_pt->at(imuon));

          h_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_nodzcut_unassociated.Fill((muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_nodzcut_unassociated->at(imuon))/muon_pt->at(imuon));

          h_muon_isosumtrackpt_vtx3D_sipcut_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_sipcut_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_sipcut_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_sipcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_sipcut_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_sipcut_unassociated.Fill((muon_isosumtrackpt_vtx3D_sipcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_sipcut_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_sipcut_unassociated->at(imuon))/muon_pt->at(imuon));

          h_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated->at(imuon))/muon_pt->at(imuon));

          h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated->at(imuon))/muon_pt->at(imuon));

          h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated.Fill((muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated->at(imuon))/muon_pt->at(imuon));

          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated->at(imuon))/muon_pt->at(imuon));

          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated->at(imuon))/muon_pt->at(imuon));

          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1->at(imuon))/muon_pt->at(imuon));
          h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated.Fill((muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1->at(imuon)+muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated->at(imuon))/muon_pt->at(imuon));
        }

      }
    }

    postprocessHist(h_muon_isosumtrackpt_vtx3D_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_unassociated);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_nodzcut_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_nodzcut_unassociated);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_sipcut_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_sipcut_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_sipcut_unassociated);

    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_nodzcut_unassociated);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_dtmuoncut_sipcut_unassociated);

    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_nodzcut_unassociated);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_0);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_associationrank_1);
    postprocessHist(h_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_sipcut_unassociated);

  }

  if (treeList.size()==2){
    std::vector<TGraph*> grlist;

    {
      TString canvasname=Form("cCompare_muoniso_muon_%s_Vtx3D", strSelectionLower.Data());
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

      std::vector<TString> hlabels; hlabels.reserve(6);
      std::vector<TH1F> plotables; plotables.reserve(6);
      plotables.emplace_back(getPlottableHistogram(list_muon_isosumtrackpt_vtx3D_unassociated.at(0)));
      plotables.back().SetLineColor(kBlack); plotables.back().SetLineWidth(2); hlabels.emplace_back("3D A, C and others (#Deltaz<0.1 cm)");
      plotables.emplace_back(getPlottableHistogram(list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated.at(0)));
      plotables.back().SetLineColor(kRed); plotables.back().SetLineWidth(2); hlabels.emplace_back("3D A, C and others (#Deltaz<0.1 cm, |t_{trk}-t_{#mu}|<3#Deltat)");
      plotables.emplace_back(getPlottableHistogram(list_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated.at(0)));
      plotables.back().SetLineColor(kBlue); plotables.back().SetLineWidth(2); hlabels.emplace_back("3D A, C and others (#Deltaz<0.1 cm, no t_{trk} or |t_{trk}-t_{#mu}|<3#Deltat)");
      plotables.emplace_back(getPlottableHistogram(list_muon_isosumtrackpt_vtx3D_unassociated.at(1)));
      plotables.back().SetLineColor(kBlack); plotables.back().SetLineWidth(2); plotables.back().SetLineStyle(7);
      plotables.emplace_back(getPlottableHistogram(list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated.at(1)));
      plotables.back().SetLineColor(kRed); plotables.back().SetLineWidth(2); plotables.back().SetLineStyle(7);
      plotables.emplace_back(getPlottableHistogram(list_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated.at(1)));
      plotables.back().SetLineColor(kBlue); plotables.back().SetLineWidth(2); plotables.back().SetLineStyle(7);
      plotables.front().GetXaxis()->SetRangeUser(0, 0.3);
      plotables.front().GetYaxis()->SetRangeUser(5e-2, 50.);

      const float legend_minX = 0.50;
      const float legend_maxY = 0.90;
      const float legend_maxX = 0.90;
      float legend_minY = legend_maxY - 0.05f*float(plotables.size());
      TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(0.04);
      legend.SetLineColor(0);
      legend.SetLineStyle(0);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(0);

      for (unsigned int ip=0; ip<plotables.size(); ip++){
        if (ip<hlabels.size()) legend.AddEntry(&(plotables.at(ip)), hlabels.at(ip), "l");
        plotables.at(ip).Draw((ip==0 ? "hist" : "histsame"));
      }
      legend.Draw("same");

      canvas.SetLogy();
      canvas.RedrawAxis();
      canvas.Modified();
      canvas.Update();
      canvas.SaveAs(canvasname+".png");
      canvas.SaveAs(canvasname+".pdf");
      foutput->WriteTObject(&canvas);
      canvas.Close();
    }

    TGraph* gr_ROC_vtx3D_unassociated = createROCFromDistributions(
      &(list_muon_isosumtrackpt_vtx3D_unassociated.at(0)),
      &(list_muon_isosumtrackpt_vtx3D_unassociated.at(1)),
      "ROC_vtx3D_unassociated"
    );
    gr_ROC_vtx3D_unassociated->SetMarkerColor(kBlack);
    gr_ROC_vtx3D_unassociated->SetMarkerStyle(28);
    gr_ROC_vtx3D_unassociated->SetMarkerSize(1.2);
    grlist.push_back(gr_ROC_vtx3D_unassociated);
    TGraph* gr_ROC_vtx3D_dtmuoncut_unassociated = createROCFromDistributions(
      &(list_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated.at(0)),
      &(list_muon_isosumtrackpt_vtx3D_dtmuoncut_unassociated.at(1)),
      "ROC_vtx3D_dtmuoncut_unassociated"
    );
    gr_ROC_vtx3D_dtmuoncut_unassociated->SetMarkerColor(kBlue);
    gr_ROC_vtx3D_dtmuoncut_unassociated->SetMarkerStyle(28);
    gr_ROC_vtx3D_dtmuoncut_unassociated->SetMarkerSize(1.2);
    grlist.push_back(gr_ROC_vtx3D_dtmuoncut_unassociated);

    TGraph* gr_ROC_vtx3D_hasdtmuon_dtmuoncut_unassociated = createROCFromDistributions(
      &(list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated.at(0)),
      &(list_muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut_unassociated.at(1)),
      "ROC_vtx3D_hasdtmuon_dtmuoncut_unassociated"
    );
    gr_ROC_vtx3D_hasdtmuon_dtmuoncut_unassociated->SetMarkerColor(kRed);
    gr_ROC_vtx3D_hasdtmuon_dtmuoncut_unassociated->SetMarkerStyle(28);
    gr_ROC_vtx3D_hasdtmuon_dtmuoncut_unassociated->SetMarkerSize(1.2);
    grlist.push_back(gr_ROC_vtx3D_hasdtmuon_dtmuoncut_unassociated);

    for (TGraph*& gr:grlist){
      gr->GetYaxis()->SetTitle(Form("%s eff.", strSampleLabel.at(0).Data()));
      gr->GetXaxis()->SetTitle(Form("%s eff.", strSampleLabel.at(1).Data()));
      gr->SetDrawOption("p");
      foutput->WriteTObject(gr);
    }

    {
      TString canvasname=Form("c_muon_%s_iso_ROC", strSelectionLower.Data());
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

      std::vector<std::pair<TGraph*, TString>> plotables;
      plotables.emplace_back(gr_ROC_vtx3D_unassociated, "3D A, C and others (#Deltaz<0.1 cm)");
      plotables.emplace_back(gr_ROC_vtx3D_hasdtmuon_dtmuoncut_unassociated, "3D A, C and others (#Deltaz<0.1 cm, |t_{trk}-t_{#mu}|<3#Deltat)");
      plotables.emplace_back(gr_ROC_vtx3D_dtmuoncut_unassociated, "3D A, C and others (#Deltaz<0.1 cm, no t_{trk} or |t_{trk}-t_{#mu}|<3#Deltat)");

      const float legend_minY = 0.20;
      const float legend_minX = 0.40;
      const float legend_maxY = legend_minY + float(plotables.size())*0.05;
      const float legend_maxX = 0.90;
      TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(0.04);
      legend.SetLineColor(0);
      legend.SetLineStyle(0);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(0);
      for (unsigned int ip=0; ip<plotables.size(); ip++){
        plotables.at(ip).first->Draw((ip==0 ? "ap" : "psame"));
        legend.AddEntry(plotables.at(ip).first, plotables.at(ip).second, "p");
      }
      legend.Draw("same");

      canvas.RedrawAxis();
      canvas.Modified();
      canvas.Update();
      canvas.SaveAs(canvasname+".png");
      canvas.SaveAs(canvasname+".pdf");
      foutput->WriteTObject(&canvas);
      canvas.Close();
    }


  }

  foutput->Close();
  for (TChain*& tree:treeList) delete tree;
}
