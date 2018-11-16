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
#include "TF1.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

using namespace std;


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

void addPoint(TGraph*& tg, double x, double y){
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

  delete tg;

  tg = new TGraph(xorder.size(), xynew[0], xynew[1]);
  tg->SetName(strname);
  tg->SetTitle(strtitle);
  tg->GetXaxis()->SetTitle(strxtitle);
  tg->GetYaxis()->SetTitle(strytitle);
  for (unsigned int i=0; i<2; i++) delete[] xynew[i];
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

void plotRelTrackVars(){
  gStyle->SetOptStat(0);

  unsigned int const nVars=3;
  unsigned int colors[nVars] ={ kRed, kOrange, kGreen+2/*, kBlue, kViolet*/ };
  vector<TString> strVarTitle{
    "dz_3D",
    "dz_4D",
    "dt_4D"
  };
  vector<TString> strVarLabel{
    "#Deltaz wrt. 3D vtx. (cm)",
    "#Deltaz wrt. 4D vtx. (cm)",
    "#Deltat wrt. 4D vtx. (ns)"
  };

  TChain* tree_sig = new TChain("muon_tree_30");
  TChain* tree_bkg = new TChain("muon_tree_30");
  addFilesInDirectory(tree_sig, "/hadoop/cms/store/user/usarica/MTD/UPS/MuonIsolation/20181024/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/crab_DY_MuonIsolationMTD_200PU_932_HGCparam_nTracks_dz0p1_pt0p9_v7");
  addFilesInDirectory(tree_bkg, "/hadoop/cms/store/user/usarica/MTD/UPS/MuonIsolation/20181024/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/crab_TTbar_MuonIsolationMTD_200PU_932_HGCparam_nTracks_dz0p1_pt0p9_v7");

  vector<TString> strSampleTitle{ "DY","ttbar" };
  vector<TString> strSampleLabel{ "DY","t#bar{t}" };

  std::vector<TH1F> sig_dists; sig_dists.reserve(nVars);
  std::vector<TH1F> bkg_dists; bkg_dists.reserve(nVars);
  std::vector<TGraph*> rocs;
  for (unsigned int i=0; i<nVars; i++){
    if (i<2){
      sig_dists.emplace_back(Form("h_sig_%i", i), strSampleLabel.at(0), 400, -0.4, 0.4);
      bkg_dists.emplace_back(Form("h_bkg_%i", i), strSampleLabel.at(1), 400, -0.4, 0.4);
    }
    else{
      sig_dists.emplace_back(Form("h_sig_%i", i), strSampleLabel.at(0), 300, -0.3, 0.3);
      bkg_dists.emplace_back(Form("h_bkg_%i", i), strSampleLabel.at(1), 300, -0.3, 0.3);
    }
    sig_dists.back().GetXaxis()->SetTitle(strVarLabel.at(i));
    bkg_dists.back().GetXaxis()->SetTitle(strVarLabel.at(i));
  }

  std::vector<float>* tracksZTCut_t=0;
  std::vector<float>* tracksZTCut_dz=0;
  std::vector<float>* tracksZCut_dz=0;
  float vz, vtx4DZ, vtx3DZ, eta, pt, vtx4DT, vtx4DTerr, t0, t0err;
  int eventid, vtx3DIdx, vtx4DIdx;
  bool genMatchedPrompt, genMatchedJet, isLooseMuon;
  tree_sig->SetBranchAddress("tracksZTCut_t", &tracksZTCut_t);
  tree_sig->SetBranchAddress("tracksZTCut_dz", &tracksZTCut_dz);
  tree_sig->SetBranchAddress("tracksZCut_dz", &tracksZCut_dz);
  tree_sig->SetBranchAddress("event", &eventid);
  tree_sig->SetBranchAddress("vtx3DIdx", &vtx3DIdx);
  tree_sig->SetBranchAddress("vtx4DIdx", &vtx4DIdx);
  tree_sig->SetBranchAddress("vz", &vz);
  tree_sig->SetBranchAddress("t0", &t0);
  tree_sig->SetBranchAddress("t0err", &t0err);
  tree_sig->SetBranchAddress("vtx4DT", &vtx4DT);
  tree_sig->SetBranchAddress("vtx4DTerr", &vtx4DTerr);
  tree_sig->SetBranchAddress("vtx4DZ", &vtx4DZ);
  tree_sig->SetBranchAddress("vtx3DZ", &vtx3DZ);
  tree_sig->SetBranchAddress("eta", &eta);
  tree_sig->SetBranchAddress("pt", &pt);
  tree_sig->SetBranchAddress("genMatchedPrompt", &genMatchedPrompt);
  tree_sig->SetBranchAddress("genMatchedJet", &genMatchedJet);
  tree_sig->SetBranchAddress("isLooseMuon", &isLooseMuon);
  //
  tree_bkg->SetBranchAddress("tracksZTCut_t", &tracksZTCut_t);
  tree_bkg->SetBranchAddress("tracksZTCut_dz", &tracksZTCut_dz);
  tree_bkg->SetBranchAddress("tracksZCut_dz", &tracksZCut_dz);
  tree_bkg->SetBranchAddress("event", &eventid);
  tree_bkg->SetBranchAddress("vtx3DIdx", &vtx3DIdx);
  tree_bkg->SetBranchAddress("vtx4DIdx", &vtx4DIdx);
  tree_bkg->SetBranchAddress("vz", &vz);
  tree_bkg->SetBranchAddress("t0", &t0);
  tree_bkg->SetBranchAddress("t0err", &t0err);
  tree_bkg->SetBranchAddress("vtx4DT", &vtx4DT);
  tree_bkg->SetBranchAddress("vtx4DTerr", &vtx4DTerr);
  tree_bkg->SetBranchAddress("vtx4DZ", &vtx4DZ);
  tree_bkg->SetBranchAddress("vtx3DZ", &vtx3DZ);
  tree_bkg->SetBranchAddress("eta", &eta);
  tree_bkg->SetBranchAddress("pt", &pt);
  tree_bkg->SetBranchAddress("genMatchedPrompt", &genMatchedPrompt);
  tree_bkg->SetBranchAddress("genMatchedJet", &genMatchedJet);
  tree_bkg->SetBranchAddress("isLooseMuon", &isLooseMuon);

  TFile* foutput = TFile::Open("trkvtx_associationtest.root", "recreate");

  for (unsigned i = 0; i < nVars; ++i){
    bool use3D = (i==0);
    int lastEvent=-99;
    int lastVtxId=-99;
    for (int ev=0; ev<tree_sig->GetEntries(); ev++){
      tree_sig->GetEntry(ev);
      if (eventid==lastEvent && ((use3D && lastVtxId==(int) vtx3DIdx) || (!use3D && lastVtxId==(int) vtx4DIdx))) continue;
      else{
        lastVtxId = (use3D ? vtx3DIdx : vtx4DIdx);
        lastEvent=eventid;
      }
      float var=-99;
      unsigned int ntrks = (use3D ? tracksZCut_dz->size() : tracksZTCut_dz->size());
      for (unsigned int itrk=0; itrk<ntrks+1; itrk++){
        if (itrk<ntrks){
          switch (i){
          case 0:
            var = tracksZCut_dz->at(itrk);
            break;
          case 1:
            var = tracksZTCut_dz->at(itrk);
            break;
          case 2:
            var = tracksZTCut_t->at(itrk) - vtx4DT;
            break;
          }
        }
        else{
          float vtxZ = (use3D ? vtx3DZ : vtx4DZ);
          switch (i){
          case 0:
          case 1:
            var = vz - vtxZ;
            break;
          case 2:
            var = t0 - vtx4DT;
            break;
          }
        }
        sig_dists[i].Fill(var);
      }
    }
    lastEvent=-99;
    lastVtxId=-99;
    for (int ev=0; ev<tree_bkg->GetEntries(); ev++){
      tree_bkg->GetEntry(ev);
      if (eventid==lastEvent && ((use3D && lastVtxId==(int) vtx3DIdx) || (!use3D && lastVtxId==(int) vtx4DIdx))) continue;
      else{
        lastVtxId = (use3D ? vtx3DIdx : vtx4DIdx);
        lastEvent=eventid;
      }
      float var=-99;
      unsigned int ntrks = (use3D ? tracksZCut_dz->size() : tracksZTCut_dz->size());
      for (unsigned int itrk=0; itrk<ntrks+1; itrk++){
        if (itrk<ntrks){
          switch (i){
          case 0:
            var = tracksZCut_dz->at(itrk);
            break;
          case 1:
            var = tracksZTCut_dz->at(itrk);
            break;
          case 2:
            var = tracksZTCut_t->at(itrk) - vtx4DT;
            break;
          }
        }
        else{
          float vtxZ = (use3D ? vtx3DZ : vtx4DZ);
          switch (i){
          case 0:
          case 1:
            var = vz - vtxZ;
            break;
          case 2:
            var = t0 - vtx4DT;
            break;
          }
        }
        bkg_dists[i].Fill(var);
      }
    }

    cout << sig_dists[i].Integral() << endl;
    cout << bkg_dists[i].Integral() << endl;

    foutput->WriteTObject(&(sig_dists[i]));
    foutput->WriteTObject(&(bkg_dists[i]));
  }

  /*
  TString canvasname="ROCs_" + strSampleTitle.at(0) + "_vs_" + strSampleTitle.at(1);
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

  TLegend legend(0.45, 0.50-0.10/3.*float(rocs.size()), 0.90, 0.50);
  legend.SetBorderSize(0);
  legend.SetTextFont(42);
  legend.SetTextSize(0.03);
  legend.SetLineColor(1);
  legend.SetLineStyle(1);
  legend.SetLineWidth(1);
  legend.SetFillColor(0);
  legend.SetFillStyle(0);

  for (unsigned int i=0; i<nVars; i++){
    if (i==0) rocs.at(i)->Draw("alp");
    else rocs.at(i)->Draw("lpsame");
    legend.AddEntry(rocs.at(i), strVarLabel.at(i), "lp");
  }
  legend.Draw("same");

  canvas.RedrawAxis();
  canvas.Modified();
  canvas.Update();
  canvas.SaveAs(canvasname+".pdf");
  canvas.Close();
  */

  foutput->Close();
  delete tree_sig;
  delete tree_bkg;
}
