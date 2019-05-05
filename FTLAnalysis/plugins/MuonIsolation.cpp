#ifndef _FTL_MUON_ISOLATION_
#define _FTL_MUON_ISOLATION_

// system include files
#include <cmath>
#include <cstdlib>
#include <memory>
#include <unordered_map>
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "TTree.h"
#include "TRandom.h"

#include "PrecisionTiming/FTLAnalysis/interface/FTLMuonIsoTree.h"


template<typename T> bool checkListVariable(std::vector<T> const& list, T const& var){
  // Look for exact match
  for (T const& v:list){ if (v==var) return true; }
  return false;
}

struct VertexInformation;
struct TrackInformation;
struct MuonInformation;

struct VertexInformation{
  reco::Vertex const* ptr;
  bool use3D;
  bool isUsable;
  unsigned int ndof;
  float chisq;
  float vx;
  float vy;
  float vz;
  float t;
  float terr;

  std::vector<TrackInformation const*> associatedTrackInfos;
  std::vector<MuonInformation const*> associatedMuonInfos;

  bool hasTime() const{ return (terr>0.); }
  void setup(){
    isUsable = (ptr ? ptr->isValid() && !ptr->isFake() : false);
    ndof = (ptr ? ptr->ndof() : 0);
    chisq = (ptr ? ptr->chi2() : 0);
    vx = (ptr ? ptr->x() : 0);
    vy = (ptr ? ptr->y() : 0);
    vz = (ptr ? ptr->z() : 0);
    t = (ptr && !use3D ? ptr->t() : 0);
    terr = (ptr && !use3D ? ptr->tError() : 0);
  }
  void addAssociatedTrackInfo(TrackInformation const& info_){ if (!checkListVariable(associatedTrackInfos, &info_)) associatedTrackInfos.push_back(&info_); }
  void addAssociatedMuonInfo(MuonInformation const& info_){ if (!checkListVariable(associatedMuonInfos, &info_)) associatedMuonInfos.push_back(&info_); }

  VertexInformation(reco::Vertex const& ref_, bool use3D_) :
    ptr(&ref_),
    use3D(use3D_)
  { setup(); }
  VertexInformation(reco::Vertex const* ptr_, bool use3D_) :
    ptr(ptr_),
    use3D(use3D_)
  { setup(); }
  VertexInformation() :
    ptr(nullptr)
  { setup(); }

};
struct TrackInformation{
  reco::TrackBaseRef ref;
  edm::ESHandle<TransientTrackBuilder> const* theTTBuilder;
  VertexInformation const* associatedVertex3D;
  VertexInformation const* associatedVertex4D;
  int vtx3DAssociationRank;
  int vtx4DAssociationRank;
  unsigned int nVtx3DWgts;
  unsigned int nVtx4DWgts;
  float px;
  float py;
  float pz;
  float pt;
  float eta;
  float phi;
  float vx;
  float vy;
  float vz;
  float t;
  float terr;

  bool hasTransientTrack;
  reco::TransientTrack ttrk;

  bool isNonnull() const{ return ref.isNonnull(); }
  bool isHighPurity() const{ return ref->quality(reco::TrackBase::highPurity); }
  bool hasTime() const{ return (terr>0.); }
  bool hasAssociatedVertex3D() const{ return (associatedVertex3D!=nullptr); }
  bool hasAssociatedVertex4D() const{ return (associatedVertex4D!=nullptr); }
  void setup(){
    bool refIsNonnull = this->isNonnull();
    px = (refIsNonnull ? ref->px() : 0);
    py = (refIsNonnull ? ref->py() : 0);
    pz = (refIsNonnull ? ref->pz() : 0);
    pt = (refIsNonnull ? ref->pt() : 0);
    eta = (refIsNonnull ? ref->eta() : 0);
    phi = (refIsNonnull ? ref->phi() : 0);
    vx = (refIsNonnull ? ref->vx() : 0);
    vy = (refIsNonnull ? ref->vy() : 0);
    vz = (refIsNonnull ? ref->vz() : 0);
#if _useTrackTime_ == 0
    t = 0;
    terr = 0;
#else
    t = (refIsNonnull ? ref->t0() : 0);
    terr = (refIsNonnull ? ref->tError() : 0);
#endif

    if (theTTBuilder && theTTBuilder->isValid() && refIsNonnull){
      hasTransientTrack=true;
      ttrk = (*theTTBuilder)->build(*ref);
    }
    else hasTransientTrack=false;
  }
  void getImpactPointAndError(VertexInformation const& vtx, GlobalPoint& b, GlobalError& berr) const{
    b = GlobalPoint(0, 0, 0);
    berr = GlobalError(0, 0, 0, 0, 0, 0);

    if (!vtx.ptr) return;

    AnalyticalImpactPointExtrapolator extrapolator(ttrk.field());
    const TrajectoryStateOnSurface tsos = extrapolator.extrapolate(ttrk.impactPointState(), RecoVertex::convertPos(vtx.ptr->position()));
    if (!tsos.isValid()) return;
    else{
      b = tsos.globalPosition();
      berr = tsos.cartesianError().position();
    }
  }

  TrackInformation(reco::TrackBaseRef const& ref_, edm::ESHandle<TransientTrackBuilder> const* theTTBuilder_=nullptr) :
    ref(ref_),
    theTTBuilder(theTTBuilder_),
    associatedVertex3D(nullptr),
    associatedVertex4D(nullptr),
    vtx3DAssociationRank(-1),
    vtx4DAssociationRank(-1),
    nVtx3DWgts(0),
    nVtx4DWgts(0)
  { setup(); }
  TrackInformation() :
    ref(),
    theTTBuilder(nullptr),
    associatedVertex3D(nullptr),
    associatedVertex4D(nullptr),
    vtx3DAssociationRank(-1),
    vtx4DAssociationRank(-1),
    nVtx3DWgts(0),
    nVtx4DWgts(0)
  { setup(); }

};
struct MuonInformation{
  reco::Muon const* ptr;
  TrackInformation const* trkinfo;
  VertexInformation const* associatedVertex3D;
  VertexInformation const* associatedVertex4D;
  int trkIndex;
  float px;
  float py;
  float pz;
  float pt;
  float eta;
  float phi;

  bool isNonnull() const{ return (ptr!=nullptr); }
  bool hasTrackInfo() const{ return (trkinfo!=nullptr); }
  bool hasTime() const{ return (hasTrackInfo() && trkinfo->terr>0.); }

  bool isTrackerMuon() const{ return (ptr && ptr->isTrackerMuon()); }
  bool isGlobalMuon() const{ return (ptr && ptr->isGlobalMuon()); }
  bool isPFMuon() const{ return (ptr && ptr->isPFMuon()); }
  float normalizedChi2() const{ return (ptr && ptr->globalTrack().isNonnull() ? ptr->globalTrack()->normalizedChi2() : -1.); }
  int numberOfValidMuonHits() const{ return (ptr && ptr->globalTrack().isNonnull() ? ptr->globalTrack()->hitPattern().numberOfValidMuonHits() : -1); }
  int numberOfMatchedStations() const{ return (ptr ? ptr->numberOfMatchedStations() : -1); }
  int numberOfValidPixelHits() const{ return (ptr && ptr->innerTrack().isNonnull() ? ptr->innerTrack()->hitPattern().numberOfValidPixelHits() : -1); }
  int trackerLayersWithMeasurement() const{ return (ptr && ptr->innerTrack().isNonnull() ? ptr->innerTrack()->hitPattern().trackerLayersWithMeasurement() : -1); }

  void setup(){
    bool refIsNonnull = this->isNonnull();
    px = (refIsNonnull ? ptr->px() : 0);
    py = (refIsNonnull ? ptr->py() : 0);
    pz = (refIsNonnull ? ptr->pz() : 0);
    pt = (refIsNonnull ? ptr->pt() : 0);
    eta = (refIsNonnull ? ptr->eta() : 0);
    phi = (refIsNonnull ? ptr->phi() : 0);
  }

  MuonInformation(reco::Muon const* ptr_, std::vector<TrackInformation> const& trkinfos) :
    ptr(ptr_),
    trkinfo(nullptr),
    associatedVertex3D(nullptr),
    associatedVertex4D(nullptr),
    trkIndex(-1)
  {
    reco::TrackRef muontrackref = ptr->track();
    reco::TrackBaseRef muontrackbaseref = reco::TrackBaseRef(muontrackref);
    int itrk=0;
    for (TrackInformation const& trkinfo_:trkinfos){
      if (trkinfo_.ref == muontrackbaseref){
        this->trkinfo = &trkinfo_;
        trkIndex=itrk;
        break;
      }
      itrk++;
    }
    setup();
  }
  MuonInformation(reco::Muon const& ref_, std::vector<TrackInformation> const& trkinfos) :
    ptr(&ref_),
    trkinfo(nullptr),
    associatedVertex3D(nullptr),
    associatedVertex4D(nullptr),
    trkIndex(-1)
  {
    reco::TrackRef muontrackref = ptr->get<reco::TrackRef>();
    int itrk=0;
    for (TrackInformation const& trkinfo_:trkinfos){
      if (trkinfo_.ref.castTo<reco::TrackRef>() == muontrackref){
        this->trkinfo = &trkinfo_;
        trkIndex=itrk;
        break;
      }
      itrk++;
    }
    setup();
  }
  MuonInformation() :
    ptr(nullptr),
    trkinfo(nullptr),
    associatedVertex3D(nullptr),
    associatedVertex4D(nullptr),
    trkIndex(-1)
  {
    setup();
  }

};


bool testTrackUsedInVertexFit(VertexInformation const& vtx, TrackInformation const& trk, float* wgt=nullptr){
  float w=0;
  bool res=(vtx.isUsable && trk.isNonnull());
  if (res) w = vtx.ptr->trackWeight(trk.ref);
  res &= (w>0.);
  if (wgt) *wgt=w;
  return res;
}
bool computeIPVals(
  VertexInformation const& vtx, TrackInformation const& trk,
  float& IP, float& d_IP
){
  IP=0; d_IP=0;
  if (vtx.isUsable && trk.isNonnull() && trk.hasTransientTrack){
    std::pair<bool, Measurement1D> IP_Measurement = IPTools::signedImpactParameter3D(
      trk.ttrk,
      GlobalVector(trk.px, trk.py, trk.pz),
      *(vtx.ptr)
    );
    IP = IP_Measurement.second.value();
    d_IP = IP_Measurement.second.error();
    return true;
  }
  else return false;
}
bool computeRelTimeVals(
  VertexInformation const& vtx, TrackInformation const& trk,
  float& dt, float& dterr
){
  dt=0; dterr=0;
  if (vtx.isUsable && vtx.hasTime() && trk.isNonnull() && trk.hasTime()){
    dt = trk.t - vtx.t;
    dterr = sqrt(pow(trk.terr, 2)+pow(vtx.terr, 2));
    return true;
  }
  else return false;
}


//
// class declaration
//
class FTLMuonIsolation : public edm::EDAnalyzer{
public:

  typedef edm::Association<reco::VertexCollection> CandToVertex;
  typedef edm::ValueMap<int> CandToVertexQuality;
  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>, ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;

  explicit FTLMuonIsolation(const edm::ParameterSet&);
  ~FTLMuonIsolation() {};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {};

  //---member data
  edm::EDGetTokenT<genXYZ>                            genXYZToken_;
  edm::Handle<genXYZ>                                 genXYZHandle_;
  edm::EDGetTokenT<float>                             genT0Token_;
  edm::Handle<float>                                  genT0Handle_;
  edm::EDGetTokenT<vector<SimVertex> >                simVtxToken_;
  edm::Handle<vector<SimVertex> >                     simVtxHandle_;
  edm::EDGetTokenT<reco::VertexCollection>            vtx3DToken_;
  edm::Handle<reco::VertexCollection>                 vtx3DHandle_;
  edm::EDGetTokenT<reco::VertexCollection>            vtx4DToken_;
  edm::Handle<reco::VertexCollection>                 vtx4DHandle_;
  edm::EDGetTokenT<reco::MuonCollection>              muonsToken_;
  edm::Handle<reco::MuonCollection>                   muonsHandle_;
  edm::EDGetTokenT<edm::View<reco::Track> >           tracksToken_;
  edm::Handle<edm::View<reco::Track> >                tracksHandle_;
  edm::EDGetTokenT<edm::ValueMap<float> >             timeToken_;
  edm::Handle<edm::ValueMap<float> >                  timeHandle_;
  edm::EDGetTokenT<edm::ValueMap<float> >             timeResToken_;
  edm::Handle<edm::ValueMap<float> >                  timeResHandle_;
  // edm::EDGetTokenT<std::vector<std::vector<float> > > ebtimeToken_;
  // edm::Handle<std::vector<std::vector<float> > >      ebtimeHandle_;    
  edm::EDGetTokenT<reco::GenParticleCollection>       genPartToken_;
  edm::Handle<reco::GenParticleCollection>            genPartHandle_;
  edm::EDGetTokenT<vector<reco::GenJet> >             genJetToken_;
  edm::Handle<vector<reco::GenJet> >                  genJetHandle_;

  //---
  vector<double> targetResolutions_;
  double dzCut_;

  //---I/O
  int iEvent_;
  edm::Service<TFileService> fs;
  map<double, FTLMuonIsoTree> outTrees_;

  //---options
  bool           useMCTruthPV_;
  bool recordTrackInfo_;
  bool recordVertexInfo_;
  double isoConeSize_;
  double isoTimeScale_;

  std::vector<std::pair<int, int>> findDuplicates(const std::vector<float>* fourvector, std::vector<int> id, std::vector<int> status);

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FTLMuonIsolation::FTLMuonIsolation(const edm::ParameterSet& pSet) :
  genXYZToken_(consumes<genXYZ>(pSet.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
  genT0Token_(consumes<float>(pSet.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
  simVtxToken_(consumes<vector<SimVertex> >(pSet.getUntrackedParameter<edm::InputTag>("genVtxTag"))),
  vtx3DToken_(consumes<std::vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtxTag3D"))),
  vtx4DToken_(consumes<std::vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtxTag4D"))),
  muonsToken_(consumes<reco::MuonCollection>(pSet.getUntrackedParameter<edm::InputTag>("muonsTag"))),
  tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("tracksTag"))),
  timeToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("timeTag"))),
  timeResToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("timeResTag"))),
  genPartToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genPartTag"))),
  genJetToken_(consumes<std::vector<reco::GenJet> >(pSet.getUntrackedParameter<edm::InputTag>("genJetsTag"))),
  targetResolutions_(pSet.getUntrackedParameter<vector<double> >("targetResolutions")),
  dzCut_(pSet.getUntrackedParameter<double>("dzCut"))
{
  iEvent_ = 0;
  for (auto& res : targetResolutions_) outTrees_[res] = FTLMuonIsoTree(
    (pSet.getUntrackedParameter<string>("treeName")+"_"+to_string(int(res*1000))).c_str(),
    "Muon tree for FTL studies"
  );
  useMCTruthPV_ = pSet.getUntrackedParameter<bool>("useMCTruthPV");
  recordTrackInfo_ = pSet.getUntrackedParameter<bool>("recordTrackInfo");
  recordVertexInfo_ = pSet.getUntrackedParameter<bool>("recordVertexInfo");
  isoConeSize_ = pSet.getUntrackedParameter<double>("isoConeSize");
  isoTimeScale_ = pSet.getUntrackedParameter<double>("isoTimeScale");
}

//
// member functions
//

// ------------ method called for each event  ------------
void FTLMuonIsolation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //---get input collections
  iEvent.getByToken(genXYZToken_, genXYZHandle_);
  iEvent.getByToken(genT0Token_, genT0Handle_);
  iEvent.getByToken(muonsToken_, muonsHandle_);
  iEvent.getByToken(tracksToken_, tracksHandle_);
  iEvent.getByToken(timeToken_, timeHandle_);
  iEvent.getByToken(timeResToken_, timeResHandle_);
  iEvent.getByToken(simVtxToken_, simVtxHandle_);
  iEvent.getByToken(vtx4DToken_, vtx4DHandle_);
  iEvent.getByToken(vtx3DToken_, vtx3DHandle_);
  iEvent.getByToken(genPartToken_, genPartHandle_);
  iEvent.getByToken(genJetToken_, genJetHandle_);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theTTBuilder);

  ++iEvent_;

  //---get truth PV
  SimVertex simPV;
  if (simVtxHandle_.isValid()) simPV = simVtxHandle_.product()->at(0);
  else{
    auto xyz = genXYZHandle_.product();
    auto t = *genT0Handle_.product();
    auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
    simPV = SimVertex(v, t);
  }

  // Determine the smeared track reco. times
  for (auto& iRes : targetResolutions_){
    //---reset output
    outTrees_[iRes].Reset();

    // Apply basic checks
    if (
      vtx3DHandle_->empty()
      || vtx4DHandle_->empty()
      || tracksHandle_->empty()
      || muonsHandle_->empty()
      ) continue;

    //---fill global info            
    outTrees_[iRes].event = iEvent.id().event();
    outTrees_[iRes].lumi = iEvent.id().luminosityBlock();
    outTrees_[iRes].run = iEvent.id().run();

    //---fill gen vtx info
    outTrees_[iRes].simPVX = simPV.position().x();
    outTrees_[iRes].simPVY = simPV.position().y();
    outTrees_[iRes].simPVZ = simPV.position().z();
    outTrees_[iRes].simPVT = simPV.position().t();

    // Collection of tracks with all the information needed
    std::vector<TrackInformation> trackInfoList;
    trackInfoList.reserve(tracksHandle_->size());
    for (unsigned i = 0; i < tracksHandle_->size(); ++i){
      auto ref = tracksHandle_->refAt(i);
      // Add track to trk info list
      trackInfoList.emplace_back(ref, &theTTBuilder);

      float time = (*timeHandle_)[ref];
      float timeReso = (*timeResHandle_)[ref] != 0.f ? (*timeResHandle_)[ref] : 0.170f;
      float extra_smearing = std::sqrt(iRes*iRes - timeReso*timeReso);
#if _useTrackTime_ == 0
      float t0 = time;
      float t0err = timeReso;
#else
      float t0 = ref->t0();
      float t0err = ref->t0Error();
#endif

      float timeextra = gRandom->Gaus(0., extra_smearing);
      time += timeextra;
      timeReso = std::sqrt(timeReso*timeReso + extra_smearing*extra_smearing);
      if (t0err>0.){
        t0 += timeextra;
        t0err = sqrt(pow(t0err, 2)+pow(extra_smearing, 2));
      }

      // Change time of the tracks to smeared time
      trackInfoList.back().t = t0;
      trackInfoList.back().terr = t0err;
    }

    // Collection of vertices with all the information needed
    std::vector<VertexInformation> vtx3DInfoList;
    vtx3DInfoList.reserve(vtx3DHandle_->size());
    for (unsigned int ivtx=0; ivtx<vtx3DHandle_->size(); ivtx++){
      const auto& vtx = vtx3DHandle_->at(ivtx);
      if (!vtx.isValid() || vtx.isFake()) continue;
      vtx3DInfoList.emplace_back(vtx, true);

      auto const& vtxInfo = vtx3DInfoList.back();
      if (recordVertexInfo_){
        outTrees_[iRes].vtx3D_vx->push_back(vtxInfo.vx);
        outTrees_[iRes].vtx3D_vy->push_back(vtxInfo.vy);
        outTrees_[iRes].vtx3D_vz->push_back(vtxInfo.vz);
        outTrees_[iRes].vtx3D_ndof->push_back(vtxInfo.ndof);
        outTrees_[iRes].vtx3D_chisq->push_back(vtxInfo.chisq);
        outTrees_[iRes].vtx3D_ntrks->push_back(vtxInfo.ptr->nTracks(0.));
      }
    }
    outTrees_[iRes].nVtx3D = vtx3DHandle_->size();

    std::vector<VertexInformation> vtx4DInfoList;
    vtx4DInfoList.reserve(vtx4DHandle_->size());
    for (unsigned int ivtx=0; ivtx<vtx4DHandle_->size(); ivtx++){
      const auto& vtx = vtx4DHandle_->at(ivtx);
      if (!vtx.isValid() || vtx.isFake()) continue;
      vtx4DInfoList.emplace_back(vtx, false);

      auto const& vtxInfo = vtx4DInfoList.back();
      if (recordVertexInfo_){
        outTrees_[iRes].vtx4D_vx->push_back(vtxInfo.vx);
        outTrees_[iRes].vtx4D_vy->push_back(vtxInfo.vy);
        outTrees_[iRes].vtx4D_vz->push_back(vtxInfo.vz);
        outTrees_[iRes].vtx4D_t->push_back(vtxInfo.t);
        outTrees_[iRes].vtx4D_terr->push_back(vtxInfo.terr);
        outTrees_[iRes].vtx4D_ndof->push_back(vtxInfo.ndof);
        outTrees_[iRes].vtx4D_chisq->push_back(vtxInfo.chisq);
        outTrees_[iRes].vtx4D_ntrks->push_back(vtxInfo.ptr->nTracks(0.));
      }
    }
    outTrees_[iRes].nVtx4D = vtx4DHandle_->size();

    // Determine vertex - track associations
    for (TrackInformation& trkInfo : trackInfoList){
      // Loop over 3D vertices
      int chosenVtx3D=-1;
      for (unsigned int ivtx=0; ivtx<vtx3DInfoList.size(); ivtx++){
        auto const& vtxInfo = vtx3DInfoList.at(ivtx);
        if (testTrackUsedInVertexFit(vtxInfo, trkInfo)){
          if (chosenVtx3D<0) chosenVtx3D = ivtx;
          trkInfo.nVtx3DWgts++;
        }
      }
      if (chosenVtx3D<0){
        float minSIP=std::numeric_limits<float>::max();
        for (unsigned int ivtx=0; ivtx<vtx3DInfoList.size(); ivtx++){
          auto const& vtxInfo = vtx3DInfoList.at(ivtx);
          float IP_Vtx = 0, dIP_Vtx = 0;
          if (!computeIPVals(vtxInfo, trkInfo, IP_Vtx, dIP_Vtx)) continue;
          const float valSIP = (dIP_Vtx!=0. ? fabs(IP_Vtx)/dIP_Vtx : 0);
          if (minSIP>valSIP){
            minSIP = valSIP;
            chosenVtx3D = ivtx;
          }
        }
        if (chosenVtx3D>=0) trkInfo.vtx3DAssociationRank = 1;
      }
      else trkInfo.vtx3DAssociationRank = 0;
      if (chosenVtx3D<0) chosenVtx3D=0;
      trkInfo.associatedVertex3D = &(vtx3DInfoList.at(chosenVtx3D));
      vtx3DInfoList.at(chosenVtx3D).addAssociatedTrackInfo(trkInfo);

      // Loop over 4D vertices
      int chosenVtx4D=-1;
      for (unsigned int ivtx=0; ivtx<vtx4DInfoList.size(); ivtx++){
        auto const& vtxInfo = vtx4DInfoList.at(ivtx);
        if (testTrackUsedInVertexFit(vtxInfo, trkInfo)){
          if (chosenVtx4D<0) chosenVtx4D = ivtx;
          trkInfo.nVtx4DWgts++;
        }
      }
      if (chosenVtx4D<0){
        float minSIP=std::numeric_limits<float>::max();
        if (trkInfo.hasTime()){ // Search within vertices with time measurement first if the track has time. Include dt/delta_dt in SIP
          for (unsigned int ivtx=0; ivtx<vtx4DInfoList.size(); ivtx++){
            auto const& vtxInfo = vtx4DInfoList.at(ivtx);
            float IP_Vtx = 0, dIP_Vtx = 0, dt = 0, dterr = 0;
            if (!computeIPVals(vtxInfo, trkInfo, IP_Vtx, dIP_Vtx)) continue;
            if (!computeRelTimeVals(vtxInfo, trkInfo, dt, dterr)) continue;
            const float valSIP = sqrt(pow((dIP_Vtx!=0. ? IP_Vtx/dIP_Vtx : 0.), 2) + pow((dterr!=0. ? dt/dterr : 0.), 2));
            if (minSIP>valSIP){
              minSIP = valSIP;
              chosenVtx4D = ivtx;
            }
          }
        }
        if (chosenVtx4D<0){
          minSIP=std::numeric_limits<float>::max();
          for (unsigned int ivtx=0; ivtx<vtx4DInfoList.size(); ivtx++){
            auto const& vtxInfo = vtx4DInfoList.at(ivtx);
            float IP_Vtx = 0, dIP_Vtx = 0;
            if (!computeIPVals(vtxInfo, trkInfo, IP_Vtx, dIP_Vtx)) continue;
            const float valSIP = (dIP_Vtx!=0. ? fabs(IP_Vtx)/dIP_Vtx : 0);
            if (minSIP>valSIP){
              minSIP = valSIP;
              chosenVtx4D = ivtx;
            }
          }
          if (chosenVtx4D>=0) trkInfo.vtx4DAssociationRank = 2;
        }
        else trkInfo.vtx4DAssociationRank = 1;
      }
      else trkInfo.vtx4DAssociationRank = 0;
      if (chosenVtx4D<0) chosenVtx4D=0;
      trkInfo.associatedVertex4D = &(vtx4DInfoList.at(chosenVtx4D));
      vtx4DInfoList.at(chosenVtx4D).addAssociatedTrackInfo(trkInfo);
    }

    // Loop over the gen. muons
    vector<float> reco_GenParticle_FV[4];
    vector<int> reco_GenParticle_id;
    vector<int> reco_GenParticle_status;
    vector<unsigned int> reco_GenParticle_isPromptFinalState;
    vector<unsigned int> reco_GenParticle_isDirectPromptTauDecayProductFinalState;
    for (reco::GenParticle const& part:(*genPartHandle_)){
      if (
        (std::abs(part.pdgId())==13 && part.status()==1) // Generated muons
        //((part.pdgId()==25 || part.pdgId()==32) && part.status()==22) // Generated Higgs
        //||
        //((std::abs(part.pdgId())>=11 && std::abs(part.pdgId())<=16) && (part.status()==23 || part.status()==1)) // Generated leptons
        //||
        //((std::abs(part.pdgId())<=6 || std::abs(part.pdgId())==21) && (part.status()==21 || part.status()==23 || part.status()==22)) // Generated partons
        ){
        reco_GenParticle_FV[0].push_back(part.px());
        reco_GenParticle_FV[1].push_back(part.py());
        reco_GenParticle_FV[2].push_back(part.pz());
        reco_GenParticle_FV[3].push_back(part.energy());
        reco_GenParticle_id.push_back(part.pdgId());
        reco_GenParticle_status.push_back(part.status());
        reco_GenParticle_isPromptFinalState.push_back(part.isPromptFinalState());
        reco_GenParticle_isDirectPromptTauDecayProductFinalState.push_back(part.isDirectPromptTauDecayProductFinalState());
      }
    }
    outTrees_[iRes].genmuon_px->assign(reco_GenParticle_FV[0].begin(), reco_GenParticle_FV[0].end());
    outTrees_[iRes].genmuon_py->assign(reco_GenParticle_FV[1].begin(), reco_GenParticle_FV[1].end());
    outTrees_[iRes].genmuon_pz->assign(reco_GenParticle_FV[2].begin(), reco_GenParticle_FV[2].end());
    outTrees_[iRes].genmuon_E->assign(reco_GenParticle_FV[3].begin(), reco_GenParticle_FV[3].end());
    outTrees_[iRes].genmuon_id->assign(reco_GenParticle_id.begin(), reco_GenParticle_id.end());
    outTrees_[iRes].genmuon_status->assign(reco_GenParticle_status.begin(), reco_GenParticle_status.end());
    outTrees_[iRes].genmuon_isPromptFinalState->assign(reco_GenParticle_isPromptFinalState.begin(), reco_GenParticle_isPromptFinalState.end());
    outTrees_[iRes].genmuon_isDirectPromptTauDecayProductFinalState->assign(reco_GenParticle_isDirectPromptTauDecayProductFinalState.begin(), reco_GenParticle_isDirectPromptTauDecayProductFinalState.end());


    // Loop over reco. muons
    std::vector<MuonInformation> muonInfoList; muonInfoList.reserve(muonsHandle_->size());
    for (auto const& muon:(*muonsHandle_)){
      //---basic check
      if (muon.track().isNull()) continue;
      muonInfoList.emplace_back(muon, trackInfoList);
      auto& muonInfo = muonInfoList.back();
      TrackInformation const* trkInfo = muonInfo.trkinfo;
      if (!trkInfo) continue;

      // Vertex association variables
      float IP3DVtx3D_0=0, dIP3DVtx3D_0=0;
      float IP3DVtx4D_0=0, dIP3DVtx4D_0=0;

      int chosenVtx3D=-1;
      //int rankVtx3D=-1;
      float IP3DVtx3D=0, dIP3DVtx3D=0;

      int chosenVtx4D=-1;
      //int rankVtx4D=-1;
      float IP3DVtx4D=0, dIP3DVtx4D=0;

      // chIso pt sums
      float muon_isosumtrackpt_vtx3D=0;
      float muon_isosumtrackpt_vtx3D_dzcut=0;

      float muon_isosumtrackpt_vtx3D_nodtmuon=0;
      float muon_isosumtrackpt_vtx3D_dzcut_nodtmuon=0;

      float muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut=0;
      float muon_isosumtrackpt_vtx3D_dzcut_hasdtmuon_dtmuoncut=0;

      float muon_isosumtrackpt_vtx4D=0;
      float muon_isosumtrackpt_vtx4D_dzcut=0;


      // Loop over 3D vertices
      for (unsigned int ivtx=0; ivtx<vtx3DInfoList.size(); ivtx++){
        auto const& vtxInfo = vtx3DInfoList.at(ivtx);
        if (testTrackUsedInVertexFit(vtxInfo, *trkInfo)) chosenVtx3D = ivtx;
      }
      if (chosenVtx3D<0){
        float minSIP=std::numeric_limits<float>::max();
        for (unsigned int ivtx=0; ivtx<vtx3DInfoList.size(); ivtx++){
          auto const& vtxInfo = vtx3DInfoList.at(ivtx);
          float IP_Vtx = 0, dIP_Vtx = 0;
          if (!computeIPVals(vtxInfo, *trkInfo, IP_Vtx, dIP_Vtx)) continue;
          const float valSIP = (dIP_Vtx!=0. ? fabs(IP_Vtx)/dIP_Vtx : 0);
          if (minSIP>valSIP){
            minSIP = valSIP;
            chosenVtx3D = ivtx;
          }
        }
        //if (chosenVtx3D>=0) rankVtx3D=1;
      }
      //else rankVtx3D=0;
      if (chosenVtx3D<0) chosenVtx3D=0;
      muonInfo.associatedVertex3D = &(vtx3DInfoList.at(chosenVtx3D));
      computeIPVals(vtx3DInfoList.at(chosenVtx3D), *trkInfo, IP3DVtx3D, dIP3DVtx3D);
      computeIPVals(vtx3DInfoList.at(0), *trkInfo, IP3DVtx3D_0, dIP3DVtx3D_0);

      // Loop over 4D vertices
      for (unsigned int ivtx=0; ivtx<vtx4DInfoList.size(); ivtx++){
        auto const& vtxInfo = vtx4DInfoList.at(ivtx);
        if (testTrackUsedInVertexFit(vtxInfo, *trkInfo)) chosenVtx4D = ivtx;
      }
      if (chosenVtx4D<0){
        float minSIP=std::numeric_limits<float>::max();
        if (trkInfo->hasTime()){ // Search within vertices with time measurement first if the track has time. Include dt/delta_dt in SIP
          for (unsigned int ivtx=0; ivtx<vtx4DInfoList.size(); ivtx++){
            auto const& vtxInfo = vtx4DInfoList.at(ivtx);
            float IP_Vtx = 0, dIP_Vtx = 0, dt = 0, dterr = 0;
            if (!computeIPVals(vtxInfo, *trkInfo, IP_Vtx, dIP_Vtx)) continue;
            if (!computeRelTimeVals(vtxInfo, *trkInfo, dt, dterr)) continue;
            const float valSIP = sqrt(pow((dIP_Vtx!=0. ? IP_Vtx/dIP_Vtx : 0.), 2) + pow((dterr!=0. ? dt/dterr : 0.), 2));
            if (minSIP>valSIP){
              minSIP = valSIP;
              chosenVtx4D = ivtx;
            }
          }
        }
        if (chosenVtx4D<0){
          minSIP=std::numeric_limits<float>::max();
          for (unsigned int ivtx=0; ivtx<vtx4DInfoList.size(); ivtx++){
            auto const& vtxInfo = vtx4DInfoList.at(ivtx);
            float IP_Vtx = 0, dIP_Vtx = 0;
            if (!computeIPVals(vtxInfo, *trkInfo, IP_Vtx, dIP_Vtx)) continue;
            const float valSIP = (dIP_Vtx!=0. ? fabs(IP_Vtx)/dIP_Vtx : 0);
            if (minSIP>valSIP){
              minSIP = valSIP;
              chosenVtx4D = ivtx;
            }
          }
          //if (chosenVtx4D>=0) rankVtx4D=2;
        }
        //else rankVtx4D=1;
      }
      //else rankVtx4D=0;
      if (chosenVtx4D<0) chosenVtx4D=0;
      muonInfo.associatedVertex4D = &(vtx4DInfoList.at(chosenVtx4D));
      computeIPVals(vtx4DInfoList.at(chosenVtx4D), *trkInfo, IP3DVtx4D, dIP3DVtx4D);

      // Compute isolation sums
      for (TrackInformation const& trkInfo_ : trackInfoList){
        if (trkInfo_.ref == trkInfo->ref) continue;
        if (!trkInfo_.isHighPurity()) continue;

        float this_dr = reco::deltaR2(trkInfo_.ref->eta(), trkInfo_.ref->phi(), muonInfo.eta, muonInfo.phi);
        bool keep_dr = (this_dr <= pow(isoConeSize_, 2));
        if (!keep_dr) continue;

        // Checks on 3D vertex-like association
        GlobalPoint b3D;
        GlobalError berr3D;
        trkInfo_.getImpactPointAndError(*(muonInfo.associatedVertex3D), b3D, berr3D);
        float dz = std::abs(b3D.z() - muonInfo.associatedVertex3D->vz);
        bool keep_dz = (dz <= dzCut_);

        bool has_dtmuon = (trkInfo_.hasTime() && muonInfo.hasTime());
        float dtmuon = (has_dtmuon ? std::abs(trkInfo_.t - trkInfo->t) : 0);
        float dtmuonerr = (has_dtmuon ? sqrt(pow(trkInfo_.terr, 2) + pow(trkInfo->terr, 2)) : 0);
        float reldtmuon = (dtmuonerr>0. ? dtmuon/dtmuonerr : 0.);
        bool keep_dtmuon = (reldtmuon <= isoTimeScale_);

        if (keep_dr){
          muon_isosumtrackpt_vtx3D += trkInfo_.pt;
          if (keep_dz) muon_isosumtrackpt_vtx3D_dzcut += trkInfo_.pt;

          if (keep_dtmuon){
            if (has_dtmuon){
              muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut += trkInfo_.pt;
              if (keep_dz) muon_isosumtrackpt_vtx3D_dzcut_hasdtmuon_dtmuoncut += trkInfo_.pt;
            }
            else{
              muon_isosumtrackpt_vtx3D_nodtmuon += trkInfo_.pt;
              if (keep_dz) muon_isosumtrackpt_vtx3D_dzcut_nodtmuon += trkInfo_.pt;
            }
          }
        }

        // Checks on 4D vertex-like association
        GlobalPoint b4D;
        GlobalError berr4D;
        trkInfo_.getImpactPointAndError(*(muonInfo.associatedVertex4D), b4D, berr4D);
        dz = std::abs(b4D.z() - muonInfo.associatedVertex4D->vz);
        keep_dz = (dz <= dzCut_);

        float dt = (trkInfo_.hasTime() && muonInfo.associatedVertex4D->hasTime() ? std::abs(trkInfo_.t - muonInfo.associatedVertex4D->t) : 0);
        float dterr = (trkInfo_.hasTime() && muonInfo.associatedVertex4D->hasTime() ? sqrt(pow(trkInfo_.terr, 2) + pow(muonInfo.associatedVertex4D->terr, 2)) : 0);
        float reldt = (dterr>0. ? dt/dterr : 0.);
        bool keep_dt = (reldt <= isoTimeScale_);

        if (keep_dr && keep_dt){
          muon_isosumtrackpt_vtx4D += trkInfo_.pt;
          if (keep_dz) muon_isosumtrackpt_vtx4D_dzcut += trkInfo_.pt;
        }
      }
      
      //---fill muon info
      outTrees_[iRes].muon_pt->push_back(muonInfo.pt);
      outTrees_[iRes].muon_eta->push_back(muonInfo.eta);
      outTrees_[iRes].muon_phi->push_back(muonInfo.phi);
      outTrees_[iRes].muon_px->push_back(muonInfo.px);
      outTrees_[iRes].muon_py->push_back(muonInfo.py);
      outTrees_[iRes].muon_pz->push_back(muonInfo.pz);

      outTrees_[iRes].muon_vx->push_back(trkInfo->vx);
      outTrees_[iRes].muon_vy->push_back(trkInfo->vy);
      outTrees_[iRes].muon_vz->push_back(trkInfo->vz);
      outTrees_[iRes].muon_t->push_back(trkInfo->t);
      outTrees_[iRes].muon_terr->push_back(trkInfo->terr);
      outTrees_[iRes].isLooseMuon->push_back(muon::isLooseMuon(muon));
      outTrees_[iRes].isMediumMuon->push_back(muon::isMediumMuon(muon));
      outTrees_[iRes].isTightMuon3D->push_back(muon::isTightMuon(muon, *(vtx3DInfoList.at(chosenVtx3D).ptr)));
      outTrees_[iRes].isTightMuon4D->push_back(muon::isTightMuon(muon, *(vtx4DInfoList.at(chosenVtx4D).ptr)));

      outTrees_[iRes].isTrackerMuon->push_back(muonInfo.isTrackerMuon());
      outTrees_[iRes].isGlobalMuon->push_back(muonInfo.isGlobalMuon());
      outTrees_[iRes].isPFMuon->push_back(muonInfo.isPFMuon());
      outTrees_[iRes].normalizedChi2->push_back(muonInfo.normalizedChi2());
      outTrees_[iRes].numberOfValidMuonHits->push_back(muonInfo.numberOfValidMuonHits());
      outTrees_[iRes].numberOfMatchedStations->push_back(muonInfo.numberOfMatchedStations());
      outTrees_[iRes].numberOfValidPixelHits->push_back(muonInfo.numberOfValidPixelHits());
      outTrees_[iRes].trackerLayersWithMeasurement->push_back(muonInfo.trackerLayersWithMeasurement());
      outTrees_[iRes].muonTrkId->push_back(muonInfo.trkIndex);
      outTrees_[iRes].muonVtx3DId->push_back(chosenVtx3D);
      outTrees_[iRes].muonIP3DVtx3D->push_back(IP3DVtx3D);
      outTrees_[iRes].muondIP3DVtx3D->push_back(dIP3DVtx3D);
      outTrees_[iRes].muonIP3DVtx3D_0->push_back(IP3DVtx3D_0);
      outTrees_[iRes].muondIP3DVtx3D_0->push_back(dIP3DVtx3D_0);
      //outTrees_[iRes].muonVtx3DAssociationRank->push_back(rankVtx3D);
      outTrees_[iRes].muon_Vtx3D_vx->push_back(muonInfo.associatedVertex3D->vx);
      outTrees_[iRes].muon_Vtx3D_vy->push_back(muonInfo.associatedVertex3D->vy);
      outTrees_[iRes].muon_Vtx3D_vz->push_back(muonInfo.associatedVertex3D->vz);
      outTrees_[iRes].muon_Vtx3D_ndof->push_back(muonInfo.associatedVertex3D->ndof);
      outTrees_[iRes].muon_Vtx3D_chisq->push_back(muonInfo.associatedVertex3D->chisq);
      outTrees_[iRes].muon_Vtx3D_ntrks->push_back(muonInfo.associatedVertex3D->ptr->nTracks(0.));

      outTrees_[iRes].muonVtx4DId->push_back(chosenVtx4D);
      outTrees_[iRes].muonIP3DVtx4D->push_back(IP3DVtx4D);
      outTrees_[iRes].muondIP3DVtx4D->push_back(dIP3DVtx4D);
      outTrees_[iRes].muonIP3DVtx4D_0->push_back(IP3DVtx4D_0);
      outTrees_[iRes].muondIP3DVtx4D_0->push_back(dIP3DVtx4D_0);
      //outTrees_[iRes].muonVtx4DAssociationRank->push_back(rankVtx4D);
      outTrees_[iRes].muon_Vtx4D_vx->push_back(muonInfo.associatedVertex4D->vx);
      outTrees_[iRes].muon_Vtx4D_vy->push_back(muonInfo.associatedVertex4D->vy);
      outTrees_[iRes].muon_Vtx4D_vz->push_back(muonInfo.associatedVertex4D->vz);
      outTrees_[iRes].muon_Vtx4D_t->push_back(muonInfo.associatedVertex4D->t);
      outTrees_[iRes].muon_Vtx4D_terr->push_back(muonInfo.associatedVertex4D->terr);
      outTrees_[iRes].muon_Vtx4D_ndof->push_back(muonInfo.associatedVertex4D->ndof);
      outTrees_[iRes].muon_Vtx4D_chisq->push_back(muonInfo.associatedVertex4D->chisq);
      outTrees_[iRes].muon_Vtx4D_ntrks->push_back(muonInfo.associatedVertex4D->ptr->nTracks(0.));

      outTrees_[iRes].muon_isosumtrackpt_vtx3D->push_back(muon_isosumtrackpt_vtx3D);
      outTrees_[iRes].muon_isosumtrackpt_vtx3D_dzcut->push_back(muon_isosumtrackpt_vtx3D_dzcut);

      outTrees_[iRes].muon_isosumtrackpt_vtx3D_nodtmuon->push_back(muon_isosumtrackpt_vtx3D_nodtmuon);
      outTrees_[iRes].muon_isosumtrackpt_vtx3D_dzcut_nodtmuon->push_back(muon_isosumtrackpt_vtx3D_dzcut_nodtmuon);

      outTrees_[iRes].muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut->push_back(muon_isosumtrackpt_vtx3D_hasdtmuon_dtmuoncut);
      outTrees_[iRes].muon_isosumtrackpt_vtx3D_dzcut_hasdtmuon_dtmuoncut->push_back(muon_isosumtrackpt_vtx3D_dzcut_hasdtmuon_dtmuoncut);

      outTrees_[iRes].muon_isosumtrackpt_vtx4D->push_back(muon_isosumtrackpt_vtx4D);
      outTrees_[iRes].muon_isosumtrackpt_vtx4D_dzcut->push_back(muon_isosumtrackpt_vtx4D_dzcut);
    }
    outTrees_[iRes].nMuons = muonInfoList.size();

    //---fill trk info
    if (recordTrackInfo_){
      for (auto const& trackInfo:trackInfoList){
        outTrees_[iRes].track_pt->push_back(trackInfo.pt);
        outTrees_[iRes].track_eta->push_back(trackInfo.eta);
        outTrees_[iRes].track_phi->push_back(trackInfo.phi);
        outTrees_[iRes].track_px->push_back(trackInfo.px);
        outTrees_[iRes].track_py->push_back(trackInfo.py);
        outTrees_[iRes].track_pz->push_back(trackInfo.pz);
        outTrees_[iRes].track_vx->push_back(trackInfo.vx);
        outTrees_[iRes].track_vy->push_back(trackInfo.vy);
        outTrees_[iRes].track_vz->push_back(trackInfo.vz);
        outTrees_[iRes].track_t->push_back(trackInfo.t);
        outTrees_[iRes].track_terr->push_back(trackInfo.terr);

        int trackVtx3DId=0;
        for (auto const& vtxInfo:vtx3DInfoList){
          if (trackInfo.associatedVertex3D==&vtxInfo){
            outTrees_[iRes].trackVtx3DId->push_back(trackVtx3DId);
            float IP3DVtx3D, dIP3DVtx3D;
            computeIPVals(vtxInfo, trackInfo, IP3DVtx3D, dIP3DVtx3D);
            outTrees_[iRes].trackIP3DVtx3D->push_back(IP3DVtx3D);
            outTrees_[iRes].trackdIP3DVtx3D->push_back(dIP3DVtx3D);
            break;
          }
          trackVtx3DId++;
        }
        outTrees_[iRes].trackVtx3DAssociationRank->push_back(trackInfo.vtx3DAssociationRank);
        outTrees_[iRes].trackNVtx3DWgts->push_back(trackInfo.nVtx3DWgts);

        int trackVtx4DId=0;
        for (auto const& vtxInfo:vtx4DInfoList){
          if (trackInfo.associatedVertex4D==&vtxInfo){
            outTrees_[iRes].trackVtx4DId->push_back(trackVtx4DId);
            float IP3DVtx4D, dIP3DVtx4D;
            computeIPVals(vtxInfo, trackInfo, IP3DVtx4D, dIP3DVtx4D);
            outTrees_[iRes].trackIP3DVtx4D->push_back(IP3DVtx4D);
            outTrees_[iRes].trackdIP3DVtx4D->push_back(dIP3DVtx4D);
            break;
          }
          trackVtx4DId++;
        }
        outTrees_[iRes].trackVtx4DAssociationRank->push_back(trackInfo.vtx4DAssociationRank);
        outTrees_[iRes].trackNVtx4DWgts->push_back(trackInfo.nVtx4DWgts);
      }
    }

    //---Fill tree
    outTrees_[iRes].GetTTreePtr()->Fill();
  }
}
  
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FTLMuonIsolation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

std::vector<std::pair<int, int>> FTLMuonIsolation::findDuplicates(const std::vector<float>* fourvector, std::vector<int> id, std::vector<int> status){
  vector<pair<int, int>> duplicates;
  if (id.size()>1) duplicates.reserve(id.size()/2);

  for (unsigned int xx=0; xx<id.size(); xx++){
    TLorentzVector p_tm(fourvector[0].at(xx), fourvector[1].at(xx), fourvector[2].at(xx), fourvector[3].at(xx));
    int id_tm = id.at(xx);
    int st_tm = status.at(xx);

    if (st_tm==23 || st_tm==1){
      double msq_tm = p_tm.M2();
      for (unsigned int yy=xx+1; yy<id.size(); yy++){
        TLorentzVector p_tbm(fourvector[0].at(yy), fourvector[1].at(yy), fourvector[2].at(yy), fourvector[3].at(yy));
        int id_tbm = id.at(yy);
        int st_tbm = status.at(yy);

        if (id_tbm==id_tm && st_tbm!=st_tm && (st_tbm==1 || st_tbm==23)){
          double dot_tm = p_tbm.Dot(p_tm);
          double diff_sqmass = dot_tm - msq_tm;
          diff_sqmass = std::abs(diff_sqmass);
          if (diff_sqmass<0.005*std::abs(msq_tm)){
            int iFirst = (st_tm==23 ? xx : yy);
            int iSecond = (st_tm==23 ? yy : xx);
            duplicates.emplace_back(iFirst, iSecond);
          }
        }
      }
    }
  }

  return duplicates;
}

//define this as a plug-in
DEFINE_FWK_MODULE(FTLMuonIsolation);

#endif
