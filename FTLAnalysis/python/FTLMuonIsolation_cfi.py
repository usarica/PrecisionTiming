import FWCore.ParameterSet.Config as cms

FTLMuonIsolation = cms.EDAnalyzer(
    "FTLMuonIsolation",
    ###---Input tags
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),	    
    muonsTag = cms.untracked.InputTag("muons", "", "RECO"),	
    tracksTag = cms.untracked.InputTag("particleFlow", "", "RECO"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIMI"),
    vtxTag = cms.untracked.InputTag("offlinePrimaryVertices", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    ###---Target time resolution (assumes sample were made with 30ps track t resolution)
    targetResolutions = cms.untracked.vdouble(0.03, 0.05, 0.07, 0.09, 0.15),
    ###---I/O options
    treeName = cms.untracked.string("muon_tree"),
    ###---Vtx choice option
    useMCTruthPV = cms.untracked.bool(False),
    ###---precessing type w/ or w/o timing
    isTimingSample = cms.untracked.bool(False),
    ###---save tracks info
    saveTracksInfo = cms.untracked.bool(False),
    ###---no ETL, get time in endcap from parametrized HGC response
    HGCToySim = cms.untracked.bool(False),
    ###---Trk-vtx dz cut
    dzCut = cms.untracked.double(0.1),
    ###---Iso options
    isoConeSizes = cms.untracked.vdouble(0.3)
)

FTLMuonIsolationHGCToy = cms.EDAnalyzer(
    "FTLMuonIsolationHGC",
    ###---Input tags
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),	    
    muonsTag = cms.untracked.InputTag("muons", "", "RECO"),	
    tracksTag = cms.untracked.InputTag("particleFlow", "", "RECO"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    vtxTag = cms.untracked.InputTag("offlinePrimaryVertices", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    ###---Target time resolution (assumes sample were made with 30ps track t resolution)
    targetResolutions = cms.untracked.vdouble(0.03), #ns
    ###---I/O options
    treeName = cms.untracked.string("muon_tree"),
    ###---Vtx choice option
    useMCTruthPV = cms.untracked.bool(False),
    ###---precessing type w/ or w/o timing
    isTimingSample = cms.untracked.bool(False),
    ###---Trk-vtx dz cut
    dzCut = cms.untracked.double(0.1),
    ###---Iso options
    isoConeSizes = cms.untracked.vdouble(0.3, 0.4),
    ###---Time resolution eta and pt bins
    absEtaBins = cms.untracked.VPSet(
        cms.PSet(
            absEtaMax = cms.untracked.double(1.5),
            ptBins = cms.untracked.vdouble(2.0, 7000.),
            timeRes = cms.untracked.vdouble(0.07, 0.03) #ns
        ),
        cms.PSet(
            absEtaMax = cms.untracked.double(4.0),
            ptBins = cms.untracked.vdouble(1.0, 7000.),
            timeRes = cms.untracked.vdouble(-1, 0.03) #ns
        )
    )
)

