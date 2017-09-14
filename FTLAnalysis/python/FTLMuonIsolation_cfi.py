import FWCore.ParameterSet.Config as cms

FTLMuonIsolation = cms.EDAnalyzer(
    "FTLMuonIsolation",
    ###---Input tags
    genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
    genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),	    
    muonsTag = cms.untracked.InputTag("muons", "", "RECO"),	
    tracksTag = cms.untracked.InputTag("generalTracks", "", "RECO"),
    timeTag = cms.untracked.InputTag("trackTimeValueMapProducer",
                                     "generalTracksConfigurableFlatResolutionModel", "RECO"), 
    timeResTag = cms.untracked.InputTag("trackTimeValueMapProducer",
                                          "generalTracksConfigurableFlatResolutionModelResolution", "RECO"),
    genVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    vtxTag4D = cms.untracked.InputTag("offlinePrimaryVertices", "", "RECO"),
    vtxTag3D = cms.untracked.InputTag("offlinePrimaryVertices1D", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    ###---Target time resolution (assumes sample were made with 20ps track t resolution)
    targetResolutions = cms.untracked.vdouble(0.02, 0.03, 0.05, 0.07, 0.09),
    ###---I/O options
    treeName = cms.untracked.string("muon_tree"),
    ###---Vtx choice option
    useMCTruthPV = cms.untracked.bool(False),
    ###---Iso options
    isoConeSizes = cms.untracked.vdouble(0.3, 0.4)
)
