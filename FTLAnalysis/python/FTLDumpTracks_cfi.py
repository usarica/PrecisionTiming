import FWCore.ParameterSet.Config as cms

FTLDumpTracks = cms.EDAnalyzer(
    "FTLDumpTracks",
    genParticlesTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    tracksTag = cms.untracked.InputTag("generalTracks", "", "RECO"),
    ftlRecHitsTag = cms.untracked.InputTag("ftlUncalibratedRecHits", "FTLBarrel", "RECO"),
    # simTkTag = cms.untracked.InputTag("g4SimHits", "", "HLT"),
    # simVtxTag = cms.untracked.InputTag("g4SimHits", "", "HLT"),
    treeName = cms.untracked.string("track_tree")
)
