import FWCore.ParameterSet.Config as cms

FTLDumpPhotons = cms.EDAnalyzer(
    "FTLDumpPhotons",
    photonsTag = cms.untracked.InputTag("slimmedPhotons", "", "RECO"),
    simTkTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    simVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    mcTruthPhoEtThr = cms.untracked.double(10),
    treeName = cms.untracked.string("pho_tree")
)
