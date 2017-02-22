import FWCore.ParameterSet.Config as cms

FTLDumpPhotonsRECO = cms.EDAnalyzer(
    "FTLDumpPhotonsRECO",
    photonsTag = cms.untracked.InputTag("gedPhotons", "", "RECO"),
    simTkTag = cms.untracked.InputTag("g4SimHits", "", "HLT"),
    simVtxTag = cms.untracked.InputTag("g4SimHits", "", "HLT"),
    mcTruthPhoEtThr = cms.untracked.double(10),
    treeName = cms.untracked.string("pho_tree")
)

FTLDumpPhotonsPAT = cms.EDAnalyzer(
    "FTLDumpPhotonsPAT",
    photonsTag = cms.untracked.InputTag("slimmedPhotons", "", "RECO"),
    simTkTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    simVtxTag = cms.untracked.InputTag("g4SimHits", "", "SIM"),
    mcTruthPhoEtThr = cms.untracked.double(10),
    treeName = cms.untracked.string("pho_tree")
)
