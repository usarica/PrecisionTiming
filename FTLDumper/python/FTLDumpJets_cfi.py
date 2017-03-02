import FWCore.ParameterSet.Config as cms

FTLDumpJets = cms.EDAnalyzer(
    "FTLDumpJets",
    jetsTag = cms.untracked.InputTag("ak4PFJets", "", "RECO"),
    photonsTag = cms.untracked.InputTag("gedPhotons", "", "RECO"),
    ftlRecHitsTag = cms.untracked.InputTag("ftlRecHits", "FTLBarrel", "RECO"),
    simTkTag = cms.untracked.InputTag("g4SimHits", "", "HLT"),
    simVtxTag = cms.untracked.InputTag("g4SimHits", "", "HLT"),
    mcTruthPhoEtThr = cms.untracked.double(10),
    readFTLRecHits = cms.untracked.bool(True),
    treeName = cms.untracked.string("jet_tree")
)
