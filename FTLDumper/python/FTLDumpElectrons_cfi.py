import FWCore.ParameterSet.Config as cms

FTLDumpElectrons = cms.EDAnalyzer(
    "FTLDumpElectrons",
    electronsTag = cms.untracked.InputTag("slimmedElectrons", "", "RECO"),
    treeName = cms.untracked.string("ele_tree")
)
