import FWCore.ParameterSet.Config as cms

FTLDumper = cms.EDAnalyzer(
    "FTLDumper",
    electronsTag = cms.untracked.InputTag("slimmedElectrons", "", "RECO")
)
