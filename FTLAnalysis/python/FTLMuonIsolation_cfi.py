import FWCore.ParameterSet.Config as cms

FTLMuonIsolation = cms.EDAnalyzer(
    "FTLMuonIsolation",
    ###---Input tags
    muonsTag = cms.untracked.InputTag("muons", "", "RECO"),	
    tracksTag = cms.untracked.InputTag("generalTracks", "", "RECO"),
    timeTag = cms.untracked.InputTag("trackTimeValueMapProducer",
                                     "generalTracksConfigurableFlatResolutionModel", "RECO"), 
    timeResTag = cms.untracked.InputTag("trackTimeValueMapProducer",
                                          "generalTracksConfigurableFlatResolutionModelResolution", "RECO"),
    vtxTag = cms.untracked.InputTag("offlinePrimaryVertices4D", "", "RECO"),
    genPartTag = cms.untracked.InputTag("genParticles", "", "HLT"),
    genJetsTag = cms.untracked.InputTag("ak4GenJets", "", "HLT"),
    ###---I/O options
    treeName = cms.untracked.string("muon_tree"),
    ###---Iso options
    isoConeSizes = cms.untracked.vdouble(0.3, 0.4)
)
