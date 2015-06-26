import FWCore.ParameterSet.Config as cms

RecoFastTiming = cms.EDAnalyzer("RecoFastTiming",
    genVtxTag       = cms.untracked.InputTag("mix", "InitialVertices", "DIGI2RAW"),
    ghostVtxTag     = cms.untracked.InputTag("VtxSmeared", "", "RecoFastTiming"),
    genJetsTag      = cms.untracked.InputTag('ak5GenJets', '', 'SIM'),
    recoJetsTag     = cms.untracked.InputTag('ak5PFJetsCHS', '', 'RECO'),
    ## add fake PU vtxs
    makeGhosts      = cms.untracked.bool(False),
    ## raw time smearing in ns
    timeResSmearing = cms.untracked.double(0.03),
    ## cut on the impact parameter (assing track to vertex) in cm
    dzCut           = cms.untracked.double(0.2),
    ## track pt and pz threshold for vtx time reconstrution in GeV
    ptCut           = cms.untracked.double(2),
    pz2Cut          = cms.untracked.double(0),
    ## output file collections option
    saveParticles   = cms.untracked.bool(True),
    saveAllRecHits  = cms.untracked.bool(False),
    saveVertices    = cms.untracked.bool(True),
    readMVA         = cms.untracked.bool(False)
)
