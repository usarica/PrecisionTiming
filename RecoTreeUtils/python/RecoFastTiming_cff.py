import FWCore.ParameterSet.Config as cms

RecoFastTiming = cms.EDAnalyzer("RecoFastTiming",
    ## raw time smearing in ns
    timeResSmearing = cms.untracked.double(0.03),
    ## cut on the impact parameter (assing track to vertex) in cm
    dzCut           = cms.untracked.double(0.2),
    ## track pt threshold for vtx time reconstrution in GeV
    ptCut           = cms.untracked.double(2),
    ## output file collections option
    saveParticles   = cms.untracked.bool(True),
    saveAllRecHits  = cms.untracked.bool(False),
    saveVertices    = cms.untracked.bool(True)
)
