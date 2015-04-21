import FWCore.ParameterSet.Config as cms

BDTInputGenerator = cms.EDAnalyzer("BDTInputGenerator",
    ## analyzed particles type: photons
    particleType = cms.untracked.int32(4),
    ## matrix side
    matrixSide = cms.untracked.int32(3)
)
