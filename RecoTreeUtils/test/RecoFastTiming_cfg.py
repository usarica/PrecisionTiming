import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from FastTiming.RecoTreeUtils.RecoFastTiming_cfi import *

options = VarParsing('analysis')

options.register('sampleName',
                 'SingleGammaE50_noPU',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "sample to process")
options.parseArguments()

filesOpt = cms.PSet(
    inputFiles = cms.untracked.vstring(""),
    outputFile = cms.string("")
)

GetSampleFiles(options.sampleName, filesOpt)

process = cms.Process("RecoFastTiming")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = filesOpt.inputFiles)

process.TFileService = cms.Service("TFileService",
                            fileName = filesOpt.outputFile)

process.ft = cms.EDAnalyzer('RecoFastTiming')

process.p = cms.Path(process.ft)
