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

## Get I/O files from the list given the sample name
filesOpt = cms.PSet(
    inputFiles = cms.untracked.vstring(""),
    outputFile = cms.string("")
)

GetSampleFiles(options.sampleName, filesOpt)

##------------------------------------------------------------------------------

process = cms.Process("RecoFastTiming")

## load the SK geometry and magnetic field config
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaper_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
                            fileNames = filesOpt.inputFiles)

process.TFileService = cms.Service("TFileService",
                            fileName = filesOpt.outputFile)

process.ft = cms.EDAnalyzer('RecoFastTiming')

process.p = cms.Path(process.ft)
