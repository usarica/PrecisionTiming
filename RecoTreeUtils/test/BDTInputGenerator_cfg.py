import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag
from FastTiming.RecoTreeUtils.IOFilesHelper import *

options = VarParsing('analysis')

options.register('sampleName',
                 'SingleCharged',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "sample to process")
options.register('nameMod',
                 '_BDTInput',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "output file modifier")
options.maxEvents = -1
options.parseArguments()

## Get I/O files from the list given the sample name
filesOpt = cms.PSet(
    inputFiles = cms.untracked.vstring(""),
    outputFile = cms.string("")
)

GetSampleFiles(options.sampleName, options.nameMod, filesOpt)

##------------------------------------------------------------------------------

process = cms.Process("BDTInputGenerator")

## load the SK geometry and magnetic field config
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
## import standard RecoFT configurations
process.load("FastTiming.RecoTreeUtils.BDTInputGenerator_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
                            fileNames = filesOpt.inputFiles)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService",
                            fileName = filesOpt.outputFile)

process.ft_path = cms.Sequence(process.BDTInputGenerator)

process.path = cms.Path(process.ft_path)
