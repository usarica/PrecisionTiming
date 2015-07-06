import time
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
from FastTiming.RecoTreeUtils.IOFilesHelper import *

options = VarParsing('analysis')

options.register('sampleName',
                 'SingleGammaE50_noPU',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "sample to process")
options.maxEvents = -1
options.parseArguments()

## Get I/O files from the list given the sample name
filesOpt = cms.PSet(
    inputFiles = cms.untracked.vstring(""),
    outputFile = cms.string("")
)

GetSampleFiles(options.sampleName, "", filesOpt)

##------------------------------------------------------------------------------

process = cms.Process("RecoFastTiming")

# randomness
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   VtxSmeared = cms.PSet(
                                                       initialSeed = cms.untracked.uint32(1),#int(time.time()%100000//1)),
                                                       engineName = cms.untracked.string('TRandom3')
                                                   )
)

randHelper = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.populate()

## load the SK geometry and magnetic field config
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
## import standard RecoFT configurations
process.load("IOMC.EventVertexGenerators.GhostVtxSmearedHLLHC_cfi")
process.load("FastTiming.RecoTreeUtils.RecoFastTiming_cff")

#process.RecoFastTiming.makeGhosts = cms.untracked.bool(True);

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
                            fileNames = filesOpt.inputFiles)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService",
                            fileName = filesOpt.outputFile)

#process.ft_path = cms.Sequence(process.VtxSmeared+process.RecoFastTiming)
process.ft_path = cms.Sequence(process.RecoFastTiming)

process.path = cms.Path(process.ft_path)

process.schedule = cms.Schedule(process.path)#, process.FEVTDEBUGoutput_step)

