import FWCore.ParameterSet.Config as cms


process = cms.Process('TimingAnalysis')

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryExtended2023D20Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D20_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100 ) )

# # Input source
process.source = cms.Source("PoolSource",
#    secondaryFileNames = cms.untracked.vstring(),
#     # fileNames = cms.untracked.vstring(files)
    fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch///store/mc/PhaseIISpr18AODMiniAOD/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/AODSIM/PU200_93X_upgrade2023_realistic_v5-v1/90000/FEA5A761-EE47-E811-9C32-0CC47A4C8E38.root")
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(15))

process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.30 $'),
    annotation = cms.untracked.string('PrecisionTiming/FTLAnalysis/python/FTLMuonIsolation_cfi.py nevts:-1'),
    name = cms.untracked.string('Applications')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')


process.TFileService = cms.Service("TFileService",
    fileName = cms.string("testfile_DY.root"),
    closeFileFast = cms.untracked.bool(True)
)


from PrecisionTiming.FTLAnalysis.FTLMuonIsolation_cfi import FTLMuonIsolation
process.MuonIsolation = FTLMuonIsolation
process.MuonIsolation.recordTrackInfo = cms.untracked.bool(True)
process.MuonIsolation.recordVertexInfo = cms.untracked.bool(True)

process.path = cms.Path(process.MuonIsolation)
process.outpath = cms.EndPath(process.MuonIsolation)

