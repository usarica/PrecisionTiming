import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("FTLDumpHgg")
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Global tag
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v10')

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:step3_RAW2DIGI_L1Reco_RECO_PAT.root'
                            )
)                                

process.load('PrecisionTiming.FTLDumper.FTLDumpPhotons_cfi')

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("ftl_hgg.root"))

process.path = cms.Path(process.FTLDumpPhotons)

process.schedule = cms.Schedule(process.path)
