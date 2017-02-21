import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("FTLDumpZee")
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Global tag
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v10')

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'root://eoscms//eos/cms/store/relval/CMSSW_9_0_0_pre4/RelValWE_14TeV/MINIAODSIM/90X_upgrade2023_realistic_v3_2023D4Timing-v1/10000/581A0CA3-A9EC-E611-8EB6-0CC47A4D7628.root'
                            )
)                                

process.load('PrecisionTiming.FTLDumper.FTLDumpElectrons_cfi')

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("ftl_zee.root"))

process.path = cms.Path(process.FTLDumpElectrons)

process.schedule = cms.Schedule(process.path)
