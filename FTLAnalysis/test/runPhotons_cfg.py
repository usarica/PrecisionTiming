import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('eosdirs',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "files location(s) on EOS")
options.register('datatier',
                 'RECO',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "")
options.register('hasftl',
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Has FTL geometry and RecHit collections")
options.register('debug',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print debug messages")
options.parseArguments()

process = cms.Process("FTLDumpPhotons")
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# Geometry
process.load('Configuration.Geometry.GeometryExtended2023D8Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D8_cff')

for eosdir in options.eosdirs:
    if eosdir[-1] != '/':
        eosdir += '/'
    print('>> Creating list of files from: \n'+eosdir)
    lsCmd = subprocess.Popen(['eos', 'ls', eosdir+'*.root'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    str_files, err = lsCmd.communicate()
    files = ['root://eoscms/'+eosdir+ifile for ifile in str_files.split("\n")]
    files.pop()
    if options.debug:
        for ifile in files:
            print(ifile)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(files)
)
                            
process.load('PrecisionTiming.FTLAnalysis.FTLDumpPhotons_cfi')
process.load('PrecisionTiming.FTLAnalysis.FTLDumpJets_cfi')
FTLDumperPhotons = process.FTLDumpPhotonsRECO if options.datatier == 'RECO' else process.FTLDumpPhotonsPAT
FTLDumperPhotons.readFTLRecHits = options.hasftl

FTLDumperJets = process.FTLDumpJets

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("ftl_photons.root"))

process.path = cms.Path(FTLDumperPhotons)

process.schedule = cms.Schedule(process.path)
