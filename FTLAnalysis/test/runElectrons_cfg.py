import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('datasets',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Input dataset(s)")
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
options.maxEvents = -1
options.parseArguments()

process = cms.Process("FTLDumpElectrons")
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# import of standard configurations
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Geometry
process.load('Configuration.Geometry.GeometryExtended2023D19Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D19_cff')

files = []
for dataset in options.datasets:
    print('>> Creating list of files from: \n'+dataset)
    query = "-query='file dataset="+dataset+"'"
    if options.debug:
        print(query)
    lsCmd = subprocess.Popen(['dasgoclient '+query+' -limit=0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    str_files, err = lsCmd.communicate()
    files.extend(['root://cms-xrd-global.cern.ch/'+ifile for ifile in str_files.split("\n")])
    files.pop()

for eosdir in options.eosdirs:
    if eosdir[-1] != '/':
        eosdir += '/'
    print('>> Creating list of files from: \n'+eosdir)
    lsCmd = subprocess.Popen(['eos', 'ls', eosdir+'*.root'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    str_files, err = lsCmd.communicate()
    files.extend(['root://eoscms/'+eosdir+ifile for ifile in str_files.split("\n")])
    files.pop()

if options.debug:
    for ifile in files:
        print(ifile)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(files)
)
                            
process.load('PrecisionTiming.FTLAnalysis.FTLDumpElectrons_cfi')
FTLDumper = process.FTLDumpElectronsRECO if options.datatier == 'RECO' else process.FTLDumpElectronsPAT
FTLDumper.readFTLRecHits = options.hasftl

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("ftl_electrons.root"))

process.path = cms.Path(FTLDumper)

process.schedule = cms.Schedule(process.path)
