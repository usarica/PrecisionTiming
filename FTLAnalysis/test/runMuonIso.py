import subprocess
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('eosdirs',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "files location(s) on EOS")
options.register('localdirs',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "files location(s) on local filesystem")
options.register('datasets',
                 '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Input dataset(s)")
options.register('outname',
                 'muon_iso.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name")
options.register('usegenpv',
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Use MC-truth to define PV")
options.register('isTimingSample',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Process sample using timing and 4D vertexing")
options.register('saveTracksInfo',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Save per-track information")
options.register('runHGCToySim',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Allow to define eta and pt dependent time resolutions")
options.register('targetres',
                 [0.03, 0.05, 0.07, 0.09, 0.15],
                 VarParsing.multiplicity.list,
                 VarParsing.varType.float,
                 "Extra time resolution smearings")
options.register('dzCut',
                 0.1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 "Trk-vtx dz cut")
options.register('debug',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print debug messages")
options.maxEvents = -1
options.parseArguments()

process = cms.Process('TimingAnalysis')

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# import of standard configurations
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

files = []
for dataset in options.datasets:
    print('>> Creating list of files from: \n'+dataset)
    for instance in ['global', 'phys03']:
        query = "-query='file dataset="+dataset+"'"
        if options.debug:
            print(query)
        lsCmd = subprocess.Popen(['dasgoclient '+query+' -limit=0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        str_files, err = lsCmd.communicate()
        files.extend(['root://cms-xrd-global.cern.ch/'+ifile for ifile in str_files.split("\n")])
        files.pop()
        if len(files) > 0:
            break

for eosdir in options.eosdirs:
    if eosdir[-1] != '/':
        eosdir += '/'
    print('>> Creating list of files from: \n'+eosdir)
    lsCmd = subprocess.Popen(['eos', 'ls', eosdir+'*.root'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    str_files, err = lsCmd.communicate()
    files.extend(['root://eoscms/'+eosdir+ifile for ifile in str_files.split("\n")])
    files.pop()

for localdir in options.localdirs:
    if localdir[-1] != '/':
        localdir += '/'
    print('>> Creating list of files from: \n'+localdir)
    lsCmd = subprocess.Popen(['ls '+localdir+'*.root'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    str_files, err = lsCmd.communicate()
    files.extend(['file:'+ifile for ifile in str_files.split("\n")])
    files.pop()

for ifile in options.inputFiles:
    files.append("file:"+ifile)
    
if options.debug:
    for ifile in files:
        print(ifile)
    

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(files),
                            secondaryFileNames=cms.untracked.vstring()
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('PrecisionTiming/FTLAnalysis/python/FTLMuonIsolation_cfi.py nevts:-1'),
    name = cms.untracked.string('Applications')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')


process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(options.outname),
    closeFileFast = cms.untracked.bool(True)
)


if options.runHGCToySim:
    from PrecisionTiming.FTLAnalysis.FTLMuonIsolation_cfi import FTLMuonIsolationHGCToy as FTLMuonIsolation
else:
    from PrecisionTiming.FTLAnalysis.FTLMuonIsolation_cfi import FTLMuonIsolation as FTLMuonIsolation
process.MuonIsolation = FTLMuonIsolation
process.MuonIsolation.useMCTruthPV = options.usegenpv
process.MuonIsolation.isTimingSample = options.isTimingSample
process.MuonIsolation.saveTracksInfo = options.saveTracksInfo
process.MuonIsolation.dzCut = options.dzCut
if process.MuonIsolation.isTimingSample:
   process.MuonIsolation.vtxTag = cms.untracked.InputTag("offlinePrimaryVertices4D", "", "RECO")
##process.MuonIsolation.targetResolutions = cms.vdouble(res for res in options.targetres)
print(process.MuonIsolation.targetResolutions)
print(process.MuonIsolation.useMCTruthPV)
print(process.MuonIsolation.isTimingSample)

process.path = cms.Path(process.MuonIsolation)

# Path and EndPath definitions

# Schedule definition
process.schedule = cms.Schedule(process.path)


