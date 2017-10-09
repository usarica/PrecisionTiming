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
options.register('runOnMiniAOD',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Sample is in MINIAOD dataformat")
options.register('targetres',
                 [0.03, 0.05, 0.07, 0.09],
                 VarParsing.multiplicity.list,
                 VarParsing.varType.float,
                 "Extra time resolution smearings")
options.register('debug',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Print debug messages")
options.parseArguments()

process = cms.Process('TimingAnalysis')

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

files = []
for dataset in options.datasets:
    print('>> Creating list of files from: \n'+dataset)
    for instance in ['global', 'phys03']:
        query = "--query='file instance=prod/"+instance+" dataset="+dataset+"'"
        if options.debug:
            print(query)
        lsCmd = subprocess.Popen(['das_client.py '+query+' --limit=0'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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
                            secondaryFileNames=cms.untracked.vstring(
                                # "/store/mc/PhaseIITDRSpring17DR/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v6/110000/80071D5E-4D85-E711-8EA0-6C3BE5B5F218.root",
                                # "/store/mc/PhaseIITDRSpring17DR/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v6/110001/28F8AAA6-8885-E711-A734-001CC47D589C.root",
                                # "/store/mc/PhaseIITDRSpring17DR/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v6/110001/925C97EF-8485-E711-A4AD-B499BAAC039C.root",
                                # "/store/mc/PhaseIITDRSpring17DR/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v6/110001/AE5310F8-7785-E711-B0A9-001CC4A63C2A.root",
                                # "/store/mc/PhaseIITDRSpring17DR/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v6/110001/E05CB668-7B85-E711-A5F0-B499BAAC09BE.root",
                                # "/store/mc/PhaseIITDRSpring17DR/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v6/110001/E29563A2-8B85-E711-A929-B499BAAC04F0.root",
                                # "/store/mc/PhaseIITDRSpring17DR/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v6/110002/8CC629A3-7185-E711-94A1-001CC4A63C8E.root",
                                # "/store/mc/PhaseIITDRSpring17DR/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_14TeV_pythia8/GEN-SIM-RECO/PU200_91X_upgrade2023_realistic_v3-v6/110004/B404C769-7585-E711-97AD-B499BAAC0676.root"
                            )
)

process.options = cms.untracked.PSet(

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


if options.runOnMiniAOD:
    from PrecisionTiming.FTLAnalysis.FTLMuonIsolation_cfi import FTLMuonIsolationMiniAOD as FTLMuonIsolation
else:
    from PrecisionTiming.FTLAnalysis.FTLMuonIsolation_cfi import FTLMuonIsolation as FTLMuonIsolation
process.MuonIsolation = FTLMuonIsolation
process.MuonIsolation.useMCTruthPV = options.usegenpv
process.MuonIsolation.isTimingSample = options.isTimingSample
#process.MuonIsolation.targetResolutions = cms.vdouble(res for res in options.targetres)
print(process.MuonIsolation.targetResolutions)
print(process.MuonIsolation.useMCTruthPV)
print(process.MuonIsolation.isTimingSample)

process.path = cms.Path(process.MuonIsolation)

# Path and EndPath definitions

# Schedule definition
process.schedule = cms.Schedule(process.path)

