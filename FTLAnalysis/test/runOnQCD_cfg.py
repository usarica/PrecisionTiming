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

process = cms.Process("FTLDumpHgg")
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
    
if len(options.inputFiles) > 0:
    files = options.inputFiles

if options.debug:
    for ifile in files:
        print(ifile)

files = [
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_1.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_10.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_11.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_12.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_13.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_14.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_15.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_16.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_17.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_18.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_19.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_2.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_20.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_21.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_22.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_23.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_24.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_25.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_26.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_27.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_28.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_29.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_3.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_30.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_31.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_32.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_33.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_34.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_35.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_36.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_37.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_38.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_39.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_4.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_40.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_41.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_42.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_43.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_44.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_45.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_46.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_47.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_48.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_49.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_5.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_50.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_51.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_52.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_53.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_54.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_55.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_56.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_57.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_58.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_59.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_6.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_60.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_61.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_62.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_63.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_64.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_65.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_66.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_67.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_68.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_69.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_7.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_70.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_71.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_72.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_73.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_74.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_75.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_76.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_77.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_78.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_79.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_8.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_80.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_81.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_82.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_83.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_84.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_85.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_86.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_87.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_88.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_89.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_9.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_90.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_91.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_92.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_93.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_94.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_95.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_96.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_97.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_98.root',
    '/store/group/upgrade/timing/QCD_Pt-100To1000_14TeV/crab_QCD_Pt-100To1000_TuneCUETP8M1_Flat_14TeV-pythia8_RECO_D8FTL/170314_153118/0000/step3_RAW2DIGI_L1Reco_RECO_PAT_99.root'
]
        
# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(files)
)
                            
process.load('PrecisionTiming.FTLAnalysis.FTLDumpJets_cfi')
FTLDumperJets = process.FTLDumpJets

# Output TFile
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("ftl_qcd.root"))

process.path = cms.Path(FTLDumperJets)

process.schedule = cms.Schedule(process.path)
