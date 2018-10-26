from WMCore.Configuration import Configuration

#sample  = "QCD"
#sample = "TTbar"
sample = "DY"

tag = "nTracks_dz0p1_pt0p9_v3"

if sample == "QCD":
    in_dataset = '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
elif sample == "TTbar":
    # in_dataset = '/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIITDRFall17DR-PU200_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
    in_dataset = '/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/PhaseIISpr18AODMiniAOD-PU200_93X_upgrade2023_realistic_v5_ext1-v2/AODSIM'
elif sample == "DY":
    in_dataset = '/DYToLL-M-50_0J_14TeV-madgraphMLM-pythia8/PhaseIISpr18AODMiniAOD-PU200_93X_upgrade2023_realistic_v5-v1/AODSIM'

config = Configuration()

config.section_('General')
config.General.requestName       = sample+'_MuonIsolationMTD_200PU_932_HGCparam_'+tag
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName          = '../test/runMuonIso.py'
config.JobType.pyCfgParams       = ['usegenpv=0']
config.JobType.priority          = 30

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
## 200 PU
config.Data.inputDataset         = in_dataset
config.Data.inputDBS             = 'global'
config.Data.splitting            = 'FileBased'
config.Data.unitsPerJob          = 1
config.Data.totalUnits           = -1
config.Data.publication          = False
#config.Data.ignoreLocality       = True
config.Data.allowNonValidInputDataset = True

# This string is used to construct the output dataset name
# config.Data.publishWithGroupName = True
config.Data.outLFNDirBase        = '/store/user/usarica/MTD/UPS/MuonIsolation/20181024/'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite          = 'T2_US_UCSD'
#config.Site.whitelist            = ['T1_IT_CNAF', 'T2_CH_CERN', 'T2_US_MIT', 'T2_US_Florida']
config.Site.blacklist            = ['T2_US_Nebraska','T2_DE_DESY']
