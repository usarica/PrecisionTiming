from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName   = 'VBF_HToGG_M125_14TeV_PU140BX25_MIB_V2'
config.General.transferLogs = False
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
#config.JobType.pluginName  = 'PrivateMC'
# Name of the CMSSW configuration file
config.JobType.psetName    = 'RecoFastTiming_cfg_forCRAB.py'
#config.JobType.inputFiles = ['gbrv3ele_52x.root', 'gbrv3ph_52x.root']
config.JobType.outputFiles = ['VBFH125GG_140PU.root']

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
# For instance, this dataset will be named /CrabTestSingleMu/something/USER
config.Data.inputDataset = '/VBF_HToGG_M-125_14TeV-powheg-pythia6/TP2023SHCALDR-SHCALTIME_PU140BX25_SHCalTime_PH2_1K_FB_V6-v1/GEN-SIM-RECO'
#config.Data.useParent = True
config.Data.inputDBS = 'global' #'phys03'
config.Data.splitting =  'LumiBased'
config.Data.unitsPerJob = 2
config.Data.totalUnits = 20000
config.Data.publication = False
# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase =  '/store/user/martelli/' # or '/store/group/<subdir>'   #'/store/group/dpg_ecal/alca_ecalcalib/amartell/'   #

config.section_('Site')
# Where the output files will be transmitted to
#config.Site.storageSite = 'T2_CH_CERN' #'srm-eoscms.cern.ch'    #'T2_US_Nowhere'
config.Site.storageSite = 'T3_IT_MIB' #'srm-eoscms.cern.ch'    #'T2_US_Nowhere'
#config.Site.whitelist = ['T2_*', 'T3_*']
