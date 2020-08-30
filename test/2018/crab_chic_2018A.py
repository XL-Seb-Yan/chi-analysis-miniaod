from datetime import date
today = date.today()
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Chic_%s_2018A'%today
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run-chic-miniaod_2018ABC.py'
config.JobType.outputFiles = ['rootuple.root']
config.JobType.sendExternalFolder = True
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/Charmonium/Run2018A-17Sep2018-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.lumiMask = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_MuonPhys.txt'
config.Data.outLFNDirBase = '/store/user/xuyan/'
config.Data.publication = False
config.Data.outputDatasetTag  = 'Chic_%s_2018A'%today
config.Site.storageSite = 'T3_US_FNALLPC'

