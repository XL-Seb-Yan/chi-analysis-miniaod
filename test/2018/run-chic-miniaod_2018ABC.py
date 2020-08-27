#input_filename = '/store/data/Run2018D/Charmonium/MINIAOD/PromptReco-v2/000/322/625/00000/68679BCD-B4B9-E811-B6EE-FA163E1B57DB.root'
input_filename = '/store/data/Run2018A/Charmonium/MINIAOD/PromptReco-v2/000/316/240/00000/9E9A8843-9459-E811-93DC-FA163EEAACDE.root'
ouput_filename = 'rootuple.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v13', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

import FWCore.PythonUtilities.LumiList as LumiList
lumi_list = LumiList.LumiList(filename ='Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt').getVLuminosityBlockRange()
process.source.lumisToProcess = lumi_list

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_Dimuon0_Jpsi3p5_Muon2_v*',
                                                                        'HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v*',
                                                                        'HLT_Dimuon0_Jpsi_L1_NoOS_v*',
                                                                        'HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v*',
                                                                        'HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v*',
                                                                        'HLT_Dimuon0_Jpsi_NoVertexing_v*',
                                                                        'HLT_Dimuon0_Jpsi_v*',
                                                                        'HLT_Dimuon0_LowMass_L1_0er1p5R_v*',
                                                                        'HLT_Dimuon0_LowMass_L1_0er1p5_v*',
                                                                        'HLT_Dimuon0_LowMass_L1_4R_v*',
                                                                        'HLT_Dimuon0_LowMass_L1_4_v*',
                                                                        'HLT_Dimuon0_LowMass_v*',
                                                                        'HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
                                                                        'HLT_Dimuon25_Jpsi_noCorrL1_v*',
                                                                        'HLT_Dimuon25_Jpsi_v*',
                                                                        'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v*',
                                                                        'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v*',
                                                                        'HLT_DoubleMu4_3_Jpsi_v*',
                                                                        'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
                                                                        'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
                                                                        'HLT_DoubleMu4_Jpsi_Displaced_v*',
                                                                        'HLT_DoubleMu4_Jpsi_NoVertexing_v*',
                                                                        'HLT_Mu7p5_L2Mu2_Jpsi_v*',
                                                                        'HLT_Mu7p5_Track2_Jpsi_v*',
                                                                        'HLT_Mu7p5_Track3p5_Jpsi_v*',
                                                                        'HLT_Mu7p5_Track7_Jpsi_v*'
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(True)
                                        )

process.load("Ponia.OniaPhoton.slimmedMuonsTriggerMatcher_cfi")

# In MiniAOD, the PATMuons are already present. We just need to run Onia2MuMu, with a selection of muons.
process.oniaSelectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuonsWithTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    ' && (abs(eta) <= 1.4 && pt > 4.)'
   ),
   filter = cms.bool(True)
)

process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('oniaSelectedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("")
process.onia2MuMuPAT.dimuonSelection=cms.string("2.7 < mass && mass < 3.5")
process.onia2MuMuPAT.addMCTruth = cms.bool(False)

process.Onia2MuMuFiltered = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("onia2MuMuPAT"),
      singlemuonSelection = cms.string(""),             
      dimuonSelection     = cms.string("2.9 < mass && mass < 3.3 && pt > 20. && abs(y) < 1.2 && charge==0 && userFloat('vProb') > 0.01"),
      do_trigger_match    = cms.bool(True),
      HLTPaths          = cms.vstring('HLT_Dimuon0_Jpsi3p5_Muon2_v*',
                                      'HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v*',
                                      'HLT_Dimuon0_Jpsi_L1_NoOS_v*',
                                      'HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v*',
                                      'HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v*',
                                      'HLT_Dimuon0_Jpsi_NoVertexing_v*',
                                      'HLT_Dimuon0_Jpsi_v*',
                                      'HLT_Dimuon0_LowMass_L1_0er1p5R_v*',
                                      'HLT_Dimuon0_LowMass_L1_0er1p5_v*',
                                      'HLT_Dimuon0_LowMass_L1_4R_v*',
                                      'HLT_Dimuon0_LowMass_L1_4_v*',
                                      'HLT_Dimuon0_LowMass_v*',
                                      'HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*',
                                      'HLT_Dimuon25_Jpsi_noCorrL1_v*',
                                      'HLT_Dimuon25_Jpsi_v*',
                                      'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v*',
                                      'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v*',
                                      'HLT_DoubleMu4_3_Jpsi_v*',
                                      'HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*',
                                      'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
                                      'HLT_DoubleMu4_Jpsi_Displaced_v*',
                                      'HLT_DoubleMu4_Jpsi_NoVertexing_v*',
                                      'HLT_Mu7p5_L2Mu2_Jpsi_v*',
                                      'HLT_Mu7p5_Track2_Jpsi_v*',
                                      'HLT_Mu7p5_Track3p5_Jpsi_v*',
                                      'HLT_Mu7p5_Track7_Jpsi_v*'
                                      ),
)

process.DiMuonCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFiltered"),
    minNumber = cms.uint32(1),
)

process.chiProducer = cms.EDProducer('OniaPhotonProducer',
    conversions     = cms.InputTag("oniaPhotonCandidates","conversions"),
    dimuons         = cms.InputTag("Onia2MuMuFiltered"),
    pi0OnlineSwitch = cms.bool(False),
    deltaMass       = cms.vdouble(0.0,2.0),
    dzmax           = cms.double(0.5),
    triggerMatch    = cms.bool(False)  # trigger match is performed in Onia2MuMuFiltered
)

process.chiFitter = cms.EDProducer('OniaPhotonKinematicFit',
                          chi_cand = cms.InputTag("chiProducer"),
                          meson_nS_mass = cms.double(3.0969), # GeV   Y1S = 9.46030   Y2S = 10.02326    Y3S = 10.35520  J/psi(1S)=3.0969
                          product_name = cms.string("jpsi"), # means we are looking for chi_c decaying to J/psi 
                          is_Debug = cms.bool(True)
                         )
                         
process.XFitter = cms.EDProducer('XDecayTreeKinematicFit',
                          chi_cand = cms.InputTag("chiProducer"),
                          meson_nS_cand = cms.InputTag("Onia2MuMuFiltered"),
                          track = cms.InputTag("packedPFCandidates"),
                          meson_nS_mass = cms.double(3.0969), # GeV  J/psi(1S)=3.0969
                          product_name = cms.string("XCharged"), # means we are looking for charged partner of X3872
                          is_Debug = cms.bool(True)
                         )

process.XSequence = cms.Sequence(
                                   process.triggerSelection *
                                   process.slimmedMuonsWithTriggerSequence *
                                   process.oniaSelectedMuons *
                                   process.onia2MuMuPAT *
                                   process.Onia2MuMuFiltered *
                                   process.DiMuonCounter *
                                   process.chiProducer *
                                   process.chiFitter *
                                   process.XFitter
				   )

process.rootuple = cms.EDAnalyzer('chicRootupler',
                          chi_cand = cms.InputTag("chiProducer"),
                          meson_nS_cand = cms.InputTag("Onia2MuMuFiltered"),
                          refit1P  = cms.InputTag("chiFitter","jpsi"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
                          isMC = cms.bool(False),
                          GenParticles = cms.InputTag("prunedGenParticles")
                         )

process.p = cms.Path(process.XSequence*process.rootuple)