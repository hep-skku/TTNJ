import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("Configuration.StandardSequences.Geometry_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.GlobalTag.globaltag = ''
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_7_4_0_pre8/RelValProdTTbar_13/MINIAODSIM/MCRUN2_74_V7-v1/00000/7CBFDA90-55BD-E411-B15C-0025905964BA.root',
  ),
)

process.load("RecoMuon.MuonIdentification.Identification.cutBasedMuonId_MuonPOG_V0_cff")
process.load("RecoMuon.MuonIsolation.muonPFIsolationCitk_cff")
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDMuonIdProducer(process, DataFormat.MiniAOD)
setupVIDMuonSelection(process, process.cutBasedMuonId_MuonPOG_V0_loose)
setupVIDMuonSelection(process, process.cutBasedMuonId_MuonPOG_V0_medium)
setupVIDMuonSelection(process, process.cutBasedMuonId_MuonPOG_V0_tight)
setupVIDMuonSelection(process, process.cutBasedMuonId_MuonPOG_V0_soft)
setupVIDMuonSelection(process, process.cutBasedMuonId_MuonPOG_V0_highpt)

process.load("RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff")
process.load("RecoEgamma.EgammaIsolationAlgos.egmGedGsfElectronPFIsolation_cfi")
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
setupVIDElectronSelection(process, process.cutBasedElectronID_PHYS14_PU20bx25_V2_standalone_loose)
setupVIDElectronSelection(process, process.cutBasedElectronID_PHYS14_PU20bx25_V2_standalone_medium)
setupVIDElectronSelection(process, process.cutBasedElectronID_PHYS14_PU20bx25_V2_standalone_tight)

process.load("TTNJ.EventSelection.topCommonObjects_cff")
process.load("TTNJ.EventSelection.topDileptonProducer_cfi")
process.goodOfflinePrimaryVertices.src = "offlineSlimmedPrimaryVertices"
process.offlinePrimaryVertices = cms.EDFilter("VertexSelector", src = cms.InputTag("offlineSlimmedPrimaryVertices"), cut = cms.string(""))

process.muonPFPileUpIsolation.srcToIsolate = "slimmedMuons"
process.muonPFNoPileUpIsolation.srcToIsolate = "slimmedMuons"
process.egmGedGsfElectronPFPileUpIsolation.srcToIsolate = "slimmedElectrons"
process.egmGedGsfElectronPFNoPileUpIsolation.srcToIsolate = "slimmedElectrons"

process.muonPFPileUpIsolation.srcForIsolationCone = "packedPFCandidates"
process.muonPFNoPileUpIsolation.srcForIsolationCone = "packedPFCandidates"
process.egmGedGsfElectronPFPileUpIsolation.srcForIsolationCone = "packedPFCandidates"
process.egmGedGsfElectronPFNoPileUpIsolation.srcForIsolationCone = "packedPFCandidates"

process.p = cms.Path(
    process.topDileptonObjects
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('f.root')
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("o.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_top*_*_ANA",
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
)
process.outPath = cms.EndPath(process.out)

