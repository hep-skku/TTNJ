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

process.load("TTNJ.ProcessEvent.partonTop_cfi")
process.load("TTNJ.ProcessEvent.topDileptonProducer_cfi")
from TTNJ.ProcessEvent.setupMiniAOD_cff import setupMiniAOD
setupMiniAOD(process)
process.load("TTNJ.ProcessEvent.ntuple_cff")

process.p = cms.Path(
    process.partonTop
  + process.topDileptonObjects
  * process.ntuple
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('f.root')
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("o.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_ANA",
    ),
)
#process.outPath = cms.EndPath(process.out)
