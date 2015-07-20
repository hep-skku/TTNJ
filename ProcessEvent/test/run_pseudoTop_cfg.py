import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("Configuration.StandardSequences.Geometry_cff")

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_7_4_4/RelValTTbar_13/MINIAODSIM/PU50ns_MCRUN2_74_V8_38Tbis-v1/00000/42193EE1-3F09-E511-94E8-0025905B8596.root',
        '/store/relval/CMSSW_7_4_4/RelValTTbar_13/MINIAODSIM/PU50ns_MCRUN2_74_V8_38Tbis-v1/00000/C47C531D-4009-E511-A0A0-0025905B85EE.root',
    ),
)

process.load("TTNJ.ProcessEvent.partonTop_cfi")
process.load("TopQuarkAnalysis.TopEventProducers.producers.pseudoTop_cfi")
process.load("TTNJ.ProcessEvent.topGenInfo_cfi")

process.ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("error"), # choose one among keep/skip/error
    eventCounters = cms.vstring(), #"nEventsTotal", "nEventsClean", "nEventsPAT"),
    int = cms.PSet(
        parton_ch = cms.PSet(src = cms.InputTag("topGenInfo", "partonChannel")),
        pseudo_ch = cms.PSet(src = cms.InputTag("topGenInfo", "pseudoChannel")),
        parton_ac = cms.PSet(src = cms.InputTag("topGenInfo", "partonAccept")),
    ),
    ints = cms.PSet(
        parton_modes = cms.PSet(src = cms.InputTag("partonTop", "modes")),
    ),
    cands = cms.PSet(
        parton = cms.PSet(
            src = cms.InputTag("partonTop"),
            exprs = cms.untracked.PSet(
                pt = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m = cms.string("mass"),
                pdgId = cms.string("pdgId"),
            ),
        ),
        pseudo = cms.PSet(
            src = cms.InputTag("pseudoTop"),
            exprs = cms.untracked.PSet(
                pt = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m = cms.string("mass"),
                pdgId = cms.string("pdgId"),
            ),
        )
    ),
)

process.p = cms.Path(
    process.partonTop
#    process.pseudoTop
#  * process.topGenInfo
#  * process.ntuple
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('f.root')
)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("o.root"),
#    outputCommands = cms.untracked.vstring(
#        "drop *",
#        "keep *_*_*_ANA",
#    ),
#)
#process.outPath = cms.EndPath(process.out)
