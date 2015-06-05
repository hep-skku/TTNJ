import FWCore.ParameterSet.Config as cms

ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("error"), # choose one among keep/skip/error
    eventCounters = cms.vstring(), #"nEventsTotal", "nEventsClean", "nEventsPAT"),
    int = cms.PSet(
        mode = cms.PSet(src = cms.InputTag("topDileptonObjects", "mode")),
        nVertex = cms.PSet(src = cms.InputTag("topDileptonObjects", "nVertex")),
        step = cms.PSet(src = cms.InputTag("topDileptonObjects", "step")),
    ),
    double = cms.PSet(
        mLL = cms.PSet(src = cms.InputTag("topDileptonObjects", "mLL")),
    ),
    doubles = cms.PSet(
    ),
    cands = cms.PSet(
        leptons = cms.PSet(
            src = cms.InputTag("topDileptonObjects", "leptons"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                q = cms.string("charge"),
            ),
            selections = cms.untracked.PSet(),
        ),
        jets = cms.PSet(
            src = cms.InputTag("topDileptonObjects", "jets"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m = cms.string("mass"),
                cisv2 = cms.string("bDiscriminator('combinedInclusiveSecondaryVertexV2BJetTags')"),
            ),
            selections = cms.untracked.PSet(),
        ),
    ),
)
