import FWCore.ParameterSet.Config as cms

ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("keep"), # choose one among keep/skip/error
    #failureMode = cms.untracked.string("error"), # choose one among keep/skip/error
    eventCounters = cms.vstring(), #"nEventsTotal", "nEventsClean", "nEventsPAT"),
    int = cms.PSet(
        mode = cms.PSet(src = cms.InputTag("topDileptonObjects", "mode")),
        nVertex = cms.PSet(src = cms.InputTag("topDileptonObjects", "nVertex")),
        step = cms.PSet(src = cms.InputTag("topDileptonObjects", "step")),
        mc_channel = cms.PSet(src = cms.InputTag("partonTop", "channel")),
    ),
    vint = cms.PSet(
        mc_modes = cms.PSet(src = cms.InputTag("partonTop", "modes")),
    ),
    double = cms.PSet(
        mLL = cms.PSet(src = cms.InputTag("topDileptonObjects", "mLL")),
    ),
    vdouble = cms.PSet(
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
        ),
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
    ),
)
