import FWCore.ParameterSet.Config as cms

topGenInfo = cms.EDProducer("TopGenInfoProducer",
    parton = cms.InputTag("partonTop"),
    pseudo = cms.InputTag("pseudoTop"),
)
