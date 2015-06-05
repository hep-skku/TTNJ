import FWCore.ParameterSet.Config as cms

topDileptonObjects = cms.EDFilter("TopDileptonProducer",
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    jet = cms.PSet(
        src = cms.InputTag("slimmedJets"),
        minPt = cms.untracked.double(30),
        maxEta = cms.untracked.double(2.4),
    ),
    met = cms.InputTag("slimmedMETs"),
    muon = cms.PSet(
        src = cms.InputTag("slimmedMuons"),
        minPt = cms.untracked.double(20),
        maxEta = cms.untracked.double(2.4),
        maxMaxRiso = cms.untracked.double(0.2),
        id = cms.InputTag("muoMuonIDs", "cutBasedMuonId-MuonPOG-V0-loose"),
        chIso = cms.InputTag("muonPFNoPileUpIsolation", "h+-DR030-ThresholdVeto000-ConeVeto000"),
        nhIso = cms.InputTag("muonPFNoPileUpIsolation", "h0-DR030-ThresholdVeto050-ConeVeto000"),
        phIso = cms.InputTag("muonPFNoPileUpIsolation", "gamma-DR030-ThresholdVeto050-ConeVeto000"),
        puIso = cms.InputTag("muonPFPileUpIsolation", "h+-DR040-ThresholdVeto050-ConeVeto000"),
    ),
    electron = cms.PSet(
        src = cms.InputTag("slimmedElectrons"),
        minPt = cms.untracked.double(20),
        maxEta = cms.untracked.double(2.4),
        maxMaxRiso = cms.untracked.double(0.2),
        id = cms.InputTag("egmGsfElectronIDs", "cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
        chIso = cms.InputTag("egmGedGsfElectronPFNoPileUpIsolation", "h+-DR030-BarVeto000-EndVeto001"),
        nhIso = cms.InputTag("egmGedGsfElectronPFNoPileUpIsolation", "h0-DR030-BarVeto000-EndVeto000"),
        phIso = cms.InputTag("egmGedGsfElectronPFNoPileUpIsolation", "gamma-DR030-BarVeto000-EndVeto008"),
        puIso = cms.InputTag("egmGedGsfElectronPFPileUpIsolation", "h+-DR030-BarVeto000-EndVeto001"),
    ),
    cutStepToAccept = cms.int32(3),
)

