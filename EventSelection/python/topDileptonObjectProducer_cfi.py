import FWCore.ParameterSet.Config as cms

topDileptonObjects = cms.EDFilter("TopDileptonObjectProducer",
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    jet = cms.InputTag("goodJets"),
    met = cms.InputTag("patMETs"),
    muon = cms.PSet(
        src = cms.InputTag("patMuons"),
        id = cms.InputTag("muoMuonIDs", "cutBasedMuonId-MuonPOG-V0-loose"),
        chIso = cms.InputTag("muonPFNoPileupIsolation", "h+-DR030-ThresholdVeto000-ConeVeto000"),
        nhIso = cms.InputTag("muonPFNoPileupIsolation", "h0-DR030-ThresholdVeto050-ConeVeto000"),
        phIso = cms.InputTag("muonPFNoPileupIsolation", "gamma-DR030-ThresholdVeto050-ConeVeto000"),
        puIso = cms.InputTag("muonPFPileUpIsolation", "h+-DR040-ThresholdVeto050-ConeVeto000"),
    ),
    electron = cms.PSet(
        src = cms.InputTag("patElectrons"),
        id = cms.InputTag("eleElectronIDs", "cutBasedElectronId-ElectronPOG-V0-loose"),
        chIso = cms.InputTag("electronPFNoPileupIsolation", "h+-DR030-ThresholdVeto000-ConeVeto000"),
        nhIso = cms.InputTag("electronPFNoPileupIsolation", "h0-DR030-ThresholdVeto050-ConeVeto000"),
        phIso = cms.InputTag("electronPFNoPileupIsolation", "gamma-DR030-ThresholdVeto050-ConeVeto000"),
        puIso = cms.InputTag("electronPFPileUpIsolation", "h+-DR040-ThresholdVeto050-ConeVeto000"),
    ),
    cutStepToAccept = cms.int32(5), 
)

