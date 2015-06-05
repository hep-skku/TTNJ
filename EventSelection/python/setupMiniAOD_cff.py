import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

def setupMiniAOD(process):
    process.load("RecoMuon.MuonIdentification.Identification.cutBasedMuonId_MuonPOG_V0_cff")
    process.load("RecoMuon.MuonIsolation.muonPFIsolationCitk_cff")
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

    process.muonPFPileUpIsolation.srcToIsolate = "slimmedMuons"
    process.muonPFNoPileUpIsolation.srcToIsolate = "slimmedMuons"
    process.egmGedGsfElectronPFPileUpIsolation.srcToIsolate = "slimmedElectrons"
    process.egmGedGsfElectronPFNoPileUpIsolation.srcToIsolate = "slimmedElectrons"

    process.muonPFPileUpIsolation.srcForIsolationCone = "packedPFCandidates"
    process.muonPFNoPileUpIsolation.srcForIsolationCone = "packedPFCandidates"
    process.egmGedGsfElectronPFPileUpIsolation.srcForIsolationCone = "packedPFCandidates"
    process.egmGedGsfElectronPFNoPileUpIsolation.srcForIsolationCone = "packedPFCandidates"

    if hasattr(process, 'goodOfflinePrimaryVertices'):
        process.goodOfflinePrimaryVertices.src = "offlineSlimmedPrimaryVertices"
    # A trick to switch input vertex collection to slimmed-.
    process.offlinePrimaryVertices = cms.EDFilter("VertexSelector",
        src = cms.InputTag("offlineSlimmedPrimaryVertices"),
        cut = cms.string("")
    )


