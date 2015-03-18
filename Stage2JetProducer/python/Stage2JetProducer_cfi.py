import FWCore.ParameterSet.Config as cms

# More robust but less complete version
Stage2JetProducer=cms.EDProducer(
    "Stage2JetProducer",

    #Trigger tower input
    towerToken = cms.InputTag("caloStage2Digis", "MP"), # Added to make compatible with emulator
    towerInput = cms.InputTag("L1CaloTowerProducer"),
    mhtThreshold = cms.double(30),
    htThreshold = cms.double(30),
    jetEtaCut = cms.int32(28)


    )

