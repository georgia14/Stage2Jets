import FWCore.ParameterSet.Config as cms

# More robust but less complete version
Stage2JetProducer=cms.EDProducer(
    "Stage2JetProducer",

    #Trigger tower input
    towerInput = cms.InputTag("L1CaloTowerProducer"),
    mhtThreshold = cms.double(10),
    htThreshold = cms.double(30),
    jetEtaCut = cms.double(3.0),


    )

