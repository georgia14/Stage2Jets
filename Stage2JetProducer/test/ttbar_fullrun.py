import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'WARNING'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(1)
#        )
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# import of standard configurations
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
                                                                    
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use

                            fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_2_1/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_72_V1-v1/00000/0AC448E5-295D-E411-B7D5-00261894389C.root')

#    fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_3_0_pre3/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_73_V5-v1/00000/4E090269-2776-E411-8205-02163E00E909.root')

#/store/relval/CMSSW_7_3_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU50ns_MCRUN2_73_V6-v1/00000/0C7DABD6-4181-E411-B86E-0025905964BA.root')

#                            fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_3_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/0496C828-4981-E411-8142-0025905A60B8.root')

#/store/relval/CMSSW_7_3_0_pre1/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_LS172_V15-v1/00000/08DC28C8-E059-E411-B9C6-0025905B85AA.root')

#fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_2_0_pre5/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/POSTLS172_V3-v1/00000/323F8170-4A30-E411-82CB-0025905A610A.root')


#    fileNames = cms.untracked.vstring(

 #       '/store/relval/CMSSW_7_1_0/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/POSTLS171_V15-v1/00000/24ADC675-6BFB-E311-9444-0025905A611E.root'
#'/store/relval/CMSSW_7_2_0_pre5/RelValTTbar_13/GEN-SIM-DIGI-RAW-HLTDEBUG/POSTLS172_V3-v1/00000/323F8170-4A30-E411-82CB-0025905A610A.root'
  #  )
#  skipEvents=cms.untracked.uint32(54500)
)

process.GlobalTag.globaltag = 'MCRUN2_71_V0::All'

# Raw to digi
process.load('Configuration.StandardSequences.RawToDigi_cff')

#Load stuff for making the stage2 jets
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")
process.load("Stage2Jets.Stage2JetProducer.Stage2JetProducer_cfi")

# upgrade calo stage 2
process.load('L1Trigger.L1TCalorimeter.L1TCaloStage2_PPFromRaw_cff')
#process.load('L1Trigger.L1TCalorimeter.l1tStage2CaloAnalyzer_cfi')


process.makeTestTree = cms.EDAnalyzer('MakeTestTree',
    )

process.p = cms.Path(
        process.RawToDigi +
        process.ecalDigis +
        process.hcalDigis +
        process.L1TCaloStage2_PPFromRaw +

        #Stage 2 test jets step
        process.SLHCCaloTrigger +
        process.Stage2JetProducer +
        process.makeTestTree
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('ttbar_test.root')
                                   )

# Output definition
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('ttbar_test_skim.root'),
    outputCommands = cms.untracked.vstring(
      'keep *',
      )
    )

#process.out = cms.EndPath(process.output)
