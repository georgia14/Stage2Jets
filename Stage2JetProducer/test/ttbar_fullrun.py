import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'WARNING'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(1)
#        )
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
                                                                    
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:/afs/cern.ch/work/j/jad/public/TriggerTowers_TTbar_13TeV_PU40_25ns_12345_skimmed.root',
    #'file:/afs/cern.ch/work/j/jad/public/TriggerTowers_TTbar_13TeV_PU40_25ns_678910_skimmed.root'
    #'file:/afs/cern.ch/work/j/jad/public/TriggerTowers_TTbar_13TeV_PU40_25ns_1-10_skimmed_v2.root'
    #'file:/afs/cern.ch/work/j/jad/public/TriggerTowers_TTbar_13TeV_noPU_test2_v3.root'
    #'file:/afs/cern.ch/work/j/jad/public/TriggerTowers_TTbar_13TeV_PU40_50ns_1-10_skimmed_v3.root'    
    #'file:/afs/cern.ch/work/a/aelwood/trigger/jets/jadjet/CMSSW_6_2_0/src/TtBar_BigSkim_12May14.root'
 #'file:/afs/cern.ch/work/a/aelwood/ntuples/trigger/jets/jadjet/TriggerTowers_TTbar_13TeV_PU40_50ns_1-10_skimmed_v3.root'

# '/store/user/aelwood/140725_jJets_ngun_PU40bx50/L1SkimRun_40_1_uha.root'

'root://eoscms.cern.ch//eos/cms/store/user/aelwood/140725_jJets_ttbar_PU40bx50/L1SkimRun_48_1_Wg8.root'

#'file://l1_skim_adam.root'

#'/store/user/aelwood/140801_jJets_ttbar_PU40bx25/L1SkimRun_48_1_pDL.root'

    #'file:/afs/cern.ch/work/m/mcitron/public/TriggerTowers_NeutrinoGun_13TeV_PU40_50ns_1-10_skimmed_v3.root'
    ),
    #skipEvents=cms.untracked.uint32(54500)
)

process.load("Stage2Jets.Stage2JetProducer.Stage2JetProducer_cfi")

process.makeTestTree = cms.EDAnalyzer('MakeTestTree',
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('ttbar_test.root')
                                   #fileName = cms.string('histograms-ttbar_13TeV_PU40_50ns_test_v3.root')
                                   )

process.p = cms.Path(
    process.Stage2JetProducer
    +process.makeTestTree
    )

# Output definition
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('ttbar_test_skim.root'),
    outputCommands = cms.untracked.vstring(
      'keep *',
      )
    )

process.out = cms.EndPath(process.output)
