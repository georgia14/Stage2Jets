# Stage2Jets
Private producer used to prototype Stage2 jet and energy sum algorithms [originally developed by Jad, Adam]

# Checkout the standalone code, for porting to CMSSW_73X:
> cd CMSSW_7_3_0_pre1/src/
> cmsenv
> git clone https://github.com/georgia14/Stage2Jets.git

# To find RelVal samples for a particular CMSSW release (e.g. RelVal TTbar):
> eos ls /eos/cms/store/relval/CMSSW_7_1_0/RelValTTbar_13/

# To find the global tag to use, go to:
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_2012_data_taking

# Added following line in cfg:
process.GlobalTag.globaltag = 'MCRUN2_71_V0::All'

# To run the standalone code for stage2 Jets:
> cd CMSSW_7_3_0_pre1/src/Stage2Jets/Stage2JetProducer/test/
> cmsRun ttbar_fullrun.py
