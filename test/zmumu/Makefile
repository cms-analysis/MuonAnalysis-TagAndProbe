create: tp_from_skim_Zmumu tp_from_skim_Wmunu tp_from_skim_TTbar
status: tp_from_skim_Zmumu/res/CMSSW_1.stdout tp_from_skim_Wmunu/res/CMSSW_1.stdout  tp_from_skim_TTbar/res/CMSSW_1.stdout
getoutput: tp_from_skim_Zmumu/res/CMSSW_1.stderr tp_from_skim_Wmunu/res/CMSSW_1.stderr  tp_from_skim_TTbar/res/CMSSW_1.stderr
fit: fit_result_GlbFromCal.root


fit_result_GlbFromCal.root: histo_dataonly_GlbFromCal.root histo_dataonly_GlbFromTk.root  histo_dataonly_HltFromGlb.root histo_dataonly_TkFromSta.root
	cmsRun fit.py

histo_dataonly_GlbFromTk.root: tp_from_skim_Zmumu/res/histo_dataonly_GlbFromTk_1.root tp_from_skim_Wmunu/res/histo_dataonly_GlbFromTk_1.root tp_from_skim_TTbar/res/histo_dataonly_GlbFromTk_1.root
	hadd -f histo_dataonly_GlbFromTk.root tp_from_skim_Zmumu/res/histo_dataonly_GlbFromTk_*.root tp_from_skim_Wmunu/res/histo_dataonly_GlbFromTk_*.root tp_from_skim_TTbar/res/histo_dataonly_GlbFromTk_*.root

histo_dataonly_GlbFromCal.root: tp_from_skim_Zmumu/res/histo_dataonly_GlbFromCal_1.root tp_from_skim_Wmunu/res/histo_dataonly_GlbFromCal_1.root tp_from_skim_TTbar/res/histo_dataonly_GlbFromCal_1.root
	hadd -f histo_dataonly_GlbFromCal.root tp_from_skim_Zmumu/res/histo_dataonly_GlbFromCal_*.root tp_from_skim_Wmunu/res/histo_dataonly_GlbFromCal_*.root tp_from_skim_TTbar/res/histo_dataonly_GlbFromCal_*.root

histo_dataonly_TkFromSta.root: tp_from_skim_Zmumu/res/histo_dataonly_TkFromSta_1.root tp_from_skim_Wmunu/res/histo_dataonly_TkFromSta_1.root tp_from_skim_TTbar/res/histo_dataonly_TkFromSta_1.root
	hadd -f histo_dataonly_TkFromSta.root tp_from_skim_Zmumu/res/histo_dataonly_TkFromSta_*.root tp_from_skim_Wmunu/res/histo_dataonly_TkFromSta_*.root tp_from_skim_TTbar/res/histo_dataonly_TkFromSta_*.root

histo_dataonly_HltFromGlb.root: tp_from_skim_Zmumu/res/histo_dataonly_HltFromGlb_1.root tp_from_skim_Wmunu/res/histo_dataonly_HltFromGlb_1.root tp_from_skim_TTbar/res/histo_dataonly_HltFromGlb_1.root
	hadd -f histo_dataonly_HltFromGlb.root tp_from_skim_Zmumu/res/histo_dataonly_HltFromGlb_*.root tp_from_skim_Wmunu/res/histo_dataonly_HltFromGlb_*.root tp_from_skim_TTbar/res/histo_dataonly_HltFromGlb_*.root

tp_from_skim_Zmumu/res/CMSSW_1.stdout:
	crab -status -c tp_from_skim_Zmumu

tp_from_skim_Wmunu/res/CMSSW_1.stdout:
	crab -status -c tp_from_skim_Wmunu

tp_from_skim_TTbar/res/CMSSW_1.stdout:
	crab -status -c tp_from_skim_TTbar

tp_from_skim_Zmumu/res/CMSSW_1.stderr:
	crab -getoutput -c tp_from_skim_Zmumu || /bin/true

tp_from_skim_Wmunu/res/CMSSW_1.stderr:
	crab -getoutput -c tp_from_skim_Wmunu || /bin/true

tp_from_skim_TTbar/res/CMSSW_1.stderr:
	crab -getoutput -c tp_from_skim_TTbar || /bin/true

tp_from_skim_Zmumu: tp_from_skim_Zmumu.crab
	rm -rf tp_from_skim_Zmumu; crab -create -submit -cfg tp_from_skim_Zmumu.crab

tp_from_skim_Wmunu: tp_from_skim_Wmunu.crab
	rm -rf tp_from_skim_Wmunu; crab -create -submit -cfg tp_from_skim_Wmunu.crab

tp_from_skim_TTbar: tp_from_skim_TTbar.crab
	rm -rf tp_from_skim_TTbar; crab -create -submit -cfg tp_from_skim_TTbar.crab

