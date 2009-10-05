create: tp_from_skim_JPsiMuMu tp_from_skim_ppMuX
status: tp_from_skim_JPsiMuMu/res/CMSSW_1.stdout tp_from_skim_ppMuX/res/CMSSW_1.stdout
getoutput: tp_from_skim_JPsiMuMu/res/CMSSW_1.stderr tp_from_skim_ppMuX/res/CMSSW_1.stderr
fit: fit_result_GlbFromCal.root


fit_result_GlbFromCal.root: histo_output_GlbFromCal.root histo_output_GlbFromTk.root  histo_output_HltFromGlb.root histo_output_TkFromSta.root
	cmsRun fit.py

histo_output_GlbFromCal.root: tp_from_skim_JPsiMuMu/res/histo_output_GlbFromCal_1.root tp_from_skim_ppMuX/res/histo_output_GlbFromCal_1.root
	hadd -f histo_output_GlbFromCal.root tp_from_skim_JPsiMuMu/res/histo_output_GlbFromCal_*.root tp_from_skim_ppMuX/res/histo_output_GlbFromCal_*.root 

histo_output_GlbFromTk.root: tp_from_skim_JPsiMuMu/res/histo_output_GlbFromTk_1.root tp_from_skim_ppMuX/res/histo_output_GlbFromTk_1.root
	hadd -f histo_output_GlbFromTk.root tp_from_skim_JPsiMuMu/res/histo_output_GlbFromTk_*.root tp_from_skim_ppMuX/res/histo_output_GlbFromTk_*.root

histo_output_HltFromGlb.root: tp_from_skim_JPsiMuMu/res/histo_output_HltFromGlb_1.root tp_from_skim_ppMuX/res/histo_output_HltFromGlb_1.root
	hadd -f histo_output_HltFromGlb.root tp_from_skim_JPsiMuMu/res/histo_output_HltFromGlb_*.root tp_from_skim_ppMuX/res/histo_output_HltFromGlb_*.root

histo_output_TkFromSta.root: tp_from_skim_JPsiMuMu/res/histo_output_TkFromSta_1.root tp_from_skim_ppMuX/res/histo_output_TkFromSta_1.root
	hadd -f histo_output_TkFromSta.root tp_from_skim_JPsiMuMu/res/histo_output_TkFromSta_*.root tp_from_skim_ppMuX/res/histo_output_TkFromSta_*.root

tp_from_skim_JPsiMuMu/res/CMSSW_1.stdout:
	crab -status -c tp_from_skim_JPsiMuMu

tp_from_skim_ppMuX/res/CMSSW_1.stdout:
	crab -status -c tp_from_skim_ppMuX

tp_from_skim_JPsiMuMu/res/CMSSW_1.stderr:
	crab -getoutput -c tp_from_skim_JPsiMuMu 

tp_from_skim_ppMuX/res/CMSSW_1.stderr:
	crab -getoutput -c tp_from_skim_ppMuX

tp_from_skim_JPsiMuMu: tp_from_skim_JPsiMuMu.crab
	rm -rf tp_from_skim_JPsiMuMu; crab -create -submit -cfg tp_from_skim_JPsiMuMu.crab

tp_from_skim_ppMuX: tp_from_skim_ppMuX.crab
	rm -rf tp_from_skim_ppMuX; crab -create -submit -cfg tp_from_skim_ppMuX.crab

