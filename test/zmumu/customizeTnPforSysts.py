import FWCore.ParameterSet.Config as cms

######### Define here the systematics ############
ALLSYSTS = [ 
    ("Mass_binning"     , [30, 50]                ), ## N. bins for mass fit (default: 40) 
    ("Mass_range"       , [[60, 130], [70, 120]]  ), ## [low mass, high mass] (default: [70, 130]) 
    ("Signal_shape"     , ["voigtPlusExpo"]       ), ## (default: vpvPlusExpo) 
    ("Background_shape" , ["vpvPlusCheb"]         ), ## (default: vpvPlusExpo) 
    ("Background_level" , [[25, 0.12], [20, 999]] )  ## [min tag pT, max tag iso] (default: [21, 0.2]) 
    ] 
################################################## 


def customizeTnPforSysts(process): 
    for modname in process.analyzerNames().split(): 
        if modname == "TnP_MuonID": continue 
        mod = getattr(process, modname) 
        if mod.type_() == "TagProbeFitTreeAnalyzer": 
            setattr(process, modname, mod)        
            setattr(process, "run_"+modname, cms.Path(mod))

            for SYST,SYSTVEC in ALLSYSTS: 
                systcnt = 0 
                for SYSTI in SYSTVEC: 
                    newoutname = mod.OutputFileName.value()
                    newoutname = newoutname.replace( ".root", str("_%s_%s.root" % (SYST,str(systcnt))) )
                    modnew = mod.clone(OutputFileName = newoutname)
                    if SYST == "Mass_binning": 
                        modnew.binnedFit = cms.bool(True) 
                        modnew.binsForFit = cms.uint32(SYSTI) 
                    if SYST == "Mass_range": 
                        modnew.Variables.mass[1] = str(SYSTI[0]) 
                        modnew.Variables.mass[2] = str(SYSTI[1]) 
                    if SYST == "Signal_shape" or SYST == "Background_shape": 
                        for effi in modnew.Efficiencies.parameters_().keys(): 
                            effpset = getattr(modnew.Efficiencies, effi) 
                            if isinstance(effpset, cms.PSet): 
                                effpset.BinToPDFmap = cms.vstring(SYSTI) 
                    if SYST == "Background_level": 
                        for effi in modnew.Efficiencies.parameters_().keys(): 
                            effpset = getattr(modnew.Efficiencies, effi) 
                            if isinstance(effpset, cms.PSet): 
                                if hasattr(effpset.BinnedVariables, "tag_pt"): 
                                    effpset.BinnedVariables.tag_pt[0] = SYSTI[0] 
                                else: 
                                    effpset.BinnedVariables.tag_pt = cms.vdouble(SYSTI[0], 500) 
                                if hasattr(effpset.BinnedVariables, "tag_combRelIsoPF04dBeta"): 
                                    effpset.BinnedVariables.tag_combRelIsoPF04dBeta[1] = SYSTI[1] 
                                else: 
                                    effpset.BinnedVariables.tag_combRelIsoPF04dBeta = cms.vdouble(-0.5, SYSTI[1]) 

                    newmodname = modname + "_" + SYST + "_" + str(systcnt) 
                    setattr(process, newmodname, modnew)        
                    setattr(process, "run_"+newmodname, cms.Path(modnew))
                    systcnt = systcnt+1

    return process
