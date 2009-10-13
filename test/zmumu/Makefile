### DEFINE THE SAMPLES
SAMPLES = Zmumu Wmunu TTbar QCD_Pt15 QCD_Pt30 InclusiveMu15

### DEFINE THE MEASUREMENTS TO DO
FITS = GlbFromCal GlbFromTk HltFromGlb TkFromSta HltFromWmunu MuFromTkWmunu ## MuFromTkZmumu

##     ____ ____      _    ____        _          __  __ 
##    / ___|  _ \    / \  | __ )   ___| |_ _   _ / _|/ _|
##   | |   | |_) |  / _ \ |  _ \  / __| __| | | | |_| |_ 
##   | |___|  _ <  / ___ \| |_) | \__ \ |_| |_| |  _|  _|
##    \____|_| \_\/_/   \_\____/  |___/\__|\__,_|_| |_|  
##                                                       
##   

### DEFINE THE NAME OF THE CRAB DIRS
### Simply add 'tp_from_skim_' before the name of the sample.
CRABDIRS = $(patsubst %,tp_from_skim_%,$(SAMPLES))


### ------------------------------------
### 'make create' will create and submit the jobs
create: $(CRABDIRS)

## each crabdir has prerequisites the crab file and the python file, and is created with crab -create -submit
$(CRABDIRS): %: %.crab tp_from_skim_dataOnly.py
	rm -rf $@ || /bin/true
	crab -cfg $@.crab -create

submit: $(CRABDIRS)
	for D in $(CRABDIRS); do crab -c $$D -submit 2>&1 | tee status.$$D; done


### ------------------------------------
### 'make status' will issue crab -status
STDOUTS = $(patsubst %,%/res/CMSSW_1.stdout, $(CRABDIRS))
status: $(STDOUTS)

allstatus:
	for D in $(CRABDIRS); do crab -c $$D -status 2>&1 | tee status.$$D; done

## this rule pretends that *.stdout is made using crab -c ... -status
$(STDOUTS): %/res/CMSSW_1.stdout: %
	crab -c $< -status 2>&1 | tee status.$(patsubst %,%.status,$<)

### ------------------------------------
### 'make getoutput' will issue crab -status
STDERRS = $(patsubst %,%/res/CMSSW_1.stderr, $(CRABDIRS))
getoutput: $(STDERRS)

## this rule defines that *.stderr is made using crab -c ... -get
$(STDERRS): %/res/CMSSW_1.stderr: %
	crab -c $< -getoutput 

##    _____ ___ _____       _          __  __ 
##   |  ___|_ _|_   _|  ___| |_ _   _ / _|/ _|
##   | |_   | |  | |   / __| __| | | | |_| |_ 
##   |  _|  | |  | |   \__ \ |_| |_| |  _|  _|
##   |_|   |___| |_|   |___/\__|\__,_|_| |_|  
##                                            
##   
HISTOS  = $(patsubst %,histo_dataonly_%.root, $(FITS))
RESULTS = $(patsubst %,fit_result_%.root, $(FITS))

## 'make fit' just means make all the results
fit: $(RESULTS)

## each fit result depends on histo and fit.py
$(RESULTS): fit_result_%.root: histo_dataonly_%.root fit.py
	cmsRun fit.py

## Although it's not really correct, we make each histo depend on all other histos, not just those of its fit
## We use 'wildcard' so that we depend only on histograms which do exist 
## (this allows to test the fit even on a partial result)
ALLCRABHISTS = $(wildcard $(patsubst %,%/res/histo_*.root,$(CRABDIRS)))
$(HISTOS): histo_%.root: $(ALLCRABHISTS)
	hadd -f $@ $(wildcard $(patsubst %,%/res/histo_$*_*.root,$(CRABDIRS)))

histos: $(HISTOS)  $(ALLCRABHISTS)
