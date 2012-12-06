#!/usr/bin/env python
from sys import argv
from re import *
import os, os.path

print argv
if len(argv) != 4: raise RuntimeError, "Usage: mkjobs.py py crabcfg datasets"
pyname   = argv[1]
py       = "".join([l for l in open(argv[1],'r')])
cfgname  = argv[2]
cfg      = "".join([l for l in open(argv[2],'r')])
datasets = open(argv[3])

jsondir = os.environ['CMSSW_BASE']+"/src/MuonAnalysis/TagAndProbe/test/zmumu/production/";

def cfgsub(cfg, dict):
    for key, value in dict.iteritems():
        cfg = sub(compile("^(\\s*)%s\\s*=\\s*\\S+" % key, MULTILINE), "\\1%s = %s" % (key,value), cfg)
    return cfg

for row in datasets:
    if row[0] == "#": continue
    (nick, dataset, tag, json, jobs) = row.split()
    pyout = pyname.replace(".py", ".%s.py" % nick)
    print "Will prepare python %s to run with global tag %s" % (pyout, tag)
    mypy = cfgsub(py, {"process.GlobalTag.globaltag":("'%s::All'"%tag)} )
    pyhandle = open(pyout, "w")
    pyhandle.write(mypy)
    pyhandle.close()    
    cfgout = os.path.basename(cfgname.replace(".cfg", ".%s.cfg" % nick))
    print "Will prepare crab cfg %s to run %s on dataset %s (%s jobs, json %s)" % (cfgout, pyout, dataset, jobs, json)
    mycfg = cfgsub(cfg, {'number_of_jobs':jobs, 'datasetpath':dataset, 'lumi_mask':jsondir+json, 'pset':pyout} )
    if json in ["-","none","MC"]:
        mycfg = sub(compile("^lumi_mask\s*=\s*\S+",MULTILINE), "", mycfg.replace("total_number_of_lumis","total_number_of_events"))
    cfghandle = open(cfgout, "w")
    cfghandle.write(mycfg)
    cfghandle.close()    
