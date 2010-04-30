## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsi_cfg import cms, process

## Loosen the cut on one muon: remove any requirement # the trigger requirement
#process.goldenMuons.cut = "isGlobalMuon && pt > %f && abs(eta) < %f" % (ptTag, etaTag)       # no trigger
process.tags1Mu.cut = "isGlobalMuon"                                                     # no cuts
process.tags2Mu.cut = "isGlobalMuon"                                                     # no cuts
#process.goldenMuons.cut = "isGlobalMuon || numberOfMatches('SegmentAndTrackArbitration') > 0" # also TM 

process.Check_Collisions = cms.Path(process.collisionFilter)

## Widen the mass cuts
process.jpsiMu.cut  = "2.5 < mass < 4"
process.jpsiTk.cut  = "2.5 < mass < 4"
process.jpsiSta.cut = "2 < mass < 20"
## and keep the di-muon objects for further inspection
#process.out.outputCommands += [ 'keep *_jpsiTk_*_*', 'keep *_jpsiMu_*_*', 'keep *_jpsiSta_*_*' ]


## Test on background
## file: /store/mc/Summer09/ppMuXLoose/GEN-SIM-RECO/STARTUP3X_V8I_900GeV-v1/0005/F83A17E0-9AEB-DE11-844A-00E08178C015.root
##  aka  root://pcmssd12.cern.ch//data/mc/900GeV/ppMuXLoose_STARTUP3X_V8I_RECO/F83A17E0-9AEB-DE11-844A-00E08178C015.root
##  aka  file:/data/mc/900GeV/ppMuXLoose_STARTUP3X_V8I_RECO/F83A17E0-9AEB-DE11-844A-00E08178C015.root 
## process.source.fileNames = [ 'file:/data/mc/900GeV/ppMuXLoose_STARTUP3X_V8I_RECO/F83A17E0-9AEB-DE11-844A-00E08178C015.root' ]
#process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/mc/900GeV/ppMuXLoose_STARTUP3X_V8I_RECO/F83A17E0-9AEB-DE11-844A-00E08178C015.root' ]
process.source.fileNames =['root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/BpToJPsiMuMu_Spring10_START3X_V26_S09_GEN-SIM-RECO_0003_E477AF59-764D-DF11-9839-002618943902.root']
#
process.jpsiSkimOut.fileName = "skimJPsiLoose_ppxMuXLoose900GeV.root"

