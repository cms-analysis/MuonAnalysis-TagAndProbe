## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

# note: both are accessible from anywhere in CERN
#process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/DiMuonSkim_ppMuX_Spring10_REDIGI_START3X_V26_S09_GEN-SIM-RECO.root' ]  
process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/JPsiMuMu_Summer10_REDIGI_START37_V5_S09_GEN-SIM-RECO_AE944BC3-5185-DF11-8906-0022198E201A.root' ]  
process.jpsiSkimOut.fileName = "skimJPsiLoose_Summer10_37X.root"

from  MuonAnalysis.TagAndProbe.skim.skimJPsi_cff import changeTriggerProcessName;
changeTriggerProcessName(process, 'REDIGI37X', process.patTrigger.processName.value())

#setattr(process, '_Process__name', 'ReSkim')
