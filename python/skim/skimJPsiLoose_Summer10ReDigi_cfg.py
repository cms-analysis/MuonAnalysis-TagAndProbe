## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

# note: both are accessible from anywhere in CERN
#process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/DiMuonSkim_ppMuX_Spring10_REDIGI_START3X_V26_S09_GEN-SIM-RECO.root' ]  
process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/JPsiMuMu_Summer10_REDIGI_START36_V9_S09_GEN-SIM-RECODEBUG_52819337-FF7B-DF11-B0F5-E0CB4E1A118E.root' ]  
process.jpsiSkimOut.fileName = "skimJPsiLoose_Summer10ReDigi.root"

from  MuonAnalysis.TagAndProbe.skim.skimJPsi_cff import Summer10ReDigi_Trigger;
Summer10ReDigi_Trigger(process)

#setattr(process, '_Process__name', 'ReSkim')
