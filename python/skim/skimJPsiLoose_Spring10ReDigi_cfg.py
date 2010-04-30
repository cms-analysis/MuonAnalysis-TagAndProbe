## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/DiMuonSkim_ppMuX_Spring10_REDIGI_START3X_V26_S09_GEN-SIM-RECO.root' ]  # accessible from anywhere in CERN
process.jpsiSkimOut.fileName = "skimJPsiLoose_Summer10ReDigi.root"

from  MuonAnalysis.TagAndProbe.skim.skimJPsi_cff import Spring10ReDigi_Trigger;
Spring10ReDigi_Trigger(process)

setattr(process, '_Process__name', 'ReSkim')
