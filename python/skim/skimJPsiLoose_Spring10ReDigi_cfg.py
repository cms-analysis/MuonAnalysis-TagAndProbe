## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

# note: both are accessible from anywhere in CERN
#process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/DiMuonSkim_ppMuX_Spring10_REDIGI_START3X_V26_S09_GEN-SIM-RECO.root' ]  
process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/ppMuX_Spring10_REDIGI_START3X_V26_S09_GEN-SIM-RECO_C0AC7DEB-8144-DF11-A1E1-00304867D838.root' ]  
process.jpsiSkimOut.fileName = "skimJPsiLoose_Spring10ReDigi.root"

from  MuonAnalysis.TagAndProbe.skim.skimJPsi_cff import Spring10ReDigi_Trigger;
Spring10ReDigi_Trigger(process)

#setattr(process, '_Process__name', 'ReSkim')
