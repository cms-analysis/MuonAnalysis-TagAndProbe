## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

#process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/MuonSkimTight_Dec19th_RECO.root' ] # accessible from anywhere in CERN, ~120Mb, 1674 events
#process.source.fileNames = [ 'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/DATA/Cocktail_runs1326xx.root' ]
#process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/Cocktail_runs1326xx.root' ]
process.source.fileNames = [ '/store/data/Run2010A/Mu/RECO/v2/000/136/119/C25E9FFB-2C68-DF11-9134-0030487CD7CA.root' ]
process.jpsiSkimOut.fileName = "skimJPsiLoose_Data.root"

# Replace collision filter for MC with the one for Data
process.globalReplace('collisionFilter', process.collisionFilterData)

# Take away MC-related stuff matching
process.slimAOD.remove(process.genMuons)
