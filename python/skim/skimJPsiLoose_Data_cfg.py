## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

#process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/MuonSkimTight_Dec19th_RECO.root' ] # accessible from anywhere in CERN, ~120Mb, 1674 events
#process.source.fileNames = [ 'rfio:/castor/cern.ch/user/g/gpetrucc/7TeV/DATA/Cocktail_runs1326xx.root' ]
#process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/Cocktail_runs1326xx.root' ]
process.source.fileNames = [ '/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/149/D6383A41-F85A-DF11-A2A7-0030487C8CBE.root' ]
process.jpsiSkimOut.fileName = "skimJPsiLoose_Data.root"

# Replace collision filter for MC with the one for Data
process.globalReplace('collisionFilter', process.collisionFilterData)

# Take away MC-related stuff matching
process.slimAOD.remove(process.genMuons)

#### Discard immediately events which don't have a global muon.
#### This speeds up the skim significantly, about x5, when running on datasets like MinimumBias that are low on muons.
# declare the filter
process.fastPreFilter = cms.EDFilter("TrackCountFilter", src = cms.InputTag("globalMuons"), minNumber = cms.uint32(1))
for pn in process.paths_().keys():          # loop on all paths names
    p = getattr(process,pn)                 # get path object
    p._seq = process.fastPreFilter + p._seq # prepend prefilter
