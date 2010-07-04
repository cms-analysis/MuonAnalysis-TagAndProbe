## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

process.source.fileNames = [ 'root://pcmssd12.cern.ch//data/gpetrucc/7TeV/jpsi/CS_Onia-Jun14thSkim_v1_RAW-RECO_run136082_443584C2-B27E-DF11-9E13-0017A477001C.root' ]
process.jpsiSkimOut.fileName = "skimJPsiLoose_Data_RAW.root"

# Replace collision filter for MC with the one for Data
process.globalReplace('collisionFilter', process.collisionFilterData)

from MuonAnalysis.TagAndProbe.skim.skimJPsi_cff import Add_CSCTF_Matching
Add_CSCTF_Matching(process, isRealData=True)

## Now, for this we only use Mu+Mu events
process.jpsiSkimOut.SelectEvents.SelectEvents = [ "Skim_jpsiMu" ]
process.patMuonSequence.remove(process.mergedMuons)
process.muonL1Info.src = 'muons'
process.patMuonsWithoutTrigger.muonSource = 'muons'
process.muonL1MatchExtended.muons = "muons"

# Take away MC-related stuff matching
process.slimAOD.remove(process.genMuons)
