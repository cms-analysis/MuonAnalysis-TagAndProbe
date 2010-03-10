## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

from HeavyFlavorAnalysis.Onia2MuMu.goodLumiSectionList_cfi import goodLumisToProcess
good900GeVLumis = cms.untracked.VLuminosityBlockRange()
good2TeVLumis   = cms.untracked.VLuminosityBlockRange()
good900GeVLumis += [ x for x in goodLumisToProcess if x.find('124120') == -1]
good2TeVLumis   += [ x for x in goodLumisToProcess if x.find('124120') != -1]

#### The Run 124120 at 2.36 TeV
#process.source.fileNames = [
#    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/Muon_skim-Dec19thSkim_341_v2/0004/F6796FFE-94ED-DE11-838F-0026189438B1.root',
#    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/Muon_skim-Dec19thSkim_341_v2/0004/0ED9D2F6-94ED-DE11-A5D2-002618943944.root',
#]
#process.out.fileName = "skimJPsiLoose_Data_R124120.root"
#process.source.lumisToProcess = good2TeVLumis

#### The 900 GeV runs (from a subskim of Muon_skim-Dec19thSkim requiring BPTX AND and one reco::Muon with a tracker track)
process.source.fileNames = [ 'root:://pcmssd12.cern.ch//data/gpetrucc/MuonSkimTight_Dec19th_RECO.root' ] # accessible from anywhere in CERN, ~120Mb, 1674 events
process.source.lumisToProcess = good900GeVLumis
process.out.fileName = "skimJPsiLoose_Data_900GeV.root"

# Replace collision filter for MC with the one for Data
process.globalReplace('collisionFilter', process.collisionFilterData)

# Take away MC-related stuff matching
from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff import onia2MuMu_isNotMC
onia2MuMu_isNotMC(process)
process.slimAOD.remove(process.genMuons)
