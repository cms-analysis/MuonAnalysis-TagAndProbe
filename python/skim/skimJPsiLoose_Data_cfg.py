## Load everything from the main cfg file
from MuonAnalysis.TagAndProbe.skim.skimJPsiLoose_cfg import cms, process

# Starting from the muon skim /MinimumBias/BeamCommissioning09-Muon_skim-Dec19thSkim_341_v2/RAW-RECO, run 124120
process.source.fileNames = [
    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/Muon_skim-Dec19thSkim_341_v2/0004/F6796FFE-94ED-DE11-838F-0026189438B1.root',
    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/Muon_skim-Dec19thSkim_341_v2/0004/0ED9D2F6-94ED-DE11-A5D2-002618943944.root',
]
process.out.fileName = "skimJPsiLoose_Data_R124120.root"

# Replace collision filter for MC with the one for Data
process.globalReplace('collisionFilter', process.collisionFilterData)

# Take away MC-related stuff matching
from HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cff import onia2MuMu_isNotMC
onia2MuMu_isNotMC(process)
process.slimAOD.remove(process.genMuons)
