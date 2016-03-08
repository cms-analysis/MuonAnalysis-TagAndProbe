import FWCore.ParameterSet.Config as cms

CentralityVariables = cms.PSet(
	hiBin	= cms.InputTag("centralityBinInfo"),
	hiHF	= cms.InputTag("centralityInfo","HFtowers"),
	hiHFplus	= cms.InputTag("centralityInfo","HFtowersPlus"),
	hiHFminus	= cms.InputTag("centralityInfo","HFtowersMinus"),
	hiHFeta4	= cms.InputTag("centralityInfo","HFtowersTrunc"),
	hiHFplusEta4	= cms.InputTag("centralityInfo","HFtowersPlusTrunc"),
	hiHFminusEta4	= cms.InputTag("centralityInfo","HFtowersMinusTrunc"),
	hiHFhit	= cms.InputTag("centralityInfo","HFhits"),
	hiNpix	= cms.InputTag("centralityInfo","PixelHits"),
	hiNpixelTracks	= cms.InputTag("centralityInfo","PixelTracks"),
	hiNtracks	= cms.InputTag("centralityInfo","Tracks"),
	hiEB	= cms.InputTag("centralityInfo","EB"),
	hiEE	= cms.InputTag("centralityInfo","EE"),
	hiET	= cms.InputTag("centralityInfo","ET"),
)
