import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Isolation.tools_cfi import isoDepositReplace

def load_muonPFiso_sequence(proc, seq_name, algo, coneR, src, src_charged_hadron='', src_neutral_hadron='', src_photon='', src_charged_pileup='', 
                            veto_charged_hadron='Threshold(0.0)', veto_neutral_hadron='Threshold(0.5)', veto_photon='Threshold(0.5)', veto_charged_pileup='Threshold(0.5)' ):

    doCH, doNH, doPh, doPU = False, False, False, False
    if src_charged_hadron != '': doCH = True
    if src_neutral_hadron != '': doNH = True
    if src_photon         != '': doPh = True
    if src_charged_pileup != '': doPU = True

    iso_seq = cms.Sequence()

    if doCH:
        setattr(proc, 'muPFIsoDepositCH'+algo, isoDepositReplace(src, src_charged_hadron))
        iso_seq += getattr(proc, 'muPFIsoDepositCH'+algo)

    if doNH:
        setattr(proc, 'muPFIsoDepositNH'+algo, isoDepositReplace(src, src_neutral_hadron))
        iso_seq += getattr(proc, 'muPFIsoDepositNH'+algo)

    if doPh:
        setattr(proc, 'muPFIsoDepositPh'+algo, isoDepositReplace(src, src_photon))
        iso_seq += getattr(proc, 'muPFIsoDepositPh'+algo)

    if doPU:
        setattr(proc, 'muPFIsoDepositPU'+algo, isoDepositReplace(src, src_charged_pileup))
        iso_seq += getattr(proc, 'muPFIsoDepositPU'+algo)

    iso_vals_seq = cms.Sequence()

    if doCH:
        setattr(proc, 'muPFIsoValueCH'+algo,
          cms.EDProducer('CandIsolatorFromDeposits',
            deposits = cms.VPSet(
              cms.PSet(
                src = cms.InputTag('muPFIsoDepositCH'+algo),
                deltaR = cms.double(coneR),
                weight = cms.string('1'),
                vetos = cms.vstring('0.0001',veto_charged_hadron),
                skipDefaultVeto = cms.bool(True),
                mode = cms.string('sum')
              )
            )
          )
        )
        iso_vals_seq += getattr(proc, 'muPFIsoValueCH'+algo)

    if doNH:
        setattr(proc, 'muPFIsoValueNH'+algo,
          cms.EDProducer('CandIsolatorFromDeposits',
            deposits = cms.VPSet(
              cms.PSet(
                src = cms.InputTag('muPFIsoDepositNH'+algo),
                deltaR = cms.double(coneR),
                weight = cms.string('1'),
                vetos = cms.vstring('0.01', veto_neutral_hadron ),
                skipDefaultVeto = cms.bool(True),
                mode = cms.string('sum')
              )
            )
          )
        )
        iso_vals_seq += getattr(proc, 'muPFIsoValueNH'+algo)

    if doPh:
        setattr(proc, 'muPFIsoValuePh'+algo,
          cms.EDProducer('CandIsolatorFromDeposits',
            deposits = cms.VPSet(
              cms.PSet(
                src = cms.InputTag('muPFIsoDepositPh'+algo),
                deltaR = cms.double(coneR),
                weight = cms.string('1'),
                vetos = cms.vstring('0.01', veto_photon),
                skipDefaultVeto = cms.bool(True),
                mode = cms.string('sum')
              )
            )
          )
        )
        iso_vals_seq += getattr(proc, 'muPFIsoValuePh'+algo)

    if doPU:
        setattr(proc, 'muPFIsoValuePU'+algo,
          cms.EDProducer('CandIsolatorFromDeposits',
            deposits = cms.VPSet(
              cms.PSet(
                src = cms.InputTag('muPFIsoDepositPU'+algo),
                deltaR = cms.double(coneR),
                weight = cms.string('1'),
                vetos = cms.vstring('0.01', veto_charged_pileup ),
                skipDefaultVeto = cms.bool(True),
                mode = cms.string('sum')
              )
            )
          )
        )
        iso_vals_seq += getattr(proc, 'muPFIsoValuePU'+algo)

    iso_seq *= iso_vals_seq

    setattr(proc, seq_name, iso_seq)


def load_pfParticle_sequence(proc, withLeptons):

    from CommonTools.PileupAlgos.Puppi_cff import puppi


    if withLeptons:
        

        puppiR05 = puppi.clone()
        setattr(proc, "puppiR05", puppiR05)


        pfAllPhotonsPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 22"))
        pfAllNeutralHadronsPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
        pfAllChargedHadronsPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05"), cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212"))
        setattr(proc, "pfAllPhotonsPuppi", pfAllPhotonsPuppi)
        setattr(proc, "pfAllNeutralHadronsPuppi", pfAllNeutralHadronsPuppi)
        setattr(proc, "pfAllChargedHadronsPuppi", pfAllChargedHadronsPuppi)
        
        pfParticlesForPUPPI = cms.Sequence()
        pfParticlesForPUPPI += getattr(proc, "puppiR05")
        pfParticlesForPUPPI += getattr(proc, "pfAllPhotonsPuppi")
        pfParticlesForPUPPI += getattr(proc, "pfAllNeutralHadronsPuppi")
        pfParticlesForPUPPI += getattr(proc, "pfAllChargedHadronsPuppi")

        return pfParticlesForPUPPI

    else :

        packedPFCandidatesNoLep = cms.EDFilter("CandPtrSelector", src = cms.InputTag("particleFlow"), cut = cms.string("abs(pdgId) != 13 && abs(pdgId) != 11 "))
        setattr(proc, "packedPFCandidatesNoLep", packedPFCandidatesNoLep)

        puppiR05NoLep = puppi.clone(
                                               candName = 'packedPFCandidatesNoLep'
                                               )
        setattr(proc, "puppiR05NoLep", puppiR05NoLep)

        pfAllPhotonsPuppiNoLep = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05NoLep"), cut = cms.string("pdgId == 22"))
        pfAllNeutralHadronsPuppiNoLep = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05NoLep"), cut = cms.string("pdgId == 111 || pdgId == 130 || pdgId == 310 || pdgId == 2112"))
        pfAllChargedHadronsPuppiNoLep = cms.EDFilter("CandPtrSelector", src = cms.InputTag("puppiR05NoLep"), cut = cms.string("pdgId == 211 || pdgId == -211 || pdgId == 321 || pdgId == -321 || pdgId == 999211 || pdgId == 2212 || pdgId == -2212"))
        setattr(proc, "pfAllPhotonsPuppiNoLep", pfAllPhotonsPuppiNoLep)
        setattr(proc, "pfAllNeutralHadronsPuppiNoLep", pfAllNeutralHadronsPuppiNoLep)
        setattr(proc, "pfAllChargedHadronsPuppiNoLep", pfAllChargedHadronsPuppiNoLep)

        pfParticlesForPUPPINoLep = cms.Sequence()
        pfParticlesForPUPPINoLep += getattr(proc, "packedPFCandidatesNoLep")
        pfParticlesForPUPPINoLep += getattr(proc, "puppiR05NoLep")
        pfParticlesForPUPPINoLep += getattr(proc, "pfAllPhotonsPuppiNoLep")
        pfParticlesForPUPPINoLep += getattr(proc, "pfAllNeutralHadronsPuppiNoLep")
        pfParticlesForPUPPINoLep += getattr(proc, "pfAllChargedHadronsPuppiNoLep")

        return pfParticlesForPUPPINoLep

def load_fullPFpuppiIsolation(proc):
    print "coucou"
    proc.pfParticlesPUPPILeptons = load_pfParticle_sequence(proc, True)

    muon_src, cone_size = 'probeMuons', 0.4
    load_muonPFiso_sequence(proc, 'MuonPFIsoSequencePUPPI', algo = 'R04PUPPI',
                            src = muon_src,
                            src_charged_hadron = 'pfAllChargedHadronsPuppi',
                            src_neutral_hadron = 'pfAllNeutralHadronsPuppi',
                            src_photon         = 'pfAllPhotonsPuppi',
                            veto_charged_hadron='Threshold(0.0)',
                            veto_neutral_hadron='Threshold(0.0)',
                            veto_photon='Threshold(0.0)',
                            coneR = cone_size
                            )
                            
    proc.pfParticlesPUPPINoLeptons = load_pfParticle_sequence(proc, False)
    load_muonPFiso_sequence(proc, 'MuonPFIsoSequencePUPPINoLep', algo = 'R04PUPPINoLep',
                            src = muon_src,
                            src_charged_hadron = 'pfAllChargedHadronsPuppiNoLep',
                            src_neutral_hadron = 'pfAllNeutralHadronsPuppiNoLep',
                            src_photon         = 'pfAllPhotonsPuppiNoLep',
                            veto_charged_hadron='Threshold(0.0)',
                            veto_neutral_hadron='Threshold(0.0)',
                            veto_photon='Threshold(0.0)',
                            coneR = cone_size
                            )

    proc.loadPUPPIisoInValueMaps = cms.EDProducer("AnyNumbersToValueMaps",
                                                  collection = cms.InputTag("probeMuons"),
                                                  associations = cms.VInputTag(cms.InputTag("muPFIsoValueCHR04PUPPI"), cms.InputTag("muPFIsoValueNHR04PUPPI"), cms.InputTag("muPFIsoValuePhR04PUPPI"),cms.InputTag("muPFIsoValueCHR04PUPPINoLep"), cms.InputTag("muPFIsoValueNHR04PUPPINoLep"), cms.InputTag("muPFIsoValuePhR04PUPPINoLep")),
    )
    proc.fullPuppIsolationSequence = cms.Sequence(proc.pfParticlesPUPPILeptons+proc.MuonPFIsoSequencePUPPI+proc.pfParticlesPUPPINoLeptons+proc.MuonPFIsoSequencePUPPINoLep+proc.loadPUPPIisoInValueMaps )
 
    return proc.fullPuppIsolationSequence
