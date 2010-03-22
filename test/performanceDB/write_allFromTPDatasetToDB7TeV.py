import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.CondDBCommon.connect = 'sqlite_file:MuonPhysicsPerformance7TeV.db'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource")

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    timetype = cms.untracked.string("runnumber"),                                          
    toPut = cms.VPSet(
    cms.PSet(
    record = cms.string('GLBMUJPSI_TEST7TEV_TABLE'),
    tag = cms.string('GLBMUJPSI_TEST7TEV_TABLE'),
    label = cms.string('GLBMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUJPSI_TEST7TEV_WP'),
    tag = cms.string('GLBMUJPSI_TEST7TEV_WP'),
    label = cms.string('GLBMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUJPSI_TEST7TEV_TABLE'),
    tag = cms.string('TRKEFFMUJPSI_TEST7TEV_TABLE'),
    label = cms.string('TRKEFFMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUJPSI_TEST7TEV_WP'),
    tag = cms.string('TRKEFFMUJPSI_TEST7TEV_WP'),
    label = cms.string('TRKEFFMUJPSI_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRGMUJPSI_TEST7TEV_TABLE'),
    tag = cms.string('TRGMUJPSI_TEST7TEV_TABLE'),
    label = cms.string('TRGMUJPSI_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMUJPSI_TEST7TEV_WP'),
    tag = cms.string('TRGMUJPSI_TEST7TEV_WP'),
    label = cms.string('TRGMUJPSI_TEST7TEV_WP')
    ),    
    cms.PSet(
    record = cms.string('GLBMUJPSICAL_TEST7TEV_TABLE'),
    tag = cms.string('GLBMUJPSICAL_TEST7TEV_TABLE'),
    label = cms.string('GLBMUJPSICAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUJPSICAL_TEST7TEV_WP'),
    tag = cms.string('GLBMUJPSICAL_TEST7TEV_WP'),
    label = cms.string('GLBMUJPSICAL_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('GLBMUZ_TEST7TEV_TABLE'),
    tag = cms.string('GLBMUZ_TEST7TEV_TABLE'),
    label = cms.string('GLBMUZ_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUZ_TEST7TEV_WP'),
    tag = cms.string('GLBMUZ_TEST7TEV_WP'),
    label = cms.string('GLBMUZ_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUZ_TEST7TEV_TABLE'),
    tag = cms.string('TRKEFFMUZ_TEST7TEV_TABLE'),
    label = cms.string('TRKEFFMUZ_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRKEFFMUZ_TEST7TEV_WP'),
    tag = cms.string('TRKEFFMUZ_TEST7TEV_WP'),
    label = cms.string('TRKEFFMUZ_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRGMUZ_TEST7TEV_TABLE'),
    tag = cms.string('TRGMUZ_TEST7TEV_TABLE'),
    label = cms.string('TRGMUZ_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMUZ_TEST7TEV_WP'),
    tag = cms.string('TRGMUZ_TEST7TEV_WP'),
    label = cms.string('TRGMUZ_TEST7TEV_WP')
    ),    
    cms.PSet(
    record = cms.string('GLBMUZCAL_TEST7TEV_TABLE'),
    tag = cms.string('GLBMUZCAL_TEST7TEV_TABLE'),
    label = cms.string('GLBMUZCAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('GLBMUZCAL_TEST7TEV_WP'),
    tag = cms.string('GLBMUZCAL_TEST7TEV_WP'),
    label = cms.string('GLBMUZCAL_TEST7TEV_WP')
    ),

    #PAG-specific selections
    cms.PSet(
    record = cms.string('TRGMUJPSI_JPSIANAL_TEST7TEV_TABLE'),
    tag = cms.string('TRGMUJPSI_JPSIANAL_TEST7TEV_TABLE'),
    label = cms.string('TRGMUJPSI_JPSIANAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMUJPSI_JPSIANAL_TEST7TEV_WP'),
    tag = cms.string('TRGMUJPSI_JPSIANAL_TEST7TEV_WP'),
    label = cms.string('TRGMUJPSI_JPSIANAL_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('MUJPSI_JPSIGLBANAL_TEST7TEV_TABLE'),
    tag = cms.string('MUJPSI_JPSIGLBANAL_TEST7TEV_TABLE'),
    label = cms.string('MUJPSI_JPSIGLBANAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('MUJPSI_JPSIGLBANAL_TEST7TEV_WP'),
    tag = cms.string('MUJPSI_JPSIGLBANAL_TEST7TEV_WP'),
    label = cms.string('MUJPSI_JPSIGLBANAL_TEST7TEV_WP')
    ),    
    cms.PSet(
    record = cms.string('MUJPSI_JPSITKMANAL_TEST7TEV_TABLE'),
    tag = cms.string('MUJPSI_JPSITKMANAL_TEST7TEV_TABLE'),
    label = cms.string('MUJPSI_JPSITKMANAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('MUJPSI_JPSITKMANAL_TEST7TEV_WP'),
    tag = cms.string('MUJPSI_JPSITKMANAL_TEST7TEV_WP'),
    label = cms.string('MUJPSI_JPSITKMANAL_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('MUJPSI_JPSIPLUSMUANAL_TEST7TEV_TABLE'),
    tag = cms.string('MUJPSI_JPSIPLUSMUANAL_TEST7TEV_TABLE'),
    label = cms.string('MUJPSI_JPSIPLUSMUANAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('MUJPSI_JPSIPLUSMUANAL_TEST7TEV_WP'),
    tag = cms.string('MUJPSI_JPSIPLUSMUANAL_TEST7TEV_WP'),
    label = cms.string('MUJPSI_JPSIPLUSMUANAL_TEST7TEV_WP')
    ),    
    cms.PSet(
    record = cms.string('MUJPSI_BEXCLANAL_TEST7TEV_TABLE'),
    tag = cms.string('MUJPSI_BEXCLANAL_TEST7TEV_TABLE'),
    label = cms.string('MUJPSI_BEXCLANAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('MUJPSI_BEXCLANAL_TEST7TEV_WP'),
    tag = cms.string('MUJPSI_BEXCLANAL_TEST7TEV_WP'),
    label = cms.string('MUJPSI_BEXCLANAL_TEST7TEV_WP')
    ),
    cms.PSet(
    record = cms.string('TRGMUZ_WMUNUANAL_TEST7TEV_TABLE'),
    tag = cms.string('TRGMUZ_WMUNUANAL_TEST7TEV_TABLE'),
    label = cms.string('TRGMUZ_WMUNUANAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('TRGMUZ_WMUNUANAL_TEST7TEV_WP'),
    tag = cms.string('TRGMUZ_WMUNUANAL_TEST7TEV_WP'),
    label = cms.string('TRGMUZ_WMUNUANAL_TEST7TEV_WP')
    ),
    cms.PSet(    
    record = cms.string('MUZ_WMUNUANAL_TEST7TEV_TABLE'),
    tag = cms.string('MUZ_WMUNUANAL_TEST7TEV_TABLE'),
    label = cms.string('MUZ_WMUNUANAL_TEST7TEV_TABLE')
    ),
    cms.PSet(
    record = cms.string('MUZ_WMUNUANAL_TEST7TEV_WP'),
    tag = cms.string('MUZ_WMUNUANAL_TEST7TEV_WP'),
    label = cms.string('MUZ_WMUNUANAL_TEST7TEV_WP')
    )))    
        
                      



process.mywriter = cms.EDFilter("PhysicsPerformanceDBWriterFromTPDataset",

                                # For each table to be loaded, set the name of the input T/P file, histogram, algorithm, and cut (if any)
                                inputHistoFiles = cms.vstring('myplots.root'),
                                
                                inputHistogramNames = cms.vstring('fit_eff_pt_eta'),
                                
                                inputAlgorithmNames = cms.vstring('GlobalMuonFromTrackerTrackJpsi'),
                                
                                inputDiscriminatorCuts = cms.vdouble(0.0),
                                
                                
                                # For each table to be loaded, set the payload and working point record names as
                                # defined above in the PoolDBOutputService
                                RecordPayloads = cms.vstring('GLBMUJPSI_TEST7TEV_TABLE'),
                                
                                RecordWPs = cms.vstring('GLBMUJPSI_TEST7TEV_WP'),

                                # Set the type of data to be stored, the binning variables, and the IOV
                                # These are currently assumed to be the same for all tables loaded in a single job
                                inputResultTypes = cms.vstring('efficiency','efficiency_symerr'),
                                inputBinningVariables = cms.vstring('pt','eta'),                                
                                IOVBegin = cms.uint64(1),
                                IOVEnd = cms.uint64(50)
                                )


process.p = cms.Path(process.mywriter)


