#!/usr/bin/env python

import os, string, sys, posix, tokenize, array, getopt

def main(argv):

    tables = ['GLBMUJPSI_OCTXTEST_TABLE','TRKEFFMUJPSI_OCTXTEST_TABLE','TRGMUJPSI_OCTXTEST_TABLE','GLBMUJPSICAL_OCTXTEST_TABLE',
              'GLBMUZ_OCTXTEST_TABLE','TRKEFFMUZ_OCTXTEST_TABLE','TRGMUZ_OCTXTEST_TABLE','GLBMUZCAL_OCTXTEST_TABLE',
              'TRGMUJPSI_JPSIANAL_OCTXTEST_TABLE', 'MUJPSI_JPSIGLBANAL_OCTXTEST_TABLE', 'MUJPSI_JPSITKMANAL_OCTXTEST_TABLE',
              'MUJPSI_JPSIPLUSMUANAL_OCTXTEST_TABLE', 'MUJPSI_BEXCLANAL_OCTXTEST_TABLE']

    wps = ['GLBMUJPSI_OCTXTEST_WP','TRKEFFMUJPSI_OCTXTEST_WP','TRGMUJPSI_OCTXTEST_WP','GLBMUJPSICAL_OCTXTEST_WP',
           'GLBMUZ_OCTXTEST_WP','TRKEFFMUZ_OCTXTEST_WP','TRGMUZ_OCTXTEST_WP','GLBMUZCAL_OCTXTEST_WP',
           'TRGMUJPSI_JPSIANAL_OCTXTEST_WP', 'MUJPSI_JPSIGLBANAL_OCTXTEST_WP', 'MUJPSI_JPSITKMANAL_OCTXTEST_WP',
           'MUJPSI_JPSIPLUSMUANAL_OCTXTEST_WP', 'MUJPSI_BEXCLANAL_OCTXTEST_WP']
    
    effs = ['GlobalMuonFromTrackerTrackJpsi',
            'TrackerTrackFromStandaloneMuonJpsi',
            'TriggerMuonFromGlobalMuonJpsi',
            'GlobalMuonFromCaloMuonJpsi',
            'GlobalMuonFromTrackerTrackZ',
            'TrackerTrackFromStandaloneMuonZ',
            'TriggerMuonFromGlobalMuonZ',
            'GlobalMuonFromCaloMuonZ',
            'TriggerMuonFromGlobalMuonJpsi_JpsiAnal',
            'MuonFromTrackerTrackJpsi_JpsiGlbAnal',
            'MuonFromTrackerTrackJpsi_JpsiTkMAnal',
            'MuonFromTrackerTrackJpsi_JpsiPlusMuAnal',
            'MuonFromTrackerTrackJpsi_BExclAnal'
            ]
    
    curdir = '/tmp/jjhollar'
    files = os.listdir(curdir)

    dbfilename = 'MuonPhysicsPerformance.db'

    i = 0
    for table in tables:
        eff = effs[i]
        metadataprefixname = 'TNPEff_' + str(eff) + '_Table'
        metadatafilename = 'TNPEff_' + str(eff) + '_Table.txt'
        newfile = open(metadatafilename, 'w')
        newfile.write('destDB oracle://cms_orcoff_prep/CMS_COND_PHYSICSTOOLS\n')
        newfile.write('inputtag\n')
        newfile.write('tag ' + table + '\n')
        newfile.write('since\n')
        newfile.write('till\n')
        newfile.write('usertext First test of October exercise Tag&Probe efficiencies ' + eff + ' table\n')
        newfile.close()
        cmd = "uuidgen -t"
        fin,fout = os.popen4(cmd)
        uuid = str(fout.read()).rstrip()

        movedfilename = metadataprefixname + '@' + uuid + '.txt'
        movecommand = 'mv ' + metadatafilename + ' ' + movedfilename
        os.system(movecommand)

        copydbfilename = metadataprefixname + '@' + uuid + '.db'
        copycommand = 'cp ' + dbfilename + ' ' + copydbfilename
        os.system(copycommand)
        print movedfilename
        print copydbfilename
        i = i+1
        
    i = 0
    for wp in wps:
        eff = effs[i]
        metadataprefixname = 'TNPEff_' + str(eff) + '_WP'                 
        metadatafilename = 'TNPEff_' + str(eff) + '_WP.txt'
        newfile = open(metadatafilename, 'w')
        newfile.write('destDB oracle://cms_orcoff_prep/CMS_COND_PHYSICSTOOLS\n')
        newfile.write('inputtag\n')
        newfile.write('tag ' + wp + '\n')
        newfile.write('since\n')
        newfile.write('till\n')
        newfile.write('usertext First test of October exercise Tag&Probe efficiencies ' + eff + ' working point\n')
        newfile.close()                                                                
        cmd = "uuidgen -t"
        fin,fout = os.popen4(cmd)
        uuid = str(fout.read()).rstrip()

        movedfilename = metadataprefixname + '@' + uuid + '.txt'
        movecommand = 'mv ' + metadatafilename + ' ' + movedfilename
        os.system(movecommand)

        copydbfilename = metadataprefixname + '@' + uuid + '.db'
        copycommand = 'cp ' + dbfilename + ' ' + copydbfilename
        os.system(copycommand)
        print movedfilename
        print copydbfilename
        i = i+1
        
if __name__ == "__main__":
    main(sys.argv[1:])
    
