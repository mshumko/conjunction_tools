# This script processes all of the RBSP-FIREBIRD conjunctions

import os
import sys
import glob
import itertools
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/mike/research/conjunction-tools/')
import conjunctiontools_v2

baseDir = '/home/mike/research/conjunction-tools/proc_all'
fb_dir = os.path.join(baseDir, 'firebird_camp_magephem')
dsx_dir = os.path.join(baseDir, 'dsx_camp_magephem')
#mission = 'FIREBIRD'
camp = 26 # If process all campaigns, camp = ''

dLArr = [1]
dMLTArr = [1]

debug = False

for dL, dMLT in zip(dLArr, dMLTArr):
    for fb_id in [3, 4]: # 'arase', rbspa', 'rbspb', 
        # Find the mission files.
        fbFiles = sorted(glob.glob(
                    os.path.join(fb_dir, 'FU{}_camp{}*'.format(fb_id, camp))
                    ))
        dsxFiles = sorted(glob.glob(
                    os.path.join(dsx_dir, 'DSX_camp{}*'.format(camp))
                    ))

        for fb_f, rb_f in zip(fbFiles, dsxFiles):
            print('Processing files {} and {}'.format(
                os.path.basename(fb_f), (os.path.basename(rb_f))))

            m = conjunctiontools_v2.MagneticConjunctions(
                'FIREBIRD', 'DSX', fb_f, rb_f, Lthresh=dL, MLTthresh=dMLT)
            m.calcConjunctions()
            if debug:
                m.plot_overview()
            else:
                if camp == '':
                    m.saveData(os.path.join(baseDir, 'conjunctions', 
                        'FU{}_DSX_conjunctions_dL{}_dMLT{}.txt'.format(
                            fb_id, int(10*dL), int(10*dMLT) )
                        ), mode='a')
                else:
                    m.saveData(os.path.join(baseDir, 'conjunctions', 
                    'FU{}_DSX_camp{}_conjunctions_dL{}_dMLT{}.txt'.format(
                        fb_id, camp, int(10*dL), int(10*dMLT))
                        ), mode='a')
