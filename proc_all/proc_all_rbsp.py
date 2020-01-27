# This script processes all of the RBSP-FIREBIRD conjunctions

import os
import sys
import glob
import itertools

sys.path.insert(0, '/home/mike/research/conjunction-tools/')
import conjunctiontools_v2

baseDir = '/home/mike/research/conjunction-tools/proc_all'
fb_dir = os.path.join(baseDir, 'firebird_camp_magephem')
rbsp_dir = os.path.join(baseDir, 'rbsp_camp_magephem')
mission = 'FIREBIRD'
camp = '' # If process all campaigns, camp = ''

dLArr = [1]
dMLTArr = [1]

for dL, dMLT in zip(dLArr, dMLTArr):
    for fb_id, rb_id in itertools.product([3, 4], ['rbspa', 'rbspb']): # 'arase', rbspa', 'rbspb', 
        # Find the mission files.
        fbFiles = sorted(glob.glob(
                    os.path.join(fb_dir, 'FU{}_camp{}*'.format(fb_id, camp))
                    ))
        rbFiles = sorted(glob.glob(
                    os.path.join(rbsp_dir, '{}_camp{}*'.format(rb_id, camp))
                    ))

        for fb_f, rb_f in zip(fbFiles, rbFiles):
            print('Processing files {} and {}'.format(
                os.path.basename(fb_f), (os.path.basename(rb_f))))

            m = conjunctiontools_v2.MagneticConjunctions(
                mission, mission, fb_f, rb_f, Lthresh=dL, MLTthresh=dMLT)
            m.calcConjunctions()
            if camp == '':
                m.saveData(os.path.join(baseDir, 'conjunctions', 
                    'FU{}_{}_conjunctions_dL{}_dMLT{}.txt'.format(
                        fb_id, rb_id.upper(), int(10*dL), int(10*dMLT) )
                    ), mode='a')
            else:
                m.saveData(os.path.join(baseDir, 'conjunctions', 
                'FU{}_{}_camp{}_conjunctions_dL{}_dMLT{}.txt'.format(
                    fb_id, rb_id.upper(), camp, int(10*dL), int(10*dMLT))
                    ), mode='a')
