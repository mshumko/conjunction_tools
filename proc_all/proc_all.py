# This script processes all of the RBSP-FIREBIRD conjunctions

import os
import sys
import glob

sys.path.insert(0, '/home/mike/research/conjunction-tools/')
import conjunctiontools_v2

baseDir = '/home/mike/research/conjunction-tools/proc_all'
fb_dir = os.path.join(baseDir, 'firebird_camp_magephem')
rbsp_dir = os.path.join(baseDir, 'rbsp_camp_magephem')
mission = 'FIREBIRD'
dL = 1
dMLT = 0.5

for fb_id in [3, 4]:
    for rb_id in ['a', 'b']:
        # Find the mission files.
        fbFiles = sorted(glob.glob(
                    os.path.join(fb_dir, 'FU{}_camp*'.format(fb_id))
                    ))
        rbFiles = sorted(glob.glob(
                    os.path.join(rbsp_dir, 'rbsp{}_camp*'.format(rb_id))
                    ))

        for fb_f, rb_f in zip(fbFiles, rbFiles):
            print('Processing files {} and {}'.format(
                os.path.basename(fb_f), (os.path.basename(rb_f))))

            m = conjunctiontools_v2.MagneticConjunctions(
                mission, mission, fb_f, rb_f, Lthresh=dL, MLTthresh=dMLT)
            m.calcConjunctions()
            m.saveData(os.path.join(baseDir, 'conjunctions', 
                'FU{}_RBSP{}_conjunctions_dL{}_dMLT{}.txt'.format(
                    fb_id, rb_id.upper(), int(10*dL), int(10*dMLT) )
                ))