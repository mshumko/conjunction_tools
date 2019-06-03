# This script processes all of the FIREBIRD conjunctions with itself.

import os
import sys
import glob
import itertools

sys.path.insert(0, '/home/mike/research/conjunction-tools/')
import conjunctiontools_v2

baseDir = '/home/mike/research/conjunction-tools/proc_all'
fb_dir = os.path.join(baseDir, 'firebird_camp_magephem')
mission = 'FIREBIRD'
camp = 22 # If process all campaigns, camp = ''

dL = 1
dMLT = 2

# Find the mission files.
fu3Files = sorted(glob.glob(
            os.path.join(fb_dir, 'FU3_camp{}*'.format(camp))
            ))
fu4Files = sorted(glob.glob(
            os.path.join(fb_dir, 'FU4_camp{}*'.format(camp))
            ))

for fu3_f, fu4_f in zip(fu3Files, fu4Files):
    print('Processing files {} and {}'.format(
        os.path.basename(fu3_f), (os.path.basename(fu4_f))))

    m = conjunctiontools_v2.MagneticConjunctions(
        mission, mission, fu3_f, fu4_f, Lthresh=dL, MLTthresh=dMLT)
    m.calcConjunctions()
    
    for L, dLi, dMLTi in zip(m.magA['L'], m.dL, m.dMLT):
        if dLi < 1 and dMLTi < 2:
            print(L, dLi, dMLTi)
    
    if camp == '':
        m.saveData(os.path.join(baseDir, 'conjunctions', 
            'FU3_FU4_conjunctions_dL{}_dMLT{}.txt'.format(
                int(10*dL), int(10*dMLT) )
            ), mode='a')
    else:
        m.saveData(os.path.join(baseDir, 'conjunctions', 
        'FU3_FU4_camp{}_conjunctions_dL{}_dMLT{}.txt'.format(
            camp, int(10*dL), int(10*dMLT))
            ), mode='a')
