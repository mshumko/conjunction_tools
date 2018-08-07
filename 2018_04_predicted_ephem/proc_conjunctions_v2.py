# Proc conjunctions
import glob
import os
import sys
import time
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname( __file__ ), '..'))
import conjunctiontools_v2

dL = 1
dMLT = 1
lowerL = 3
FDIR = '/home/mike/research/conjunction-tools/2018_04_predicted_ephem'
missionA = 'FIREBIRD' # For data formatting keys.
missionB = 'FIREBIRD'

for fb_id in [3, 4]:
    for rb_id in ['A', 'B']:
        FNAME_A = ('FU{}_2018-04-18_2018-05-31_magephem_pre.txt'.format(fb_id))
        FNAME_B = ('rbsp{}_2018-04-18_2018-05-31_magephem_pre.txt'.format(
                    rb_id.lower()))
        SAVE_NAME = ('FU{}_RBSP{}_camp15_conjunctions_pre_v2.txt'.format(
                    fb_id, rb_id))
        magA = os.path.join(FDIR, FNAME_A)
        magB = os.path.join(FDIR, FNAME_B)
        m = conjunctiontools_v2.MagneticConjunctions(missionA, missionB, magA, magB)
        m.calcConjunctions()
        m.saveData(os.path.join(FDIR, SAVE_NAME))
        #m.testPlots()
        #plt.show()
