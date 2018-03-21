# Proc conjunctions
import glob
import os
import sys
import time
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname( __file__ ), '..'))
import conjunctiontoolkit

dL = 1
dMLT = 1
lowerL = 3
FDIR = '/home/mike/research/conjunction-tools/2018_03_predicted_ephem'

for fb_id in [3, 4]:
    for rb_id in ['A', 'B']:
        FNAME_A = ('FU{}_2018-02-20_2018-03-31_magephem_pre.txt'.format(fb_id))
        FNAME_B = ('rbsp{}_2018-02-25_2018-03-31_magephem_pre.txt'.format(
                    rb_id.lower()))
        SAVE_NAME = ('FU{}_RBSP{}_camp14_conjunctions_pre.txt'.format(
                    fb_id, rb_id))

        cCalc = conjunctiontoolkit.MagneticConjunctionCalc(
            FDIR, FNAME_A, FDIR, FNAME_B,
            mission_id_A='FU{}'.format(fb_id), 
            mission_id_B='RBSP{}'.format(rb_id), verbose=True, 
            Lthresh=dL,  MLTthresh=dMLT, lowerL=lowerL) 
            
        cCalc.magephemB[cCalc.timeKeyB] = np.array([i.replace(tzinfo=None)  
            for i in cCalc.magephemB[cCalc.timeKeyB]])
        
        cCalc.periodic_data_indicies() 
        #cCalc.non_periodic_data_indicies()
        cCalc.calc_magnetic_seperation()
        cCalc.calc_L_MLT_conjunction_delay(0)
        cCalc.lower_L_bound()
        #try: # If no conjunctions are found, log it, and move on.
        cCalc.calc_conjunction_duration()
        # except AssertionError as err:
        #     print(err)
        #     logging.info('{} \n No conjunctions found on day {} \n {} \n'.format(
        #         datetime.now(), np.array(datesA)[i], str(err)))
        #     continue
        #cCalc.plotLMLTandDLDMLT(A_id='FU'+str(fb_id), B_id='RBSP'+str(rb_id))
        cCalc.save_to_file(save_name = SAVE_NAME, save_dir = FDIR)
