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
FDIR = '/home/mike/research/conjunction-tools/2017_11_predicted_ephem'

for fb_id in [3, 4]:
    for rb_id in ['A', 'B']:
        FNAME_A = ('FU{}_SGP4_LLA_2017-11-13_to_2017-12-22_gen_w_'
            '2017-11-13_TLE_magephem.txt'.format(fb_id))
        FNAME_B = ('RBSP{}_2017-11-13_2017-12-23_magephem.'
            'txt'.format(rb_id))
        SAVE_NAME = ('2017_nov_predicted_FU{}_RBSP{}'
            '_conjunctions.txt'.format(fb_id, rb_id))

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
        #cCalc.plotLMLTandDLDMLT( 'FU' + str(fb_id), 'RBSP' + str(rbsp_id))
        cCalc.save_to_file(save_name = SAVE_NAME, save_dir = FDIR)