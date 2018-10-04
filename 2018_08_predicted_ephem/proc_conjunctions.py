# Proc conjunctions
import glob
import os
import sys
import time
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname( __file__ ), '..'))
import conjunctiontools

baseDir = '/home/mike/research/conjunction-tools/proc_all'
saveDir = '/home/mike/research/conjunction-tools/2018_08_predicted'
fb_dir = os.path.join(baseDir, 'firebird_camp_magephem')
rbsp_dir = os.path.join(baseDir, 'rbsp_camp_magephem')
mission = 'FIREBIRD'
dL = 1
dMLT = 1

for fb_id in [3, 4]:
    for rb_id in ['A', 'B']:
        # Find the mission files.
        fbFiles = sorted(glob.glob(
                    os.path.join(fb_dir, 'FU{}_camp17*'.format(fb_id))
                    ))
        rbFiles = sorted(glob.glob(
                    os.path.join(rbsp_dir, 'rbsp{}_camp17*'.format(rb_id.lower()))
                    ))
        SAVENAME = 'FU{}_RBSP{}_conjunctions_dL{}_dMLT{}_old.txt'.format(
                    fb_id, rb_id.upper(), int(10*dL), int(10*dMLT))
        print(fbFiles, rbFiles)
        FNAME_A = os.path.basename(fbFiles[0])
        FNAME_B = os.path.basename(rbFiles[0])
        FDIR_A = os.path.dirname(fbFiles[0])
        FDIR_B = os.path.dirname(rbFiles[0])

        cCalc = conjunctiontools.MagneticConjunctionCalc(
            FDIR_A, FNAME_A, FDIR_B, FNAME_B,
            mission_id_A='FU{}'.format(fb_id), 
            mission_id_B='RBSP{}'.format(rb_id), verbose=True, 
            Lthresh=dL,  MLTthresh=dMLT, lowerL=3) 
            
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
        cCalc.save_to_file(save_name=SAVENAME, save_dir=saveDir)
