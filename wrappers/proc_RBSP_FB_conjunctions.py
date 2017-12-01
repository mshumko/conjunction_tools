import glob, os, sys
import time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname( __file__ ), '..'))
sys.path.insert(0, '/home/mike/research/mission-tools/misc')
import conjunctiontoolkit
import dates_in_filenames
from datetime import datetime
import re # Regular experssion matching for erro catching
import logging
import argparse

# Begin logging and keeping track of start time.
logging.basicConfig(filename='proc_RBSP_FB_conjunctions.log', level=logging.DEBUG)
startTime = time.time()

# Parse optional arguments given to the script
parser = argparse.ArgumentParser()
parser.add_argument("-dL", type=float, default = 1, nargs = '?',
                    help="delta L conjunction threshold")
parser.add_argument("-dMLT", type=float, default = 1, nargs = '?',
                    help="delta MLT conjunction threshold")
parser.add_argument("-lowerL", type=float, default = 3, nargs = '?',
                    help="Lower L conjunction threshold")                    
args = parser.parse_args()

print(args)

dL = args.dL
dMLT = args.dMLT
lowerL = args.lowerL
years = [2017]

for fb_id in ['3', '4']:
    for rbsp_id in ['A', 'B']:
        fDirA = '/home/mike/research/firebird/Datafiles/FU_{}/magephem'.format(
            fb_id)
        fDirB =os.path.join('/home/mike/research/rbsp/magephem/rbsp{}'.format(
            rbsp_id.lower()))

        # Now inteligently combine the txt files.
        pathsA, datesA = dates_in_filenames.find_dates(fDirA, 
            wildCard='*FU{}_T89*.txt'.format(fb_id))
        pathsB, datesB = dates_in_filenames.find_dates(fDirB, 
            wildCard='*T89D*.txt')
    
        for i in range(len(datesA)):
            # Ignore non daily magnetic ephemeris files.
            if len(datesA[i]) > 1: continue
          
            match = np.where(np.array(datesA)[i] == np.array(datesB))[0]
            if len(match) == 0:
                print('No matches found for {}'.format(datesA[i]))
                continue
            elif len(match) == 1:
                # Run the conjunction code here! 
                fnameA = os.path.basename(pathsA[i])
                fnameB = os.path.basename(pathsB[match[0]])
                # Do not analyze years outside of the years given by years array.
                if not int(fnameA[0:4]) in years: 
                    continue
                print('Calculating conjunctions on files:')
                print(fnameA, fnameB)
                
                # MODIFY THIS!
                saveName = ('FU{}_RBSP{}_'.format(fb_id, rbsp_id) + 
                    datesA[i][0] + '_conjunctions.txt')
                save_dir = ('/home/mike/research/conjunction-tools/data/'
                    'daily_conjunctions')
                    
                # Select the pitch angle of Lstar to use
                LcolB = -1
                    
                cCalc = conjunctiontoolkit.MagneticConjunctionCalc(
                    fDirA, fnameA, fDirB, fnameB, mission_id_A='FU'+str(fb_id), 
                    mission_id_B='RBSP'+str(rbsp_id), verbose=True, 
                    Lthresh=dL,  MLTthresh=dMLT, timeKeyB='DateTime', 
                    MLTkeyB='EDMAG_MLT', LkeyB='L', LcolB=LcolB, 
                    lowerL=lowerL) 
                    
                cCalc.magephemB[cCalc.timeKeyB] = np.array([i.replace(tzinfo=None)  
                    for i in cCalc.magephemB[cCalc.timeKeyB]])
              
                cCalc.periodic_data_indicies() 
                #cCalc.non_periodic_data_indicies()
                cCalc.calc_magnetic_seperation()
                cCalc.calc_L_MLT_conjunction_delay(0)
                cCalc.lower_L_bound()
                try: # If no conjunctions are found, log it, and move on.
                    cCalc.calc_conjunction_duration()
                except AssertionError as err:
                    print(err)
                    logging.info('{} \n No conjunctions found on day {} \n {} \n'.format(
                        datetime.now(), np.array(datesA)[i], str(err)))
                    continue
                #cCalc.plotLMLTandDLDMLT( 'FU' + str(fb_id), 'RBSP' + str(rbsp_id))
                cCalc.save_to_file(save_dir=save_dir)
                
print('Program execution:{}'.format(time.time() - startTime))
