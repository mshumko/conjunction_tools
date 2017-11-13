import glob, os, sys
import time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname( __file__ ), '..'))
sys.path.insert(0, os.path.join(os.path.dirname( __file__ ), '..', '..', 'auxiliary'))
import conjunctiontoolkit
import dates_in_filenames
from datetime import datetime
import re # Regular experssion matching for erro catching
import logging
import argparse
import read_ac_data

logging = False

# Begin logging and keeping track of start time.
if logging:
    logging.basicConfig(filename='proc_RBSP_AC_conjunctions.log',
    level=logging.DEBUG)
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

dL = args.dL
dMLT = args.dMLT
lowerL = args.lowerL

print('Running conjunction calculator with dL = {}, dMLT = {}, lowerL = {}'.format(
    dL, dMLT, lowerL))

for ac_id in ['A', 'B']:
    for rbsp_id in ['A', 'B']:
        # Define data paths.
        fDirA = '/home/ms30715/ssd_data/ac6/ac6{}/ascii/level2'.format(ac_id.lower())
        # Since RBSP magEphem was interpolated onto AC6 time stamps, 
        fDirB =os.path.join(os.sep, 'home', 'ms30715', 'ssd_data', 'shumko', 
            'data', 'RBSP', 'orbit', 'RBSP{}'.format(rbsp_id.upper()),
            'magephem_interpolated', 'AC6{}'.format(ac_id.upper()))

        # Now inteligently combine the txt files.
        pathsA, datesA = dates_in_filenames.find_dates(fDirA, 
            wildCard =  '*coords*.csv')
        pathsB, datesB = dates_in_filenames.find_dates(fDirB, 
            wildCard =  '*TS04D*AC6{}*.txt'.format(ac_id.upper()))

        for i in range(len(datesA)):
            # Ignore non daily magnetic ephemeris files.
            if len(datesA[i]) > 1: continue
                      
            match = np.where(np.array(datesA)[i] == np.array(datesB))[0]
            if len(match) == 0:
                #print('No matches found for {}'.format(datesA[i]))
                continue
            elif len(match) == 1:
                # Run the conjunction code here! 
                fnameA = os.path.basename(pathsA[i])
                fnameB = os.path.basename(pathsB[match[0]])
                #print('Calculating conjunctions on files:')
                print(fnameA, fnameB)

                ### Read in the AC6 data. ###
                magephemAC6 = read_ac_data.read_ac_data(pathsA[i])
                
                ### SAVE DIRECTORY AND FILENAME ###
                saveName = ('AC{}_RBSP{}_'.format(ac_id, rbsp_id) + 
                    datesA[i][0] + '_conjunctions.txt')
                save_dir = ('/home/ms30715/ssd_data/shumko/'
                    '0_code/conjunctiontoolkit/data/daily_conjunctions' )
                    
                # Select the pitch angle of Lstar to use
                LcolB = -2
                    
                cCalc = conjunctiontoolkit.MagneticConjunctionCalc(
                    fDirA, fnameA, fDirB, fnameB, mission_id_A = 'AC6' + str(ac_id), 
                    mission_id_B = 'RBSP' + str(rbsp_id), verbose = False, 
                    Lthresh = dL,  MLTthresh = dMLT, timeKeyA = 'dateTime', 
                    MLTkeyA = 'MLT_OPQ', LkeyA = 'Lm_OPQ', timeKeyB = 'DateTime', 
                    MLTkeyB = 'EDMAG_MLT', LkeyB = 'Lstar', LcolB = LcolB, 
                    lowerL = lowerL, magEphemA = magephemAC6) 
                cCalc.cadence = 1/60
                
                # Remove time zone info    
#                cCalc.magephemB[cCalc.timeKeyB] = np.array([i.replace(tzinfo=None)  
#                    for i in cCalc.magephemB[cCalc.timeKeyB]])

                #cCalc.bruteForceConjunctionIndicies()
              
##                cCalc.match_start_end_ind() 
                cCalc.calc_magnetic_seperation()
                cCalc.calc_L_MLT_conjunction_delay(0)
                cCalc.lower_L_bound()
                try: # If no conjunctions are found, log it, and move on.
                    cCalc.calc_conjunction_duration()
                except AssertionError as err:
                    if logging:
                        logging.info('{} \n No conjunctions found on day {}' 
                        '\n {} \n'.format(datetime.now(), np.array(datesA)[i], 
                        str(err)))
                    continue 
##                cCalc.plotLMLTandDLDMLT(A_id = 'AC6' + str(ac_id), 
##                    B_id = 'RBSP' + str(rbsp_id))
                cCalc.save_to_file(save_name = saveName, save_dir = save_dir)
print('Program execution:{}'.format(time.time() - startTime))
