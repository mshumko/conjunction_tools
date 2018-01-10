# This script will create conjunction list 

import numpy as np
import sys
import os
import itertools
import glob
import fnmatch
from datetime import datetime, timedelta

# import various personal libraries to run the code.
sys.path.insert(0, '/home/mike/research/mission-tools/rbsp')
import download_RBSP_magephem
sys.path.insert(0, '/home/mike/research/firebird/data_processing/magnetic_ephemeris')
import make_magnetic_ephemeris
sys.path.insert(0, '/home/mike/research/conjunction-tools/')
import conjunctiontoolkit
sys.path.insert(0, '/home/mike/research/conjunction-tools/data/')
import datamerge

# Set user parameters for the script
START_DATE = datetime(2017, 11, 18)
END_DATE = datetime(2017, 12, 14)
FIREBIRD_KP = 20

firebird_ephem_dir = lambda sc_id: '/home/mike/research/firebird/Datafiles/FU_{}/ephem/'.format(sc_id)
firebird_magepehem_dir = lambda sc_id: '/home/mike/research/firebird/Datafiles/FU_{}/magephem/'.format(sc_id)
rbsp_magephem_dir = lambda sc_id: ('/home/mike/research/rbsp/magephem/rbsp{}/'.format(sc_id.lower()))
DAILY_SAVE_DIR = ('/home/mike/research/conjunction-tools/data/daily_conjunctions')
MERGED_SAVE_DIR = ('/home/mike/research/conjunction-tools/data/merged_conjunctions')


def download_rbsp_magephem_wrapper():
    """Download RBSP magepehem files for the days of interest"""
    dates = [START_DATE + timedelta(days=i) for i in range((END_DATE-START_DATE).days)]
    for (sc, date) in itertools.product(['a', 'b'], dates):
        print('Downloading RBSP-{} MagEphem from {}'.format(sc, date))
        download_RBSP_magephem.saveRBSPMagEphem(sc, date, rbsp_magephem_dir(sc))

def generate_firebird_magepehem_wrapper():
    """ Process FIREBIRD magephem files"""
    for sc_id in ['3', '4']:
        # Find the correct ephem file to process
        ephemPaths = glob.glob(firebird_ephem_dir(sc_id) + '*{}*{}*'.format(
            START_DATE.date(), END_DATE.date()))
        assert len(ephemPaths) == 1, 'None or multiple FIREBIRD ephemeris file found.'
        magEphem = make_magnetic_ephemeris.Magnetic_Ephemeris('FB', sc_id = sc_id,
            ephemName=ephemPaths[0], ephemDir=firebird_ephem_dir(sc_id), 
            singleKp=FIREBIRD_KP)
        magEphem.run_mag_model(multiprocess=True, nonPeriodicFlag=False,
                                mirrorPointAlt=False)
        magEphem.saveToFile(daily=True)

### GET DAILY CONJUNCTION LISTS ###
class ConjunctionWrapper:
    def __init__(self, fb_id, rbsp_id, startDate, endDate):
        """ 
        This is a wrapper class to process magnetic ephemeris 
        for the FIREBIRD-RBSP conjunctions.
        """
        self.fDirA = firebird_magepehem_dir(fb_id)
        self.fDirB = rbsp_magephem_dir(rbsp_id)
        self.fb_id = fb_id
        self.rbsp_id = rbsp_id
        self.dates = [startDate + timedelta(days=i) for i in
            range( (endDate - startDate).days )]
        return

    def matchFiles(self):
        """
        This method will go through the requested dates to calculate the
        conjunction, and search for magephem files for each spacecraft.
        """
        self.matchedMagFile = np.zeros((len(self.dates), 3), dtype=object)
        filesA = glob.glob(self.fDirA + '*')
        filesB = glob.glob(self.fDirB + '*')
        
        for (i, date) in enumerate(self.dates):
            # Find the magephem dates from both spacecraft.
            dateStr = '*'+date.date().isoformat().replace('-', '')+'*'
            matchesA = fnmatch.filter(filesA, dateStr)
            matchesB = fnmatch.filter(filesB, dateStr)
            assert len(matchesA) == 1 and len(matchesB) == 1, (
                'None or multiple magephem files found on {}'.format(date.date()))
            self.matchedMagFile[i, :] = [dateStr, matchesA[0], matchesB[0]]
        return

    def procConjunctions(self, LcolB=-1, dL=1, dMLT=1, lowerL=3):
        """ 
        This method will run the conjunction calculator for each day.
        Lcol is an integer on which L shell column to use from RBSP 
        magephem.
        """
        for matches in self.matchedMagFile:
            saveName = ('FU{}_RBSP{}_'.format(self.fb_id, self.rbsp_id) + 
                matches[0] + '_conjunctions.txt')

            fnameA = os.path.basename(matches[1])
            fnameB = os.path.basename(matches[2])

            cCalc = conjunctiontoolkit.MagneticConjunctionCalc(
                self.fDirA, fnameA, self.fDirB, fnameB, 
                mission_id_A=self.fb_id, mission_id_B=self.rbsp_id, 
                verbose=True, Lthresh=dL,  MLTthresh=dMLT, 
                timeKeyB='DateTime', MLTkeyB='EDMAG_MLT', LkeyB='L', 
                LcolB=LcolB, lowerL=lowerL)
                
            cCalc.magephemB[cCalc.timeKeyB] = np.array([i.replace(tzinfo=None)  
                for i in cCalc.magephemB[cCalc.timeKeyB]])
            
            cCalc.periodic_data_indicies() 
            cCalc.calc_magnetic_seperation()
            cCalc.calc_L_MLT_conjunction_delay(0)
            cCalc.lower_L_bound()
            try: # If no conjunctions are found, log it, and move on.
                cCalc.calc_conjunction_duration()
            except AssertionError as err:
                print(err)
                continue
            cCalc.save_to_file(save_dir=DAILY_SAVE_DIR)

    def mergeFiles(self):
        """ This is a wrapper method to merge conjunction files """
        save_name = '{}_{}_dmlt_1_dL_1_conjunctions.txt'.format(
            self.fb_id, self.rbsp_id)
        cObj = datamerge.ConjunctionWData(self.fb_id, self.rbsp_id, 
            DAILY_SAVE_DIR, saveName=save_name)
        cObj.combine_conjunction_data()
        cObj.writeCombinedConjunctionFile(save_dir=MERGED_SAVE_DIR)

### RUN THE SCRIPT ###
if __name__ == '__main__':
    #download_rbsp_magephem_wrapper()
    #generate_firebird_magepehem_wrapper()

    for (fb_id, rbsp_id) in itertools.product(['3', '4'], ['A', 'B']):
        c = ConjunctionWrapper(fb_id, rbsp_id, START_DATE, END_DATE)
        c.matchFiles()
        c.procConjunctions()
        c.mergeFiles()
# for fb_id in ['3', '4']:
#     for rbsp_id in ['A', 'B']:
#         fDirA = firebird_magepehem_dir(fb_id)
#         fDirB = directories.rbsp_magephem_dir(fbsp_id)

#         # Now inteligently combine the txt files.
#         pathsA, datesA = dates_in_filenames.find_dates(fDirA, 
#             wildCard='*FU{}_T89*.txt'.format(fb_id))
#         pathsB, datesB = dates_in_filenames.find_dates(fDirB, 
#             wildCard='*T89D*.txt')
    
#         for i in range(len(datesA)):
#             # Ignore non daily magnetic ephemeris files.
#             if len(datesA[i]) > 1: continue
          
#             match = np.where(np.array(datesA)[i] == np.array(datesB))[0]
#             if len(match) == 0:
#                 print('No matches found for {}'.format(datesA[i]))
#                 continue
#             elif len(match) == 1:
#                 # Run the conjunction code here! 
#                 fnameA = os.path.basename(pathsA[i])
#                 fnameB = os.path.basename(pathsB[match[0]])
#                 # Do not analyze years outside of the years given by years array.
#                 if not int(fnameA[0:4]) in years: 
#                     continue
#                 print('Calculating conjunctions on files:')
#                 print(fnameA, fnameB)
                
#                 # MODIFY THIS!
#                 saveName = ('FU{}_RBSP{}_'.format(fb_id, rbsp_id) + 
#                     datesA[i][0] + '_conjunctions.txt')
#                 save_dir = ('/home/mike/research/conjunction-tools/data/'
#                     'daily_conjunctions')
                    
#                 # Select the pitch angle of Lstar to use
#                 LcolB = -1
                    
#                 cCalc = conjunctiontoolkit.MagneticConjunctionCalc(
#                     fDirA, fnameA, fDirB, fnameB, mission_id_A='FU'+str(fb_id), 
#                     mission_id_B='RBSP'+str(rbsp_id), verbose=True, 
#                     Lthresh=dL,  MLTthresh=dMLT, timeKeyB='DateTime', 
#                     MLTkeyB='EDMAG_MLT', LkeyB='L', LcolB=LcolB, 
#                     lowerL=lowerL) 
                    
#                 cCalc.magephemB[cCalc.timeKeyB] = np.array([i.replace(tzinfo=None)  
#                     for i in cCalc.magephemB[cCalc.timeKeyB]])
              
#                 cCalc.periodic_data_indicies() 
#                 #cCalc.non_periodic_data_indicies()
#                 cCalc.calc_magnetic_seperation()
#                 cCalc.calc_L_MLT_conjunction_delay(0)
#                 cCalc.lower_L_bound()
#                 try: # If no conjunctions are found, log it, and move on.
#                     cCalc.calc_conjunction_duration()
#                 except AssertionError as err:
#                     print(err)
#                     logging.info('{} \n No conjunctions found on day {} \n {} \n'.format(
#                         datetime.now(), np.array(datesA)[i], str(err)))
#                     continue
#                 #cCalc.plotLMLTandDLDMLT( 'FU' + str(fb_id), 'RBSP' + str(rbsp_id))
#                 cCalc.save_to_file(save_dir=save_dir)
