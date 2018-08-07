import numpy as np
import csv
import sys
import os
import glob
from datetime import datetime
import dateutil.parser
import re

import spacepy.datamodel

sys.path.append('/home/mike/research/mission-tools/misc')
import dates_in_filenames

CONJUNCTION_DIR = ('/home/mike/research/conjunction-tools'
                    '/proc_all/conjunctions')
HIRES_DIR = '/home/mike/research/firebird/Datafiles'

class ConjunctionHRfilter:
    def __init__(self, cPath, hrDir):
        fbRegex = re.compile(r'FU3|FU4')
        mo = fbRegex.search(cPath)
        self.fb_id = mo.group()
        self.hrDir = os.path.join(hrDir, 
                'FU_{}/hires/level2'.format(self.fb_id[-1]))
        
        # Load conjunction file
        self._load_conjunction_file(cPath)
        # Get list and dates of all HiRes files
        self._get_hr_paths()
        return
    
    def filterConjunctions(self):
        """
        This method filters the conjunction dictionary.
        """
        self.fltConjunctions = {}
        for key in self.cKeys:
            self.fltConjunctions[key] = np.array([])
            
        # Loop over the conjunctions
        d0 = datetime.min
        z = zip(self.cData['startTime'], self.cData['endTime'])
        for i, (sT, eT) in enumerate(z):
            # Load HiRes file if the conjucntion date 
            # does not match with the currently loaded 
            # HiRes file (d0)
            if sT.date() != d0.date():
                # Load HiRes
                flag = self._load_hires(sT)
                if flag == -1: # If no HiRes file was found on that day.
                    continue
                d0 = sT
                print('Processing date:', d0.date())
                
            # Check if the HiRes data exists for that time.
            idhr = np.where((self.hrTimes > sT) &
                            (self.hrTimes < eT))[0]
            if len(idhr) > 0:
                for key in self.cKeys:
                    self.fltConjunctions[key] = np.append(
                        self.fltConjunctions[key], 
                        self.cData[key][i])
                
        return
        
    def _get_hr_paths(self):
        """
        This method searched the FIREBIRD HiRes directory
        and finds all of the filenames, and saves their 
        dates.
        """
        self.hrPaths, self.hrDates = dates_in_filenames.find_dates(
            self.hrDir, dateTimeConvert=True)
        self.hrDates = np.array([t[0].date()
                        for t in self.hrDates])
        return
        
    def _load_hires(self, time):
        """
        This methods loads in a HiRes file from the 
        predefined fb_id, and time (used to extract date).
        """
        idt = np.where(self.hrDates == time.date())[0]
        if len(idt) != 1:
            return -1
        cPath = self.hrPaths[idt[0]]
        hr = spacepy.datamodel.readJSONheadedASCII(cPath)
        # Convert times
        self.hrTimes = np.array([dateutil.parser.parse(t) 
                for t in hr['Time']])
        return 1
                
    def _load_conjunction_file(self, cPath):
        """
        This method reads in the csv conjunction file
        and saves it into a self.cData dictionary.
        """
        self.cData = {}
        with open(cPath, 'r') as f:
            r = csv.reader(f)
            self.cKeys = next(r)
            data = np.array(list(r))
        # Now parse the 2d array and save it to 
        # self.cData. Convert the start/end times
        validInd = np.where(data[:, 0] != 'nan')[0]
        
        for i, key in enumerate(self.cKeys):
            if 'time' in key.lower():
               self.cData[key] = np.array(
                    [dateutil.parser.parse(t) 
                    for t in data[validInd, i]])
            else:
                self.cData[key] = data[validInd, i]
        return

# Load the conjunction files one by one
cPaths = glob.glob(os.path.join(CONJUNCTION_DIR, '*.txt'))

for f in cPaths:
    c = ConjunctionHRfilter(f, HIRES_DIR)
    c.filterConjunctions()
    
    
        
