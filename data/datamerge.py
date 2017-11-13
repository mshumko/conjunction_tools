# Combine all of the campaign conjunction data into one file. 

import sys, os, glob
import numpy as np
import spacepy.datamodel
from collections import defaultdict

class ConjunctionWData:
    def __init__(self, sc_a, sc_b, fDir, **kwargs):
        self.sc_a = sc_a
        self.sc_b = sc_b 
        self.fDir = fDir
        self.saveName = kwargs.get('saveName', None)
        self.verbose = kwargs.get('verbose', False)
        return
        
    def combine_conjunction_data(self, fPaths = None):
        """
        This function combines the data from the conjunction
        files located in fDir between two satelies, sc_a and
        sc_b. Output is the file named save_name, 
        in the same format, in the same directory. 
        """
        if fPaths is None:
            self.fPaths = sorted(glob.glob(os.path.join(self.fDir,
                '{}_{}*'.format(self.sc_a, self.sc_b))))
            if self.verbose: print('Path directory:', self.fDir)
            if self.verbose: print('Search string {}_{}*'.format(
                self.sc_a, self.sc_b))
        else:
            self.fPaths = fPaths
        if self.verbose: print('Input daily conjunction paths:', self.fPaths)
        assert len(self.fPaths) > 0, ('ERROR, no merge files found!')
        
        #for i in self.fPaths:
        #    print(os.path.basename(i))
        
        # Read in all the data into a list of dictionaries.
        self.allData = [spacepy.datamodel.readJSONheadedASCII(i) 
            for i in self.fPaths]
        print(self.sc_a, self.sc_b, self.allData[0]['startTime'][0], self.allData[-1]['startTime'][-1])
        
        # Fancy way to combine dictionaries.
        self.outDict = defaultdict(list)
        for d in self.allData:
            for key, value in d.items():
                self.outDict[key] = np.hstack((self.outDict[key], value))       
        return
         
    def writeCombinedConjunctionFile(self, save_dir = None, **kwargs):
        # Now copy/save the data
        self.saveName = kwargs.get('saveName', self.saveName)
        
        self.combinedData = spacepy.datamodel.SpaceData()
        self.combinedData.attrs = self.allData[0].attrs
        
        if save_dir is None:
            save_dir = self.fDir
        keys = list(self.outDict.keys())
        
        for key in keys:
            self.combinedData[key] = spacepy.datamodel.dmarray(
                self.outDict[key], attrs = self.allData[0][key].attrs)
        
        # Rearange the keys to allow special keys with spacecraft IDs.
        # Is there a better way to do this?
        order = keys   
        order.remove('startTime'); order.remove('endTime'); order.remove('duration')
        order.insert(0, 'duration'); order.insert(0, 'endTime'); order.insert(0, 'startTime')
        
        spacepy.datamodel.toJSONheadedASCII(os.path.join(save_dir, 
            self.saveName), self.combinedData, order = order, 
            depend0 = 'startTime')
        return
        
    def removeFiles(self):
        os.remove(self.fPaths)
        pass
        
if __name__ == '__main__':
    for sc_a in ['FU3', 'FU4']:
        for sc_b in ['RBSPA', 'RBSPB']:
            cDir = '/home/mike/research/conjunctiontoolkit/data/daily_conjunctions'
            save_name = '{}_{}_dmlt_1_dL_1_conjunctions.txt'.format(sc_a, sc_b)
            cObj = ConjunctionWData(sc_a, sc_b, cDir, saveName = save_name)
            cObj.combine_conjunction_data()
            cObj.writeCombinedConjunctionFile(save_dir = 
               '/home/mike/research/conjunctiontoolkit/data/merged_conjunctions')
            #cObj.removeFiles()
