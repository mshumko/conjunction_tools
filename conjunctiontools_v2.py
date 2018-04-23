import time
from datetime import datetime, timedelta
import spacepy.datamodel as dm
import dateutil.parser
import multiprocessing
import numpy as np
import os, sys
import matplotlib.pylab as plt
from matplotlib.dates import date2num
import operator 
import itertools
import scipy

import IRBEM

# This dictionary contains the magEphem L, MLT, lat, lon, alt, 
# time keys for reading in the data. 
keyDict = {
        'FIREBIRD':{
            'L':'McllwainL', 'MLT':'MLT', 'time':'dateTime', 
            'lat':'Rgeod_LatLon', 'latcol':0, 'lon':'Rgeod_LatLon', 
            'loncol':1, 'alt':'Rgeod_Altitude'
            },

        'AC6':{
            'L':'Lm_OPQ', 'MLT':'MLT_OPQ', 'time':'dateTime',
            'lat':None, 'lon':None, 'alt':None
            },

        'RBSP':{
            'L':'Lstar', 'Lcol':-1, 'MLT':'EDMAG_MLT', 
            'time':'DateTime', 'lat':'Rgeod_LatLon', 'latcol':0,
            'lon':'Rgeod_LatLon', 'loncol':1, 'alt':'Rgeod_Height'
            }
        }

class MagneticConjunctions(IRBEM.MagFields):
    def __init__(self, missionA, missionB, magA, magB, **kwargs):
        """
        NAME:    MagneticConjunctions()
        USE:     This is the second iteration of the magnetic conjunction calc.
                 It easily combines the two magnetic ephemeris files, and cals.
                 the nearest footprint separation with an interpolation of L
                 and MLT.

                 If fPathA/B is a tuple, the values must be in order and contain:
                 file path to a JSON headed ASCII file, time key, MLT key, L key,
                 and optional L shell column number.
        RETURNS: None 
        AUTHOR:  Mykhaylo Shumko
        MOD:     2018-4-8
        """
        # Basic conjunction params
        self.lowerL = kwargs.get('lowerL', 1)
        self.Lthresh = kwargs.get('Lthresh', 1)
        self.MLTthresh = kwargs.get('MLTthresh', 1)

        # magephem args
        #self.magArgsA = magephemAargs
        #self.magArgsB = magephemBargs
        
        self.magArgsA = kwargs.get('magArgsA', keyDict[missionA])
        self.magArgsB = kwargs.get('magArgsB', keyDict[missionB])

        # Aux params
        self.REPLACE_ERROR_VALS = kwargs.get('REPLACE_ERROR_VALS', np.nan)
        self.Kp = kwargs.get('Kp', 2) # For interplation.

        # Initializa IRBEM
        IRBEM.MagFields.__init__(self)

        ### Load data if necessary (assuming it is in a JSON headed ASCII file) ###
        if isinstance(magA, str):
            self.magA = self._load_magephem(magA, self.magArgsA)
        else:
            self.magA = magA
            
        if isinstance(magB, str):
            self.magB = self._load_magephem(magB, self.magArgsB)
        else:
            self.magB = magB
        self._find_common_times()
        return
    
    def calc_conjunctions(self, interpP=None):
        """
        This method calculates the conjunctions, and their duration.
        If interpP is a tuple with lat, lon, alt keys for spacecraft 
        A and B, as well as a footpoint alt (total of 7 params) in 
        that order, then for each conjunction, the lat/lon/alt will
        be interpolated, and this function will use OPQ model to 
        return dMLT at min(dL), as well as the smallest footpoint 
        separation at a given altitude (in whatever hemisphere it 
        is the smallest).
        """
        self.dL = np.abs(self.magA['L'] - self.magB['L'])
        self.dMLT = dmlt(self.magA['MLT'], self.magB['MLT'])

        # Calc indicies where separation meets conjunction criteria.
        idC = np.where((self.dL < self.Lthresh) & (self.dMLT < self.MLTthresh))[0]
        #idC = idC[:-1]
        #print(idC)
        self.startInd, self.endInd = self._calc_start_stop(idC)



        ### TEST CODE ###
        # for sI, eI in zip(startInd, endInd):
        #     plt.plot(self.dL[sI-2:eI+1], 'r')
        #     interpT = np.linspace(0, (eI+2-sI+1))
        #     fL = scipy.interpolate.interp1d(np.arange(eI+2-sI+2), self.dL[sI-2:eI+2], kind='quadratic')
        #     plt.scatter(interpT, fL(interpT), c='r')
        #     #plt.plot(self.dMLT[sI-1:eI+1])
        # plt.show()
        #for sI, eI in zip(startInd, endInd):
        #    print(self.dL[sI:eI], self.dMLT[sI:eI])
        #######
        # Now interpolate L and MLT around the flagged conjunctions.

        ### Calculate closest footpoint separation ###
        return

    def testPlots(self):
        """ 
        This method plots the dL and dMLT for each of the times 
        flagged for a conjunction.
        """
        for si, ei in zip(self.startInd, self.endInd):
            fig, ax = plt.subplots()
            ax.plot(self.magA['dateTime'][si-5:ei+5], self.dL[si-5:ei+5], 'r', label='dL')
            ax.plot(self.magA['dateTime'][si-5:ei+5], self.dMLT[si-5:ei+5], 'b', label='dMLT')
            ax.axvline(self.magA['dateTime'][si])
            ax.axvline(self.magA['dateTime'][ei-1])
            plt.legend()
        return

    def _load_magephem(self, fPath, args):
        """
        This helper method loads in the magnetic ephemeris data,
        and converts the times to datetimes.
        """
        mag = dm.readJSONheadedASCII(fPath)
        
        # Copy over the usefull data, and standardize the dict keys
        magFlt = {}
        magFlt['dateTime'] = np.array([dateutil.parser.parse(t).replace(tzinfo=None)
                            for t in mag[args['time']]]) 

        # Loop over all of the following keys, and save to appropriate keys in
        # the magEphem dict.
        for key in ['L', 'MLT', 'lat', 'lon', 'alt']:
            # If true, take the specified column.
            if ''.join([key.lower(), 'col']) in [x.lower() for x in args.keys()]:
                magFlt[key] = mag[args[key]][:, args[key+'col']]
            else:
                magFlt[key] = mag[args[key]]
            if key == 'L': # Fold in negative (BLC) values.
                magFlt[key] = np.abs(magFlt[key])
            magFlt[key][magFlt[key] == 1.0e+31] = self.REPLACE_ERROR_VALS
        return magFlt

    def _find_common_times(self):
        """
        This method finds the common times between the two data sets
        and shrinks them to those times.
        """
        # Convert times to numbers
        tA = date2num(self.magA['dateTime'])
        tB = date2num(self.magB['dateTime'])

        idTA = np.in1d(tA, tB)
        idTB = np.in1d(tB, tA)
        # Try-except block to not accidently try to filter the an auxillary key.
        for key in self.magA.keys():
            try:
                self.magA[key] = self.magA[key][idTA]
            except KeyError as err:
                raise
        for key in self.magB.keys():
            try:
                self.magB[key] = self.magB[key][idTB]
            except KeyError as err:
                raise
        return

    def _calc_start_stop(self, ind):
        """ 
        This method calculates the start/stop indicies for a list 
        of indicies.
        """
        startInd = np.array([], dtype=int)
        endInd = np.array([], dtype=int)

        for k, g in itertools.groupby(enumerate(ind), lambda i: i[0]-i[1]):
            ind = list(map(operator.itemgetter(1), g))
            startInd = np.append(startInd, ind[0])
            endInd = np.append(endInd, ind[-1]+1)
        return startInd, endInd

    def _interp_conjunction(self, startInd, endInd, interpP):
        """
        This helper method uses the lat/lon/alt data in the
        magephem to interpolate it, and calculate L and MLT with
        T89 with Kp = 2 (default)
        """

        return

def dmlt(a, b):
    """
    NAME:    dmlt(a, b)
    USE:     Finds the absolute value of the difference
             of two MLT arrays a and b. This function
             correctly maps the differences over midnight
             Example: MLTa = 23, MLTb = 2 => dMLT = 3, NOT
             21!
    INPUT:   Two integers, or arrays to be differenced
    RETURNS: difference in MLT. 
    AUTHOR:  Mykhaylo Shumko
    MOD:     2017-04-30
    """
    # Convert the difference to spduo-angle.
    arg = 2*np.pi*np.abs(a - b)/24 
    # Utilize the even symmetry of cos to get 
    # the correct dmlt.
    return 24/(2*np.pi)*np.arccos(np.cos(arg))

if __name__ == '__main__':
    missionA = 'FIREBIRD'
    missionB = 'RBSP'
    fNameA = '20171130_FU4_T89_MagEphem.txt'
    fNameB = 'rbspa_def_MagEphem_T89D_20171130_v1.0.0.txt'
    # fNameA = '20171129_FU4_T89_MagEphem.txt'
    # fNameB = 'rbspa_def_MagEphem_T89D_20171129_v1.0.0.txt'
    fDirA = '/home/mike/research/firebird/Datafiles/FU_4/magephem/'
    fDirB = '/home/mike/research/rbsp/magephem/rbspa/'

    #magA = (os.path.join(fDirA, fNameA), 'dateTime', 'MLT', 'McllwainL')
    #magB = (os.path.join(fDirB, fNameB), 'DateTime', 'EDMAG_MLT', 'Lstar', -1)

    magA = os.path.join(fDirA, fNameA)
    magB = os.path.join(fDirB, fNameB)

    m = MagneticConjunctions(missionA, missionB,
        magA, magB)
    m.calc_conjunctions(interpP=())
    #m.testPlots()
    #plt.show()
