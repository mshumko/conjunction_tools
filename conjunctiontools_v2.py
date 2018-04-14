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

class MagneticConjunctions(IRBEM.MagFields):
    def __init__(self, magephemAargs, magephemBargs, **kwargs):
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
        self.magArgsA = magephemAargs
        self.magArgsB = magephemBargs

        # Aux params
        self.REPLACE_ERROR_VALS = kwargs.get('REPLACE_ERROR_VALS', np.nan)
        self.Kp = kwargs.get('Kp', 2) # For interplation.

        # Initializa IRBEM
        IRBEM.MagFields.__init__(self)

        ### Load data if necessary (assuming it is in a JSON headed ASCII file) ###
        if isinstance(magephemAargs, tuple):
            self.magA = self._load_magephem(*magephemAargs)
        if isinstance(magephemBargs, tuple):
            self.magB = self._load_magephem(*magephemBargs)
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
        self.dL = np.abs(self.magA[self.magArgsA[3]] - 
                self.magB[self.magArgsB[3]])
        self.dMLT = dmlt(self.magA[self.magArgsA[2]], 
                self.magB[self.magArgsB[2]])

        # Calc indicies where separation meets conjunction criteria.
        idC = np.where((self.dL < self.Lthresh) & (self.dMLT < self.MLTthresh))[0]
        #idC = idC[:-1]
        #print(idC)
        startInd, endInd = self._calc_start_stop(idC)



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

    def _load_magephem(self, fPath, tKey, MLTkey, Lkey, Lcol=None):
        """
        This helper method loads in the magnetic ephemeris data,
        and converts the times to datetimes.
        """
        mag = dm.readJSONheadedASCII(fPath)
        # Convert times
        mag[tKey] = np.array([dateutil.parser.parse(t) for t in mag[tKey]])

        # Replace L and MLT errors values.
        for key in [Lkey, MLTkey]:
            mag[key][mag[key] == -1.0e+31] = self.REPLACE_ERROR_VALS

        if Lcol is not None:
            mag[Lkey] = mag[Lkey][:, Lcol]

        # Make all L shells positive (to avoid issues in the BLC)
        mag[Lkey] = np.abs(mag[Lkey])
        return mag

    def _find_common_times(self):
        """
        This method finds the common times between the two data sets
        and shrinks them to those times.
        """
        # Convert times to numbers
        tA = date2num(self.magA[self.magArgsA[1]])
        tB = date2num(self.magB[self.magArgsB[1]])

        idTA = np.in1d(tA, tB)
        idTB = np.in1d(tB, tA)

        for key in self.magArgsA[1:4]:
            self.magA[key] = self.magA[key][idTA]
        for key in self.magArgsB[1:4]:
            self.magB[key] = self.magB[key][idTB]
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
    fNameA = '20171130_FU4_T89_MagEphem.txt'
    fNameB = 'rbspa_def_MagEphem_T89D_20171130_v1.0.0.txt'
    # fNameA = '20171129_FU4_T89_MagEphem.txt'
    # fNameB = 'rbspa_def_MagEphem_T89D_20171129_v1.0.0.txt'
    fDirA = '/home/mike/research/firebird/Datafiles/FU_4/magephem/'
    fDirB = '/home/mike/research/rbsp/magephem/rbspa/'

    magA = (os.path.join(fDirA, fNameA), 'dateTime', 'MLT', 'McllwainL')
    magB = (os.path.join(fDirB, fNameB), 'DateTime', 'EDMAG_MLT', 'Lstar', -1)

    M = MagneticConjunctions(magA, magB)
    M.calc_conjunctions(interpP=())
