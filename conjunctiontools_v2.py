import time
from datetime import datetime, timedelta
import spacepy.datamodel as dm
import dateutil.parser
import multiprocessing
import numpy as np
import os, sys
import matplotlib.pylab as plt
from matplotlib.dates import date2num

import IRBEM

class MagneticConjunctions():
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

        ### Load data if necessary (assuming it is in a JSON headed ASCII file) ###
        if isinstance(magephemAargs, tuple):
            self.magA = self._load_magephem(*magephemAargs)
        if isinstance(magephemBargs, tuple):
            self.magB = self._load_magephem(*magephemBargs)
        self._find_common_times()
        return
    
    def calc_conjunctions(self, footpointAlt=None):
        """
        This method calculates the conjunctions, their duration, and
        if footpointAlt is set to an integer, then this method will
        calculate the smallest footpoint separation.
        """
        self.dL = np.abs(self.magA[self.magArgsA[3]] - 
                self.magB[self.magArgsB[3]])
        self.dMLT = dmlt(self.magA[self.magArgsA[2]], 
                self.magB[self.magArgsB[2]])

        # Calc indicies where separation meets conjunction criteria.
        idC = np.where((self.dL < self.Lthresh) & (self.dMLT < self.MLTthresh))[0]
        return

    def _load_magephem(self, fPath, tKey, MLTkey, Lkey, Lcol=None):
        """
        This helper method loads in the magnetic ephemeris data,
        and converts the times to datetimes.
        """
        print(fPath, tKey, MLTkey, Lkey, Lcol)
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
    M.calc_conjunctions()
