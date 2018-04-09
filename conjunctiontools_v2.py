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
        self.dL = np.abs(self.magA[self.LkeyA] - 
                self.magA[self.LkeyB])
        self.dMLT = dmlt(self.magB[self.MLTkeyA], 
                self.magB[self.MLTkeyB])

        # Calc indicies where separation meets conjunction criteria.
        idC = np.where((self.dL < self.Lthresh) & (self.dMLT < self.MLTthresh))[0]

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
            mag[Lkey] = mag[Lkey][L, Lcol]

        # Make all L shells positive (to avoid issues in the BLC)
        mag[Lkey] = np.abs(mag[Lkey])
        return

    def _find_common_times(self):
        """
        This method finds the common times between the two data sets
        and shrinks them to those times.
        """
        # Convert times to numbers
        tA = date2num(self.magA['dateTime'])
        tB = date2num(self.magB['dateTime'])

        idTA = np.in1d(tA, tB)
        idtB = np.in1d(tB, tA)

        for key in self.magA.keys():
            self.magA[key] = self.magA[key][idTA]
        for key in self.magB.keys():
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
    pass
