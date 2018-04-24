import time
from datetime import datetime, timedelta
import spacepy.datamodel as dm
import dateutil.parser
import multiprocessing
import numpy as np
import os, sys
import matplotlib.pylab as plt
from matplotlib.dates import date2num, num2date
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
        self.missionA = missionA
        self.missionB = missionB

        # magephem args       
        self.magArgsA = kwargs.get('magArgsA', keyDict[missionA])
        self.magArgsB = kwargs.get('magArgsB', keyDict[missionB])

        # Aux params
        self.REPLACE_ERROR_VALS = kwargs.get('REPLACE_ERROR_VALS', np.nan)
        self.Kp = kwargs.get('Kp', 20) # For interplation.

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
    
    def calc_conjunctions(self, interp=True):
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
        # Calc where indicies are continous
        self.startInd, self.endInd = self._calc_start_stop(idC) 
        ### Interpolate L and MLT around the flagged conjunctions. ###
        if interp:

            interpDict = self._interp_geo_pos()
        ### Calculate closest footpoint separation ###

        return

    def testPlots(self):
        """ 
        This method plots the dL and dMLT for each of the times 
        flagged for a conjunction.
        """
        for si, ei in zip(self.startInd, self.endInd):
            # Interpolate
            interpDict = self._interp_geo_pos(si-5, ei+5)

            fig, ax = plt.subplots(5, sharex=True, figsize=(6, 8))
            ax[0].plot(self.magA['dateTime'][si-5:ei+5], self.dL[si-5:ei+5], 'r', label='dL')
            ax[0].plot(self.magA['dateTime'][si-5:ei+5], self.dMLT[si-5:ei+5], 'b', label='dMLT')
            ax[0].set_ylabel('dL ad dMLT')

            ax[1].plot(self.magA['dateTime'][si-5:ei+5], self.magA['L'][si-5:ei+5],
                     'r', label='{}'.format(self.missionA))
            ax[1].plot(self.magA['dateTime'][si-5:ei+5], self.magB['L'][si-5:ei+5],
                     'b', label='{}'.format(self.missionB))
            ax[1].set_ylabel('L shell')

            ax[-3].plot(self.magA['dateTime'][si-5:ei+5], self.magA['alt'][si-5:ei+5],'r')
            ax[-3].plot(self.magA['dateTime'][si-5:ei+5], self.magB['alt'][si-5:ei+5],'b')
            ax[-3].scatter(interpDict['dateTime'], interpDict['altA'])
            ax[-3].scatter(interpDict['dateTime'], interpDict['altB'])
            ax[-3].set_ylabel('Altitude [km]')

            ax[-2].plot(self.magA['dateTime'][si-5:ei+5], self.magA['lat'][si-5:ei+5],'r')
            ax[-2].plot(self.magA['dateTime'][si-5:ei+5], self.magB['lat'][si-5:ei+5],'b')
            ax[-2].scatter(interpDict['dateTime'], interpDict['latA'])
            ax[-2].scatter(interpDict['dateTime'], interpDict['latB'])
            ax[-2].set_ylabel('latitude [deg]')

            ax[-1].plot(self.magA['dateTime'][si-5:ei+5], self.magA['lon'][si-5:ei+5],'r')
            ax[-1].plot(self.magA['dateTime'][si-5:ei+5], self.magB['lon'][si-5:ei+5],'b')
            ax[-1].scatter(interpDict['dateTime'], interpDict['lonA'])
            ax[-1].scatter(interpDict['dateTime'], interpDict['lonB'])
            ax[-1].set_ylabel('longitude [deg]')

            for a in ax:
                a.axvline(self.magA['dateTime'][si])
                a.axvline(self.magA['dateTime'][ei-1])
                a.legend(loc=1)
            plt.tight_layout()
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
            magFlt[key][np.abs(magFlt[key]) == 1.0e+31] = self.REPLACE_ERROR_VALS
        return magFlt

    def _find_common_times(self):
        """
        This method finds the common times between the two data sets
        and shrinks them to those times.
        """
        # Convert times to numbers
        self.tA = date2num(self.magA['dateTime'])
        self.tB = date2num(self.magB['dateTime'])

        idTA = np.in1d(self.tA, self.tB) # Efficient way to calculate the same times
        idTB = np.in1d(self.tB, self.tA) # between two data sets.
        
        for key in self.magA.keys(): # Filter data.
            self.magA[key] = self.magA[key][idTA]
            self.magB[key] = self.magB[key][idTB]
        self.tA = self.tA[idTA]
        self.tB = self.tB[idTB]
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

    def _calc_d_dMLT_param(self, startInd, endInd, alt=500):
        """
        This method is a wrapper for the interpolation in lat/lon/alt
        coordinates. This method also calculates the smallest footprint
        separation at an altitude alt. Lastly, it calculates dMLT when 
        L shells cross.
        """
        # get interpolated lat/lon/alt
        interpDict = self._interp_geo_pos(startInd, endInd)
        #Xa = 
        fpotNA = self.find_foot_point(X, {'Kp':self.Kp})
        return

    def _interp_geo_pos(self, startInd, endInd):
        """
        This helper method interpolates the lat/lon/alt data. The longitude
        coordinate is treated separately.
        """
        tInterp = np.linspace(self.tA[startInd], self.tA[endInd-1]) 
        interpDict = {'t':tInterp}

        flatA = scipy.interpolate.interp1d(self.tA[startInd:endInd],
            self.magA['lat'][startInd:endInd], kind='cubic')
        flatB = scipy.interpolate.interp1d(self.tB[startInd:endInd],
            self.magB['lat'][startInd:endInd], kind='cubic')
        faltA = scipy.interpolate.interp1d(self.tA[startInd:endInd],
            self.magA['alt'][startInd:endInd], kind='cubic')
        faltB = scipy.interpolate.interp1d(self.tB[startInd:endInd],
            self.magB['alt'][startInd:endInd], kind='cubic')
        interpDict['latA'] = flatA(tInterp)
        interpDict['latB'] = flatB(tInterp)
        interpDict['altA'] = faltA(tInterp)
        interpDict['altB'] = faltB(tInterp)

        # Be carefull with lons over the 180 degree boundary.
        interpDict['lonA'], interpDict['lonB'] = self._interp_lon(
                                        tInterp, startInd, endInd)
        interpDict['dateTime'] = num2date(tInterp)
        return interpDict

    def _interp_lon(self, tInterp, startInd, endInd):
        """

        """
        lonA = self.magA['lon'][startInd:endInd].copy()
        lonB = self.magB['lon'][startInd:endInd].copy()
        dLonA = lonA[1:] - lonA[:-1]
        dLonB = lonB[1:] - lonB[:-1]

        LonA_idx = np.where(dLonA > 100)[0] + 1
        LonB_idx = np.where(dLonB > 100)[0] + 1
        #+1 accounts for different array sizes
        for i in LonA_idx:
            lonA[i:] -= 360
            # dLonA[i:] -= 360
        for i in LonB_idx:
            # dLonB[i:] -= 360
            lonB[i:] -= 360
        
        flonA = scipy.interpolate.interp1d(self.tA[startInd:endInd],
            lonA, kind='cubic')
        flonB = scipy.interpolate.interp1d(self.tB[startInd:endInd],
            lonB, kind='cubic')
        flonA, flonB = flonA(tInterp), flonB(tInterp)
        while np.min(flonA) < -180:
            idx = np.where(flonA < -180)[0][0]
            flonA[idx:] += 360
        while np.min(flonB) < -180:
            idx = np.where(flonB < -180)[0][0]
            flonB[idx:] += 360
        return flonA, flonB

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
    magA = os.path.join(fDirA, fNameA)
    magB = os.path.join(fDirB, fNameB)

    m = MagneticConjunctions(missionA, missionB,
        magA, magB)
    m.calc_conjunctions(interpP=())
    m.testPlots()
    plt.show()
