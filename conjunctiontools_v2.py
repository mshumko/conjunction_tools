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
import spacepy
import csv

import IRBEM

Re=6371 # Earth radius km

# This dictionary contains the magEphem L, MLT, lat, lon, alt, 
# time keys for reading in the data. 
keyDict = {
        'FIREBIRD':{
            'L':'McllwainL', 'MLT':'MLT', 'time':'dateTime', 
            'lat':'Rgeod_LatLon', 'latcol':0, 'lon':'Rgeod_LatLon', 
            'loncol':1, 'alt':'Rgeod_Altitude'
            },
        
        'DSX':{
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
        self.lowerL = kwargs.get('lowerL', 3)
        self.Lthresh = kwargs.get('Lthresh', 1)
        self.MLTthresh = kwargs.get('MLTthresh', 1)
        self.verbose = kwargs.get('verbose', False)
        self.missionA = missionA
        self.missionB = missionB

        # magephem args       
        self.magArgsA = kwargs.get('magArgsA', keyDict[missionA])
        self.magArgsB = kwargs.get('magArgsB', keyDict[missionB])

        # Aux params
        self.REPLACE_ERROR_VALS = kwargs.get('REPLACE_ERROR_VALS', np.nan)
        #self.Kp = kwargs.get('Kp', 20) # Deault kp for interplation.

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
    
    def calcConjunctions(self, interp=True):
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
        idC = np.where((self.dL < self.Lthresh) & 
                        (self.dMLT < self.MLTthresh) &
                        (self.magA['L'] > self.lowerL) & 
                        (self.magB['L'] > self.lowerL))[0]
        # Calc where indicies are continous
        self.startInd, self.endInd = self._calc_start_stop(idC) 
        ### Interpolate L and MLT around the flagged conjunctions. ###
        ### Calculate closest footpoint separation ###
        if interp:
            self._calc_d_dMLT_param()        
        else:
            raise NotImplementedError('Simple non-interp mode not implemented.')
        return

    def saveData(self, sPath, mode='a'):
        """
        This method saves the conjunction data to a csv file.
        """
        saveData = np.stack((self.startTime, self.endTime, 
                            self.meanL, self.meanMLT, self.minMLT, 
                            self.dmin), axis=1)
        # Delete rows with error values that are -1
        invalidRows = np.where(saveData[:, 0] == -1)[0]
        saveData = np.delete(saveData, invalidRows, axis=0)
        exists = os.path.exists(sPath)
        with open(sPath, mode, newline='') as f:
            w = csv.writer(f)
            # Save header
            if not exists: # If the file is newly generated, write the header.
                w.writerow(['startTime', 'endTime', 'meanL', 
                            'meanMLT', 'minMLT', 'minD [km]'])
            w.writerows(saveData)
        return

    def testPlots(self):
        """ 
        This method plots the dL and dMLT for each of the times 
        flagged for a conjunction.
        """
        for i, (si, ei) in enumerate(zip(self.startInd, self.endInd)):
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
                a.axvline(self.magA['dateTime'][si], c='k', label='Approx bounds')
                a.axvline(self.magA['dateTime'][ei-1], c='k')
                a.axvline(self.startTime[i], c='r', label='Bounds')
                a.axvline(self.endTime[i], c='r')
                a.legend(loc=1)
            plt.tight_layout()
        return

    def plot_overview(self):
        """
        Plots the overview of the L and MLT of both spacecraft. Also plots the 
        dL and dMLT.
        """
        fig, ax = plt.subplots(3, sharex=True)
        
        ax[0].plot(self.magA['dateTime'], self.magA['L'],
                     'r', label='{}'.format(self.missionA))
        ax[0].plot(self.magB['dateTime'], self.magB['L'],
                     'b', label='{}'.format(self.missionB))
        ax[0].set_ylabel('L')
                     
        ax[1].plot(self.magA['dateTime'], self.magA['MLT'],
                     'r', label='{}'.format(self.missionA))
        ax[1].plot(self.magB['dateTime'], self.magB['MLT'],
                     'b', label='{}'.format(self.missionB))
        ax[1].set_ylabel('MLT')
                     
        ax[2].plot(self.magA['dateTime'], self.dL, 'k', label='dL')
        ax[2].plot(self.magB['dateTime'], self.dMLT,'k--', label='dMLT')
        ax[2].set_ylabel('Difference')    
        
        ax[0].legend()
        ax[2].legend()
        plt.show()
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

    def _calc_d_dMLT_param(self, alt=500, interpOffset=2):
        """
        This method is a wrapper for the interpolation in lat/lon/alt
        coordinates. This method also calculates the smallest footprint
        separation at an altitude alt. Lastly, it calculates dMLT when 
        L shells cross.
        """
        self.dmin = np.nan*np.ones_like(self.startInd)
        self.MLTmin = np.nan*np.ones_like(self.startInd)
        self.startTime = np.nan*np.ones(len(self.startInd), dtype=object)
        self.endTime = np.nan*np.ones(len(self.startInd), dtype=object)
        self.minMLT = np.nan*np.ones_like(self.startInd)
        self.meanMLT = np.nan*np.ones_like(self.startInd)
        self.meanL = np.nan*np.ones_like(self.startInd)

        for ci, (si, ei) in enumerate(zip(self.startInd, self.endInd)):
            # get interpolated lat/lon/alt
            if si < interpOffset: # If a conjunction starts at the start of file (rare but happens)
                si = interpOffset
                ei += interpOffset
            interpDict = self._interp_geo_pos(si-interpOffset, ei+interpOffset) # Expand to interpolate
            if interpDict == -1:
                return
            t0 = self.magA['dateTime'][si]
            footpointNA = np.nan*np.ones((len(interpDict['latA']), 3), dtype=float)
            footpointSA = np.nan*np.ones((len(interpDict['latA']), 3), dtype=float)
            footpointNB = np.nan*np.ones((len(interpDict['latA']), 3), dtype=float)
            footpointSB = np.nan*np.ones((len(interpDict['latA']), 3), dtype=float)

            LA = np.nan*np.ones(len(interpDict['latA']), dtype=float)
            LB = np.nan*np.ones(len(interpDict['latA']), dtype=float)
            MLTA = np.nan*np.ones(len(interpDict['latA']), dtype=float)
            MLTB = np.nan*np.ones(len(interpDict['latA']), dtype=float)

            # Run IRBEM find_foot_point()
            zA = zip(interpDict['dateTime'], interpDict['altA'], 
                    interpDict['latA'], interpDict['lonA'])
            for (i, (ti, alti, lati, loni)) in enumerate(zA):
                # Get footprint params.
                self.find_foot_point(
                    {'dateTime':ti, 'x1':alti, 'x2':lati, 'x3':loni}, 
                    {'Kp':self.getKp(t0)}, alt, 1)
                footpointNA[i] = self.find_foot_point_output['XFOOT']
                self.find_foot_point(
                    {'dateTime':ti, 'x1':alti, 'x2':lati, 'x3':loni}, 
                    {'Kp':self.getKp(t0)}, alt, -1)
                footpointSA[i] = self.find_foot_point_output['XFOOT']
                # Get L and MLT params
                self.make_lstar(
                    {'dateTime':ti, 'x1':alti, 'x2':lati, 'x3':loni}, 
                    {'Kp':self.getKp(t0)})
                LA[i] = self.make_lstar_output['Lm'][0]
                MLTA[i] = self.make_lstar_output['MLT'][0]

            zB = zip(interpDict['dateTime'], interpDict['altB'], 
                    interpDict['latB'], interpDict['lonB'])
            for (i, (ti, alti, lati, loni)) in enumerate(zB):
                self.find_foot_point(
                    {'dateTime':ti, 'x1':alti, 'x2':lati, 'x3':loni}, 
                    {'Kp':self.getKp(t0)}, alt, 1)
                footpointNB[i] = self.find_foot_point_output['XFOOT']
                self.find_foot_point(
                    {'dateTime':ti, 'x1':alti, 'x2':lati, 'x3':loni}, 
                    {'Kp':self.getKp(t0)}, alt, -1)
                footpointSB[i] = self.find_foot_point_output['XFOOT']
                # Get L and MLT params
                self.make_lstar(
                    {'dateTime':ti, 'x1':alti, 'x2':lati, 'x3':loni}, 
                    {'Kp':self.getKp(t0)})
                LB[i] = self.make_lstar_output['Lm'][0]
                MLTB[i] = self.make_lstar_output['MLT'][0]
            
            # Calc the closst separation.
            dN = greatCircleDist(footpointNA, footpointNB)
            dS = greatCircleDist(footpointSA, footpointSB)
            try:
                self.dmin[ci] = np.min(np.append(dN[dN > 0], dS[dS > 0]))
            except ValueError as err:
                if str(err) == 'zero-size array to reduction operation minimum which has no identity':
                    if self.verbose:
                        print('WARNING: dmin not found (IRBEM error values)')
                    self.dmin[ci] == -1E31
                else: 
                    raise
            # Calculate start/stop times as well as dMLT at closest L
            # shell separation.
            self.startTime[ci], self.endTime[ci], self.minMLT[ci], idx = self._find_c_bounds(
                interpDict['dateTime'], LA, MLTA, LB, MLTB
                )   
            # Save mean values of L and center MLT values at the closest approach.
            self.meanL[ci] = 0.5*(np.mean(np.abs(LA[idx])) + np.mean(np.abs(LB[idx])))
            if hasattr(idx, '__len__'):
                self.meanMLT[ci] = MLTA[idx[len(idx)//2]]
            else:
                self.meanMLT[ci] = MLTA[idx]
        return

    def _interp_geo_pos(self, startInd, endInd):
        """
        This helper method interpolates the lat/lon/alt data. The longitude
        coordinate is treated separately.
        """
        try:
            tInterp = np.linspace(self.tA[startInd], self.tA[endInd-1]) 
        except IndexError as err:
            if 'is out of bounds for axis' in str(err):
                # Small chance that the end time conjunction index is at 
                # the end of the data set.
                return -1
            else:
                raise
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
        This method interpolates the longitude assuming that the 
        spacecraft lon motion is westward.
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

    def _find_c_bounds(self, time, LA, MLTA, LB, MLTB):
        """
        This method will calculate the start/stop times of a conjunction,
        as well as dMLT when the L shells are closest (usually near a 
        crossing). Lastly, it returns the index of the LA/LB/MLTA/MLTB 
        arrays during closest approach. 
        """
        dL = np.abs(np.abs(LA)-np.abs(LB))
        #print(LA, LB)
        # Order does not matter since dmlt() does not distinguish sign.
        dMLT = dmlt(MLTA, MLTB) 
        validInd = np.where((dL < self.Lthresh) & (dMLT < self.MLTthresh))[0]
        # If the flagged conjunction was not really a conjunction when interpolated
        if len(validInd) == 0:
            return (-1, -1, -1, -1)
        startTime = time[validInd[0]].replace(tzinfo=None)
        endTime = time[validInd[-1]].replace(tzinfo=None)
        minMLT = dMLT[int(np.argmin(dL))]
        idx = int(np.argmin(dL))
        return startTime, endTime, minMLT, idx

    def getKp(self, t, default=20, kpDir='/home/mike/research/firebird/data_processing/geomag_indicies/indicies'):
        """
        This method will open up and append the kp indicies for the time range
        in the magEphem file. Given a time t, it will look for the correct kp.
        """
        if not hasattr(self, 'kp'): # Load in the files first time this is called
            years = sorted(set([t.year for t in self.magA['dateTime']]))
            self.kp = {'dateTime':np.array([]), 'kp':np.array([])}
            for year in years:
                f = spacepy.datamodel.readJSONheadedASCII(os.path.join(kpDir, '{}_kp.txt'.format(year)))
                self.kp['dateTime'] = np.append(self.kp['dateTime'], f['dateTime'])
                self.kp['kp'] = np.append(self.kp['kp'], f['kp'])
            # Convert time strings to datetime
            self.kp['dateTime'] = np.array([dateutil.parser.parse(tt) for tt in self.kp['dateTime']])
        # Now find the matching Kp
        t = t.replace(hour=(t.hour - t.hour%3), minute=0, second=0, microsecond=0) # Convert to a multiple of 3 hrs.
        idt = np.where(self.kp['dateTime'] == t)[0]

        if len(idt) == 0:
            if self.verbose:
                print('WARNING: No kp value found at time {}. '
                      'Returning default kp of {}'.format(t, default))
            return 20
        return self.kp['kp'][idt[0]]

def greatCircleDist(X1, X2):
    """
    X1 and X2 must be N*3 array of lat, lon, alt. phi = lat, lambda = lon
    """ 
    X1 = np.asarray(X1)
    X2 = np.asarray(X2)
    R = (Re+(X1[:, 0]+X2[:, 0])/2)
    s = 2*np.arcsin( np.sqrt( np.sin(np.deg2rad(X1[:, 1]-X2[:, 1])/2)**2 + np.cos(np.deg2rad(X1[:, 1]))*np.cos(np.deg2rad(X2[:, 1]))*np.sin(np.deg2rad(X1[:, 2]-X2[:, 2])/2)**2 ))
    return R*s

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
    # Find the difference (only correct if dMLT < 12!).
    dMLT = np.abs(a - b)
    # Correct the dMLT > 12 values by looking at the reciprical
    dMLT[dMLT > 12] = 24 - dMLT[dMLT > 12]
    return dMLT

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
    m.calcConjunctions()
    m.saveData('/home/mike/Desktop/test.csv')
    m.testPlots()
    plt.show()
