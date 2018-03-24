#import argparse
#parser = argparse.ArgumentParser()
#parser.add_argument('-l', '--license', help="Print a short licence", action="store_true")
#args = parser.parse_args()
#if args.license:
#    print('conjunctionToolkit  Copyright (C) 2017  Mykhaylo Shumko', 
#        'This program comes with ABSOLUTELY NO WARRANTY. This is free software',
#         ', and you are welcome to redistribute it under the license conditions')
 
import time
from datetime import datetime, timedelta
import spacepy.datamodel as dm
import dateutil.parser
import multiprocessing
import numpy as np
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname( __file__ ), 'tests'))
import duration_test

import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
#import matplotlib.patches as mpatches

class MagneticConjunctionCalc:
    def __init__(self, fDirA, fnameA, fDirB, fnameB, **kwargs):
        """
        NAME:    ConjunctionCalculator()
        USE:     A conjunction calculator class that calculates conjuntions
                 between two sets of magnetic ephemeris data with some sort
                 of L shell and MLT. The plotting function 
        INPUT:   REQUIRED ARGS: file directories and filenames of the two
                 magnetic ephemeris data sets.
                 OPTIONAL KWARGS:
                 verbose -------- Default False, prints status messages.
                 REPLACE_ERROR_VALS -Default np.inf, Eror message for
                                     proper handling of valid data.
                 mission_id_A or B - Default 'missionA or B', is the 
                                     mission id string used in saving 
                                     filenames.
                 lowerL ------------ Default 1, lower L shell of 
                                     conjunctions. Usefull for 
                                     calculating conjunctions
                                     only in the outer belt.
                 Lthresh ----------- Default is 1, threshold L seperation 
                                     for conjunction so that if 
                                     |dL| < Lthresh, its a conjunction.
                 MLTthresh --------- Default is 1, threshold MLT seperation
                                     for conjunction so that if
                                     |dMLT| < MLTthresh, its a conjunction.
                 cadence ----------- Default is 1 minute, cadence of the
                                     magnetic ephemeris data. 
                                     NOTE: the cadence of 
                                     the two data sets must be the same!
                 timeKeyA or B ----- Default 'dateTime', the key values
                                     of the time dictionary of the data.
                 LkeyA or B -------- Deault 'McllwainL', the key values
                                     of the L dictionary of the data.
                 MLTkeyA or B ------ Deault 'MLT', the key values
                                     of the MLT dictionary of the
                                     data.                     
        EXAMPLE: 
            fDir = << INSERT YOUR DIR HERE >>
            fnameA = << INSERT YOUR FILENAME A HERE >>
            fnameB = << INSERT YOUR FILENAME B HERE >>
            save_dir = << INSERT YOUR SAVE DIR HERE >>
            save_name = << INSERT YOUR SAVE DIR HERE >>
            lowerLBound = 3 # Lower L shell bound
            cCalc = ConjunctionCalculator(fDir, fnameA, fDir, fnameB, 
              mission_id_A = 'spaceshipA', L_thresh = 1, MLT_thresh = 1,
              mission_id_B = 'spaceshipB', lowerL = lowerLBound,
              verbose = True)
            
            cCalc.match_start_end_ind()
            cCalc.calc_magnetic_seperation()
            cCalc.calc_indicies_L_MLT_conjunction()
            cCalc.lower_L_bound()
            cCalc.calc_conjunction_duration()
            cCalc.save_to_file(save_dir = save_dir, save_name = save_name)
            cCalc.plotLMLTandDLDMLT()
        RETURNS: None
        AUTHOR:  Mykhaylo Shumko
        MOD:     2017-05-01
        """
        self.verbose = kwargs.get('verbose', False)
        self.REPLACE_ERROR_VALS = kwargs.get('REPLACE_ERROR_VALS', np.inf)
        self.mission_id_A = kwargs.get('mission_id_A', 'missionA')
        self.mission_id_B = kwargs.get('mission_id_B', 'missionB')
        self.lowerL = kwargs.get('lowerL', 1)
        self.L_thresh = kwargs.get('Lthresh', 1)
        self.cadence = kwargs.get('cadence', 1)
        self.MLT_thresh = kwargs.get('MLTthresh', 1)

        # Optionally, provide the data in dictionary form.
        self.magephemA = kwargs.get('magEphemA', None)
        self.magephemB = kwargs.get('magEphemB', None)

        # Optionally set the data keys
        self.timeKeyA = kwargs.get('timeKeyA', 'dateTime')
        self.timeKeyB = kwargs.get('timeKeyB', 'dateTime')
        self.LkeyA = kwargs.get('LkeyA', 'McllwainL')
        self.LkeyB = kwargs.get('LkeyB', 'McllwainL')
        self.MLTkeyA = kwargs.get('MLTkeyA', 'MLT')
        self.MLTkeyB = kwargs.get('MLTkeyB', 'MLT')
        self.LcolA = kwargs.get('LcolA', None)
        self.LcolB = kwargs.get('LcolB', None)
        
        if self.verbose:
            print('Loading Magnetic Ephemeris Data')
        
        # Assume if user does not supply the data, that it is in JSON format.
        if self.magephemA is None:
            self.magephemA = dm.readJSONheadedASCII(os.path.join(fDirA, fnameA))
        if self.magephemB is None:
            self.magephemB = dm.readJSONheadedASCII(os.path.join(fDirB, fnameB))

        if self.verbose:
            print('Converting times to datetime objects')
        
        # Use multiprocessing module to speed up the time conversion.
        with multiprocessing.Pool() as p:
            if not isinstance(self.magephemA[self.timeKeyA][0], datetime):
                self.magephemA[self.timeKeyA] = np.array(list(
                p.map(dateutil.parser.parse, self.magephemA[self.timeKeyA])))
            if not isinstance(self.magephemB[self.timeKeyB][0], datetime):
                self.magephemB[self.timeKeyB] = np.array(list(
                p.map(dateutil.parser.parse, self.magephemB[self.timeKeyB])))

        # Change the invalid magnetic coordinates to NAN, and 
        # make L shell absolute value. (This is for numpy math 
        # operations)
        for key in [self.LkeyA, self.MLTkeyA]:
            self.magephemA[key][
            self.magephemA[key] == -1.0e+31] = self.REPLACE_ERROR_VALS
            
        for key in [self.LkeyB, self.MLTkeyB]:
            self.magephemB[key][
            self.magephemB[key] == -1.0e+31] = self.REPLACE_ERROR_VALS
            
        # Optionally, load in the column with a specified pitch angle.
        if self.LcolA is not None:
            self.magephemA[self.LkeyA] = self.magephemA[self.LkeyA][:, self.LcolA]
        if self.LcolB is not None:
            self.magephemB[self.LkeyB] = self.magephemB[self.LkeyB][:, self.LcolB]
            
        # Make L all positive
        self.magephemA[self.LkeyA] = np.abs(self.magephemA[self.LkeyA])
        self.magephemB[self.LkeyB] = np.abs(self.magephemB[self.LkeyB])
        
        if self.verbose: print('Times converted. Done initializing')
            
    def downsampleData(self, dataID, sample_every = 5):
        """
        NAME:      downsampleData(self, dataID, sample_every = 5) 
        USE:          Downsample all of the data products in a magnetic ephemeris file 
                          identified with the dataID argument.
        INPUT:      dataID either 'A', 'B' (lowercase or uppercase). sample_every, 
                          default to 5, this samples that data set every sample_every points.
        RETURNS: None (Downsamples all of the magnetic ephemeis data)
        AUTHOR:  Mykhaylo Shumko
        MOD:         2017-05-30
        """
        if self.verbose: print('Downsampling data')
        
        assert dataID.upper() in ['A', 'B'], 'ERROR, dataID must be either "A" or "B".'
        if dataID.upper() == 'A':
            for key in list(self.magephemA.keys()):
                self.magephemA[key] = self.magephemA[key][::sample_every]
        elif dataID.upper() == 'B':
            for key in list(self.magephemB.keys()):
                self.magephemB[key] = self.magephemB[key][::sample_every]
        return
        
    def periodic_data_indicies(self):
        """
        NAME:      match_start_end_ind(self)
        USE:       This function will find a common start and end indicies and 
                   shorten the magnetic ephemeris files approprietely. This will 
                   allow for fast comparison. THIS FUNCTION ASSUMES THAT 
                   BOTH MAGNETIC EPHEMERIS FILES ARE AT AN 
                   IDENTICAL CADENCE! 
        INPUT:     None
        RETURNS: None (Shortens all of the magnetic ephemeis data)
        AUTHOR: Mykhaylo Shumko
        MOD:     2017-05-30
        """
        # Find common start and end times.
        if self.verbose: print('Finding common start and end times.' , 
            'Assuming cadence is the same')
                
        # Intersect1d finds all of the common times.
        sharedTimes = np.intersect1d(self.magephemA[self.timeKeyA], 
        self.magephemB[self.timeKeyB], assume_unique = True)

        assert len(sharedTimes) > 0, 'Error, no shared times found!' 

        startA = np.where(self.magephemA[self.timeKeyA] == sharedTimes[0])[0][0]
        startB = np.where(self.magephemB[self.timeKeyB] == sharedTimes[0])[0][0]
        endA = np.where(self.magephemA[self.timeKeyA] == sharedTimes[-1])[0][0]
        endB = np.where(self.magephemB[self.timeKeyB] == sharedTimes[-1])[0][0]
                
        # Now concatinate the magnetic ephemeris files for 1 to 1 comparion!
        if self.verbose: print('Reducing the data with the common times.')
        for key in list(self.magephemA.keys()):
            self.magephemA[key] = self.magephemA[key][startA:endA]
        for key in list(self.magephemB.keys()):
            self.magephemB[key] = self.magephemB[key][startB:endB]
            
        if self.verbose: print('Data reduced for 1-to-1 comparison.')
        return
    
    def calc_magnetic_seperation(self):
        """
        NAME:       calc_magnetic_seperation(self)
        USE:          This function calculates the L and MLT seperation as a function 
                           of time for the two ephemeris files. Make sure that the two 
                           magnetic ephemeris data sets are synced!
        INPUT:       None
        RETURNS: None (self.dL and self.dMLT class attributes)
        AUTHOR:   Mykhaylo Shumko
        MOD:          2017-05-30
        """
        try:
            self.dL = np.abs(self.magephemA[self.LkeyA] - 
                self.magephemB[self.LkeyB])
            self.dMLT = dmlt(self.magephemA[self.MLTkeyA], 
                self.magephemB[self.MLTkeyB])
        except ValueError as err:
            print('Number of data points in magephems do not match! Is the cadence', 
                'the same? Original error:')
            raise
        return
    
    def calc_indicies_L_conjunction(self, L_thresh):
        """
        NAME:      
        USE:          
        INPUT:      
        RETURNS: 
        AUTHOR: Mykhaylo Shumko
        MOD:     2017-05-30
        """
        
        """
        This function calculates and returns indicies where the difference in 
        L is smaller than L_thresh function arguments.
        """
        self.ind = np.where(self.dL < L_thresh)[0]
        return self.ind
        
    def calc_indicies_MLT_conjunction(self, MLT_thresh):
        """
        NAME:      
        USE:          
        INPUT:      
        RETURNS: 
        AUTHOR: Mykhaylo Shumko
        MOD:     2017-05-30
        """
        
        """
        This function calculates and returns indicies where the difference in 
        MLT is smaller than MLT_thresh function arguments.
        """
        self.ind = np.where(self.dMLT < MLT_thresh)[0]
        return self.ind
    
    def calc_indicies_L_MLT_conjunction(self, **kwargs):
        """
        NAME:  calc_indicies_L_MLT_conjunction(self, **kwargs)    
        USE:     This function calculates and returns indicies where the difference in 
                      MLT is smaller than self.MLT_thresh and L is smaller than 
                      self.L_thresh.  
        INPUT:  None (can change the self.L_thresh and self.MLT_thresh class 
                      attributes)
        RETURNS: self.ind, the indicies where the L and MLT conditions are satisfied
        AUTHOR: Mykhaylo Shumko
        MOD:        2017-05-30
        """
        self.L_thresh = kwargs.get('Lthresh', self.L_thresh)
        self.MLT_thresh = kwargs.get('MLTthresh', self.MLT_thresh)
        self.ind = np.where((self.dMLT < self.MLT_thresh) & (self.dL < self.L_thresh))[0]
        return self.ind
        
    def calc_conjunction_duration(self, ind = None, **kwargs):
        """
        NAME:      calc_conjunction_duration(self, ind = None, **kwargs)
        USE:          Calculates the duration of conjunctions. Duration and 
                           self.cadence should be in minutes.
        INPUT:      Optional: ind is the data indicies of conjunctions (this function 
                          looks for consecutive numbers for the duration calculation)
        RETURNS: self.durationDict, a dictionary of conjunction start time, end time, 
                           duration, and mean L and mean MLT for both spacecraft
        AUTHOR: Mykhaylo Shumko
        MOD:     2017-05-30
        """
        if self.verbose: print('Calculating conjunction durations')
        self.cadence = kwargs.get('cadence', self.cadence)
        self.durationDict = {}
        if ind is None:
            ind = self.ind
        assert len(ind) > 0, 'No conjunctions found for specified time interval!'
                
        # Find the consecutive times in the indicies.
        #startInd, endInd = locateConsecutiveNumbers(ind)
        tDict = duration_test.time_continuity(
            self.magephemA[self.timeKeyA][ind],
            self.cadence*60)
        startInd = tDict['startTime'].astype(int)
        endInd = tDict['endTime'].astype(int)

        self.durationDict['startTime'] = np.array([datetime.min for 
        i in range(len(startInd))])
        self.durationDict['endTime'] = np.array([datetime.min for 
        i in range(len(startInd))])
        
        self.durationDict['duration'] = -9999*np.ones(len(startInd))
        self.durationDict['meanL_' + self.mission_id_A] = -9999*np.ones(len(startInd))
        self.durationDict['meanL_' + self.mission_id_B] = -9999*np.ones(len(startInd))
        self.durationDict['meanMLT_' + self.mission_id_A] = -9999*np.ones(len(startInd))
        self.durationDict['meanMLT_' + self.mission_id_B] = -9999*np.ones(len(startInd))
        
        for i in range(len(startInd)):
            # Dumb bug... FIX this later!
            if endInd[i] == startInd[i]: 
                continue
                
            # Save the start time, end time, and duration of each conjunction.
            self.durationDict['startTime'][i] = self.magephemA[self.timeKeyA][
                ind[startInd[i]:endInd[i]]][0]
                
            # The offset by the cadence here is to make the end time of the
            # conjunction to be over at that point, and not a minute after.
            self.durationDict['endTime'][i] = (self.magephemA[self.timeKeyA][
                ind[startInd[i]:endInd[i]]][-1] + 
                timedelta(minutes = self.cadence))
            self.durationDict['duration'][i] = self.cadence*(
                len(ind[startInd[i]:endInd[i]]))
                
            # Save the mean value of L and MLT.
            self.durationDict[ 'meanL_' + self.mission_id_A][i] = np.mean(
                self.magephemA[self.LkeyA][ind[startInd[i]:endInd[i]]])
            self.durationDict[ 'meanL_' + self.mission_id_B][i] = np.mean(
                self.magephemB[self.LkeyB][ind[startInd[i]:endInd[i]]])
            self.durationDict[ 'meanMLT_' + self.mission_id_A][i] = np.mean(
                self.magephemA[self.MLTkeyA][ind[startInd[i]:endInd[i]]])
            self.durationDict[ 'meanMLT_' + self.mission_id_B][i] = np.mean(
                self.magephemB[self.MLTkeyB][ind[startInd[i]:endInd[i]]])

        # Look for, and remove error values
        validInd = np.where(self.durationDict['duration'] != -9999)[0]
        #print(validInd)
        for key in list(self.durationDict.keys()):
            self.durationDict[key] = self.durationDict[key][validInd]
  
        return self.durationDict
        
    def lower_L_bound(self):
        """
        NAME:     lower_L_bound(self, lowerL)
        USE:        Applies a lower L shell cutoff to the conjunctions.        
        INPUT:    None (change the theshold with self.lowerL)
        RETURNS: self.ind, a class instance of indicies that is an intersection of
                        the previous conjunction indicies and the indicies with the lower L
                        shell threshold.
        AUTHOR: Mykhaylo Shumko
        MOD:     2017-05-30
        """
        threshInd = np.where((self.magephemA[self.LkeyA] > self.lowerL) & 
                    (self.magephemB[self.LkeyB] > self.lowerL))[0]
        self.ind = np.array(sorted(list(set(threshInd) & set(self.ind))))
        return self.ind
        
    def plotLMLTandDLDMLT(self, **kwargs):
        """
        NAME:      plotLMLTandDLDMLT(self, **kwargs)
        USE:          Survey plots of the L and MLT, dL and dMLT of both spacecraft,
                          and the times and duration of the conjunctions. 
        INPUT:      Required: None. Optional: startDate, endDate for making plots
                          cleaner. formatTimes to make the plot time stamps look more 
                          readable. A_id, and B_id is the spacecraft ID's used in labeling
                          the plots. Will assume the mission_id_A/B class instances 
                          otherwise.
        RETURNS: None
        AUTHOR:  Mykhaylo Shumko
        MOD:         2017-05-30
        """  
        startDate = kwargs.get('startDate', None)
        endDate = kwargs.get('endDate', None)
        formatTimes = kwargs.get('formatTimes', False)
        A_id = kwargs.get('A_id', self.mission_id_A)
        B_id = kwargs.get('B_id', self.mission_id_B)
        
        fig = plt.figure(figsize=(15, 10), dpi = 80, facecolor = 'grey')
        gs = gridspec.GridSpec(5,1)
        self.dtPlt = fig.add_subplot(gs[0, 0])
        self.LPlt = fig.add_subplot(gs[1:3, 0], sharex = self.dtPlt)
        self.MLTPlt = fig.add_subplot(gs[3:, 0], sharex = self.dtPlt)
        self.dMLTPlt = self.MLTPlt.twinx()
        
        # Plot duration of conjunctions
        self.dtPlt.scatter(self.durationDict['startTime'],
                           self.durationDict['duration'], marker ='x', 
                            c = 'g', s = 100)
        self.dtPlt.set_yticks(range(1, 5))
       
        # PLOT L SHELL and DL
        self.LPlt.plot(self.magephemA[self.timeKeyA], self.magephemA[self.LkeyA], 
                     '-', c = 'r', label = str(A_id) + ' L')
        self.LPlt.plot(self.magephemB[self.timeKeyB], self.magephemB[self.LkeyB], 
                     '--', c = 'r', label = str(B_id) + ' L')
                     
        self.LPlt.plot(self.magephemA[self.timeKeyA], self.dL, 'k',
        label = r'$\Delta L$')
 
        # PLOT MLT and 
        self.MLTPlt.plot(self.magephemA[self.timeKeyA], 
        self.magephemA[self.MLTkeyA], '-', c = 'b', label = str(A_id) + ' MLT')
        self.MLTPlt.plot(self.magephemB[self.timeKeyB], 
        self.magephemB[self.MLTkeyB], '--', c = 'b', label = str(B_id) + ' MLT')
        self.MLTPlt.plot(self.magephemA[self.timeKeyA], self.dMLT, 'k',
        label = r'$\Delta MLT$')             
        
        # Format Plot
        self.MLTPlt.set_xlabel('UTC')
        self.MLTPlt.set_ylabel('MLT')
        self.LPlt.set_ylabel('McLlwain L (T89, Kp = 2)')
        self.LPlt.legend(loc = 1)
        self.MLTPlt.legend(loc = 1)
        self.LPlt.set_ylim((0, 10))
        self.MLTPlt.set_ylim((0, 25))
        
        if (startDate is not None) and (endDate is not None):
            self.LPltPlt.set_xlim((startDate, endDate))
        
        if formatTimes: # Format time stamps.
            import matplotlib.dates as mdates
            days = mdates.DayLocator()
            dayFmt = mdates.DateFormatter('%D')
            self.MLTPlt.xaxis.set_major_locator(days)
            self.MLTPlt.xaxis.set_major_formatter(dayFmt)
            hour = mdates.HourLocator()
            self.MLTPlt.xaxis.set_minor_locator(hour)
            
        self.LPlt.legend(loc = 1)
        self.MLTPlt.set_xlabel('UTC')
        self.dtPlt.set_ylabel('Conjunction duration (minutes)')
        self.dtPlt.set_title(str(A_id) + '-' + str(B_id) + 
        ' conjunction ' + r' $\Delta L_{thresj} = $' + str(self.L_thresh) +  
        r', $\Delta MLT_{thresh} = $' + str(self.MLT_thresh))
        self.LPlt.set_xlim((self.magephemA[self.timeKeyA][0], 
            self.magephemA[self.timeKeyA][-1]))
        gs.tight_layout(fig)
        plt.show()
        return
        
    def save_to_file(self, **kwargs):
        """
        NAME:      save_to_file(self, **kwargs)
        USE:         Saves the data to a JSON headed ASCII file      
        INPUT:      save_dir and save_name, the directory and name of the saved file.
        RETURNS: None
        AUTHOR: Mykhaylo Shumko
        MOD:     2017-05-30
        """
        # Make these mandatory kwargs when I am done!!!
        save_dir = kwargs.get('save_dir', '/home/ms30715/ssd_data/shumko/' +
        'conjunctions/intermediate_processing')
        save_name = kwargs.get('save_name', None)
        
        globalAttrs = {'deltaL':str(self.L_thresh), 'deltaMLT':str(self.MLT_thresh), 
                         'minL':str(self.lowerL)}
        
        if save_name is None:
            # If spacecraft id's not specified, do not include them in the filename.
            if (self.mission_id_A is not 'missionA') and (self.mission_id_B is not 'missionB'):
                startDate = self.magephemA[self.timeKeyA][0].date()
                endDate = self.magephemA[self.timeKeyA][-1].date()
                save_name = ('{}_{}_{}_dL_{}_dMLT_{}_conjunctions.txt'.format(
                    self.mission_id_A, self.mission_id_B, startDate, self.L_thresh,
                    self.MLT_thresh))
                globalAttrs['missions'] = [self.mission_id_A, self.mission_id_B]
            else:
                save_name = (self.magephemA[self.timeKeyA][0].date().isoformat() + 
                '_' + self.magephemA[self.timeKeyA][-1].date().isoformat() + 
                '_' + self.mission_id_A + '_' + self.mission_id_B + 
                '_conjunctions.txt')
                
        startTimeAttrs = {"UNIT":'ISO 8601 UTC', 'DESCRIPTION':'Start time of the '+
        'conjunction', 'TITLE':'Start time'}
        endTimeAttrs = {"UNIT":'ISO 8601 UTC', 'DESCRIPTION':'End time of the '+
        'conjunction', 'TITLE':'End time'}
        durationAttrs = {'TITLE':'Duration', 'DESCRIPTION':'Duration of the ' +
        'conjunction', 'UNIT':'Minutes' }
        meanLattrs = {'TITLE':'Mean L shell', 'DESCRIPTION':'Mean value of the L shells' + 
        ' during the conjunction'}
        meanMLTattrs = {'TITLE':'Mean MLT', 'DESCRIPTION':'Mean value of the MLT' + 
        ' during the conjunction'}

        outData = dm.SpaceData(attrs = globalAttrs)
        outData['startTime'] = dm.dmarray(self.durationDict['startTime'], 
        attrs = startTimeAttrs)
        outData['endTime'] = dm.dmarray(self.durationDict['endTime'], 
        attrs = endTimeAttrs)
        outData['duration'] = dm.dmarray(self.durationDict['duration'], 
        attrs = durationAttrs)
        
        # Save the mean L and MLT data.
        outData['meanL_' + self.mission_id_A] = dm.dmarray(
            self.durationDict['meanL_' + self.mission_id_A], attrs = meanLattrs)
        outData['meanL_' + self.mission_id_B] = dm.dmarray(
            self.durationDict['meanL_' + self.mission_id_B], attrs = meanLattrs)
        outData['meanMLT_' + self.mission_id_A] = dm.dmarray(
            self.durationDict['meanMLT_' + self.mission_id_A], attrs = meanMLTattrs)
        outData['meanMLT_' + self.mission_id_B] = dm.dmarray(
            self.durationDict['meanMLT_' + self.mission_id_B], attrs = meanMLTattrs)
        
        order = ['startTime', 'endTime', 'duration',
                   'meanL_' + self.mission_id_A, 'meanL_' + self.mission_id_B,
                   'meanMLT_' + self.mission_id_A, 'meanMLT_' + self.mission_id_B]  
        if self.verbose: 
            print('Saving data to {}'.format(os.path.join(save_dir, save_name)))
        dm.toJSONheadedASCII(os.path.join(save_dir, save_name), outData, order = order, 
                             depend0 = 'startTime')             
                             
    def calc_L_MLT_conjunction_delay(self, delay, **kwargs):
        ###############################
        ###### TEST FUNCTION ########
        ###############################
        """
        NAME:      
        USE:          This function calculates conjunctions with the same L and MLT
        thresholds, but allows for a maximum time delay, in minutes.
        INPUT:      
        RETURNS: 
        AUTHOR: Mykhaylo Shumko
        MOD:     2017-05-30
        """
        N = len(self.magephemA[self.timeKeyA])
        Ndelay = delay//self.cadence
        
        # If delay is 0, run the faster code
        if delay == 0:
            self.L_thresh = kwargs.get('Lthresh', self.L_thresh)
            self.MLT_thresh = kwargs.get('MLTthresh', self.MLT_thresh)
            self.ind = np.where((self.dMLT < self.MLT_thresh) & (self.dL < self.L_thresh))[0]
            print(self.ind)
            return self.ind
            
        self.ind = np.array([])
        
        # This is not the most efficiant, but will work!
        for i in range(N):
            # If you are at the ends of the array, the delay 
            # boundaries become asymmetric.
            # IMPLEMENT THIS SOON!
            if i < Ndelay or i > N - Ndelay:
                continue
            # Arrays of conjunction indicies for L and MLT.
            cLind = np.where(np.abs(self.magephemA[self.LkeyA][i] - 
                self.magephemB[self.LkeyB][i - Ndelay:i + Ndelay + 1]) 
                < self.L_thresh)[0]
            cMLTind = np.where(np.abs(self.magephemA[self.MLTkeyA][i] - 
                self.magephemB[self.MLTkeyB][i - Ndelay:i + Ndelay + 1]) 
                < self.MLT_thresh)[0]
            commonInd = np.intersect1d(cLind, cMLTind)
            if len(commonInd) > 0:
                self.ind = np.append(self.ind, i - Ndelay + commonInd[0])
        return
        
    def non_periodic_data_indicies(self, tThresh = 0):
        ###################################
        ######## EXPERIMENTAL ##########
        ###################################
        """
        NAME:      
        USE:   
        INPUT:      
        RETURNS: 
        AUTHOR: Mykhaylo Shumko
        MOD:    2017-05-30
        """
        if self.verbose: print('Reducing the data to the common times.')
        self.matchedIdx = np.array([], dtype = int).reshape(0,2)
        
        # Loop over one of the magephem times
        for i in range(len(self.magephemA[self.timeKeyA])):
            # For each time in the first magephem list, look for a matching 
            # time in the other. Use the threshold if case time stamps are not 
            # at completely the same.
            ltBoolArr = (self.magephemA[self.timeKeyA][i] >= 
                        self.magephemB[self.timeKeyB])
            gtBoolArr = (self.magephemA[self.timeKeyA][i] <= 
                        self.magephemB[self.timeKeyB] + 
                        timedelta(seconds = tThresh))

            idx = np.where(ltBoolArr & gtBoolArr)[0]
            assert len(idx) == 1, ('ERROR: None or multiple times matched!'
                ' (check the tThresh kwarg)')
            indArr = np.array([i, idx[0]]).reshape((1,2))
            self.matchedIdx = np.concatenate((self.matchedIdx, indArr))
        print(self.matchedIdx)

        # Now filter the data
        for key in list(self.magephemA.keys()):
            self.magephemA[key] = self.magephemA[key][self.matchedIdx[:, 0]]
        for key in list(self.magephemB.keys()):
            self.magephemB[key] = self.magephemB[key][self.matchedIdx[:, 1]]
            
        if self.verbose: print('Data reduced for 1-to-1 comparison.')
        return 

    def interpConjunction(self):

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

def locateConsecutiveNumbers(x):
    """
    NAME:    locateConsecutiveNumbers(x)
    USE:       Locates consecutive numbers in an array.
    INPUT:   An sorted array of numbers.
    RETURNS: startIndex and endIndex arrays of consecutive numbers. 
                    These arrays assume that they will be used as slices. Otherwise there 
                    will be a off by 1 index problem! Use caution. 
                    
                    Example:
                    x1 = np.append(np.arange(10), [15, 16])
                    startInd, endInd = locateConsecutiveNumbers(x1)
                    for i in range(len(startInd)):
                        print(x1[startInd[i]:endInd[i]])
                    
                    OUTPUT: 
                    [0 1 2 3 4 5 6 7 8 9], [15 16]
                    
    AUTHOR:  Mykhaylo Shumko
    MOD:     2017-05-19
    """
    # The append is to fix the off by one error!
    x = np.append(x, -9999)
    conv = np.convolve([1, -1], x, mode = 'valid') - 1
    consecutiveFlag = np.where(conv != 0)[0] + 1
    startInd = np.insert(consecutiveFlag, 0, 0)
    endInd = np.insert(consecutiveFlag, len(consecutiveFlag), len(x)-1)
    return startInd[:-1], endInd[:-1] # Fix the off by 1 error!

if __name__ == '__main__':
    fNameA = '20171130_FU4_T89_MagEphem.txt'
    fNameB = 'rbspa_def_MagEphem_T89D_20171130_v1.0.0.txt'
    # fNameA = '20171129_FU4_T89_MagEphem.txt'
    # fNameB = 'rbspa_def_MagEphem_T89D_20171129_v1.0.0.txt'
    fDirA = '/home/mike/research/firebird/Datafiles/FU_4/magephem/'
    fDirB = '/home/mike/research/rbsp/magephem/rbspa/'

    cCalc = MagneticConjunctionCalc(fDirA, fNameA, fDirB, fNameB, Lthresh=1, MLTthresh=1, timeKeyB='DateTime', MLTkeyB='EDMAG_MLT', LkeyB='L', LcolB=-1, lowerL=3)

    cCalc.magephemB[cCalc.timeKeyB] = np.array([i.replace(tzinfo=None)  
                for i in cCalc.magephemB[cCalc.timeKeyB]])
            
    cCalc.periodic_data_indicies() 
    cCalc.calc_magnetic_seperation()
    cCalc.calc_L_MLT_conjunction_delay(0)
    cCalc.lower_L_bound()

    cCalc.calc_conjunction_duration()
    cCalc.plotLMLTandDLDMLT()