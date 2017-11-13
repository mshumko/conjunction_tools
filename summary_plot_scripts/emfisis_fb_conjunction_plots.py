import numpy as np
import sys, os
from datetime import datetime, timedelta
sys.path.insert(0, os.path.join(os.path.dirname( __file__ ), '..', '..', 'auxiliary'))
import plot_emfisis_spectra
import download_day_emfisis
import matplotlib.pylab as plt
from matplotlib import dates
import matplotlib.gridspec as gridspec
import dates_in_filenames
import spacepy.datamodel
import multiprocessing
import dateutil.parser

class FIREBIRD_RBSP_Conjunction_Plots:
    def __init__(self, rbsp_id, fb_id, **kwargs):
        """
        This class loads in the conjunction data between FIREBIRD and RBSP,
        and generates EMFISIS & FIREBIRD data plots for those times, as well
        as 
        """
        self.rbsp_id = rbsp_id
        self.fb_id = fb_id
        self.tPad = kwargs.get('tPad', timedelta(minutes = 2))
        self.temp_data_dir = kwargs.get('temp_save_dir', 
            '/home/ms30715/ssd_data/shumko/data/temp_data')
        self.saveDir = kwargs.get('saveDir', 
            '/home/ms30715/ssd_data/shumko/results/conjunctions/{}'.format(
                datetime.now().date().isoformat()))
        if not os.path.exists(self.saveDir):
            os.makedirs(self.saveDir)
            print('Created directory: {}'.format(self.saveDir))
        return
            
    def readCdata(self, cDataPath = None, cDataName = None, dL = 1, dMLT = 1): 
        """
        This function reads in the conjunction data. 
        """
        if (cDataPath is None) and (cDataName is None):
            cDataPath = ('/home/ms30715/ssd_data/shumko/0_code/conjunctiontoolkit/data/merged_conjunctions/2017_05_31_FB_T89_RBSP_TS04_conjunctions')
            self.cDataName = 'FU{}_RBSP{}_dmlt_{}_dL_{}_conjunctions_hr_filtered.txt'.format(
                self.fb_id, self.rbsp_id, dMLT, dL)
        
        # Read in the conjunction data
        self.cData = spacepy.datamodel.readJSONheadedASCII(
            os.path.join(cDataPath, self.cDataName))
            
        # Convert conjunction times to datetime objects
        for key in ['startTime', 'endTime']:
            self.cData[key] = np.array(list(map(dateutil.parser.parse, self.cData[key])))

        return
        
    def downloadEMFISISData(self, saveDir = None, **kwargs):
        """
        
        """
        detectorType = kwargs.get('detectorType', 'WFR')
        spectraType = kwargs.get('spectraType', 'spectral-matrix-diagonal')
        hour = kwargs.get('hour', '')
        
        if saveDir is None:
            saveDir = self.temp_data_dir
        
        conjunctionDates =  list(set([i.date() for i in self.cData['startTime']]))
        conjunctionDateTimes = [datetime.combine(i, datetime.min.time()) 
            for i in conjunctionDates]
        #print(conjunctionDates)
        for date in conjunctionDateTimes:
            download_day_emfisis.saveEMFISISSpectra(self.rbsp_id, date, saveDir,
                detectorType = detectorType, spectraType = spectraType, hour = hour)
    
        return
        
    def generatePlots(self, saveType = 'png', saveImg = True, saveUsrData = False):
        """
        Runs a loop that generates and saves plots of EMFISIS and RB data.
        """
        if saveUsrData:
            # Prep output data to be appended.
            saveData = spacepy.datamodel.SpaceData()
            observedEvents = np.zeros(len(self.cData['startTime']), dtype = object)
            
            # Copy over the conjunction data.
            for i in list(self.cData.keys()):
                saveData[i] = spacepy.datamodel.dmarray(self.cData[i])
        
        # Loop over the conjunctions.   
        for t in range(len(self.cData['startTime'])):
            print('Processing {}-{}, {} to {} conjunction'.format(self.fb_id, self.rbsp_id,
                self.cData['startTime'][t].isoformat(), self.cData['endTime'][t].isoformat()))
            tBounds = [self.cData['startTime'][t] - self.tPad, 
                self.cData['endTime'][t] + self.tPad]
            
            # Format date correctly
            date = datetime.combine(self.cData['startTime'][t].date(), datetime.min.time())
            statMsg = self.compare_emfisis_fb_col(date, tBounds = tBounds) # Create plots
            if statMsg == -9999:
                continue
            # Save plots
            if saveImg:
                # Format times into <YYYY><MM><DD>T<HH><MM><SS> format.
                startDate =  datetime.now().isoformat().replace('-', '').replace(':', '')[0:15]
                endDate = tBounds[1].isoformat().replace('-', '')
                saveName = '{}_FU{}_RBSP{}_conjunction'.format(
                    startDate, self.fb_id, self.rbsp_id)
                self.saveFig(saveName = saveName, saveType = saveType)
                
            # Get input from user.    
            if saveUsrData:
                plt.show()                
                #observedEvents[t]
                observedEvents[t] = input('Are microbursts, precipitation bands, or chorus, observed? (X,Y,Z) (Y/N) ')
                print(observedEvents[t])

        # Save the user input data.
        if saveUsrData:
            absEAttrs = {'TITLE':'Observed events', 
                'TYPES':'(Microburst, precipitation band, chorus' }
            saveData['obsEvents'] = spacepy.datamodel.dmarray(observedEvents, attrs = absEAttrs)
#            order = ['startTime', 'endTime', 'duration',
#                   'meanL_' + self.mission_id_A, 'meanL_' + self.mission_id_B,
#                   'meanMLT_' + self.mission_id_A, 'meanMLT_' + self.mission_id_B, 
#                   ]  
            
            spacepy.datamodel.toJSONheadedASCII(
                os.path.join(self.saveDir, self.cDataName.replace('.', '_events.')), 
                    saveData, depend0 = 'startTime')

        return
              
    def compare_emfisis_fb_col(self, date, tBounds = None, **kwargs):
        """
        Plots EMFISIS WFR spectra and FIREBIRD Collimated HiRes data.
        """
        hr_dir = kwargs.get('hr_dir', 
            '/home/ms30715/ssd_data/firebird/data/FU_{}/hires/level2'.format(self.fb_id))
        rbsp_magephem_dir = kwargs.get('rbsp_magephem_dir', os.path.join(os.sep, 
            'home', 'ms30715', 'ssd_data', 'magEIS', 'data_mission', 'rbsp_magephem', 
            'rbsp{}'.format(rbsp_id.lower()), 'def', '{}'.format(date.year)))
        yscale = kwargs.get('yscale', 'log')
        show = kwargs.get('show', False)
        plotLoc = kwargs.get('plotLoc', True)
   
        # Find matching HiRes file.
        hrPaths, hrDates = dates_in_filenames.find_dates(hr_dir, dateTimeConvert = True)
        hrMatchIdx = np.where(date == np.array(hrDates))[0]
        assert len(hrMatchIdx) == 1, 'Error, none or multiple hires files found!'
        hr = spacepy.datamodel.readJSONheadedASCII(hrPaths[hrMatchIdx[0]])
        
        with multiprocessing.Pool(processes = 20) as p: # Convert HiRes times
                hr['Time'] = np.array(list(
                    p.map(dateutil.parser.parse, hr['Time'])))
        
        if plotLoc: # Load in RBSP magephem
            rbMagPaths, rbMagDates = dates_in_filenames.find_dates(
                rbsp_magephem_dir , dateTimeConvert = True, wildCard = '*T89Q*.txt')
            rbMagMatchIdx = np.where(date == np.array(rbMagDates))[0]
            assert len(rbMagMatchIdx) == 1, ('Error, none or multiple RBSP ' +
                 'magephem files found!')
        
        try:
            rbspMagEph = spacepy.datamodel.readJSONheadedASCII(
                rbMagPaths[rbMagMatchIdx[0]])
        
            with multiprocessing.Pool(processes = 20) as p: # Convert RBSP MagEphem times
                    rbspMagEph['DateTime'] = np.array(list(
                        p.map(dateutil.parser.parse, rbspMagEph['DateTime'])))
            rbspMagEph['DateTime'] = np.array([i.replace(tzinfo=None)  
                            for i in rbspMagEph['DateTime']]) # Remove timezone
        except ValueError as err:
            print(str(err))
            return -9999
        
        if tBounds is not None:
            hrDataInd = np.where( (hr['Time'] > tBounds[0]) & (hr['Time'] < tBounds[1]))[0]
            magDataInd = np.where( (rbspMagEph['DateTime'] > tBounds[0]) & 
                (rbspMagEph['DateTime'] < tBounds[1]) )[0]
                
        ###########################################
        ############## PLOT STUFFS ###############
        ###########################################
        fig = plt.figure(figsize=(20, 10), dpi=80, facecolor = 'white')
        
        if plotLoc:
            gs = gridspec.GridSpec(6,1)
            ax = fig.add_subplot(gs[0:2, 0], facecolor='k')
            bx = fig.add_subplot(gs[2:4, 0], facecolor='w', sharex = ax)
            LPlt = fig.add_subplot(gs[4:6, 0], facecolor='w', sharex = ax)
            MLTPlt = LPlt.twinx()
            
            for i in [ax, bx, LPlt]:
                i.axvline(tBounds[0] + self.tPad)
                i.axvline(tBounds[1] - self.tPad)                
        else:
            gs = gridspec.GridSpec(4,1)
            ax = fig.add_subplot(gs[0:2, 0], facecolor='k')
            bx = fig.add_subplot(gs[2:4, 0], facecolor='w', sharex = ax)
        bx.set_yscale(yscale)
        
        bx.plot(hr['Time'][hrDataInd], hr['Col_counts'][hrDataInd])
        bx.set(title = 'Collimated FU{} Counts'.format(self.fb_id), ylabel = 'Counts')

        if plotLoc:    
            LPlt.plot(hr['Time'][hrDataInd], np.abs(hr['McllwainL'][hrDataInd]), '-ro', 
                label = 'FU{}'.format(fb_id))
            LPlt.plot(rbspMagEph['DateTime'][magDataInd], 
                rbspMagEph['L'][magDataInd, -2], '--ro', 
                label = 'RBSP{}'.format(rbsp_id)) 
            LPlt.yaxis.label.set_color('r')
            LPlt.spines['left'].set_color('r')
            LPlt.tick_params(axis='y', colors='r')
            LPlt.legend(loc = 2)
            LPlt.set_ylim([1, 10])
                
            MLTPlt.plot(hr['Time'][hrDataInd], hr['MLT'][hrDataInd], '-go', 
                label = 'FU{}'.format(fb_id))
            MLTPlt.plot(rbspMagEph['DateTime'][magDataInd], 
                rbspMagEph['CDMAG_MLT'][magDataInd], '--go', 
                label = 'RBSP{}'.format(rbsp_id))
            MLTPlt.yaxis.label.set_color('g')
            MLTPlt.spines['right'].set_color('g')            
            MLTPlt.tick_params(axis='y', colors='g')                
            MLTPlt.legend(loc = 1)
            
                
            LPlt.set(ylabel = 'L shell (red), solid is FU{}, dashed is RBSP{}'.format(
               self.fb_id, self.rbsp_id))
            MLTPlt.set(ylabel = 'MLT (green), solid is FU{}, dashed is RBSP{}'.format(
               self.fb_id, self.rbsp_id))
        
        #fig.get_axes()
        plot_emfisis_spectra.plotWFRSpectra(self.rbsp_id, date, 
            tBounds = tBounds, cLevels = 0.5, ax = ax, dataDir = self.temp_data_dir )
            
        # Format the time stamps
        fmtr = dates.DateFormatter('%H:%M:%S')
        
        MLTPlt.xaxis.set_major_formatter(fmtr)
        for i in [ax, bx]:
            plt.setp(i.get_xticklabels(), visible = False)  

        gs.tight_layout(fig)
        gs.update(wspace = 0.05, hspace = 0.05) # Set spacing between exes
        
        if show ==  True:
            plt.show()
        return 1
        
    def saveFig(self, **kwargs):
        saveName = kwargs.get('saveName', None)
        saveType = kwargs.get('saveType', 'png')
        print(os.path.join(self.saveDir, saveName + '.' + saveType))
        plt.savefig(os.path.join(self.saveDir, saveName + '.' + saveType), dpi = 80)
        return
    
if __name__ == '__main__':
    for rbsp_id in ['A', 'B']:
        for fb_id in [3, 4]:
            cPlt = FIREBIRD_RBSP_Conjunction_Plots(rbsp_id, fb_id)
            cPlt.readCdata(dL = '0dot1', dMLT = 3)
            #cPlt.downloadEMFISISData()
            cPlt.generatePlots(saveImg = True, saveUsrData = False)
    
    # Load 
#    date = datetime(2015, 2, 2)
#    tBounds = [datetime(2015,2,2,7,44,0), datetime(2015,2,2,8)]
#    compare_emfisis_fb_col(rbsp_id, fb_id, date, tBounds)
#    
