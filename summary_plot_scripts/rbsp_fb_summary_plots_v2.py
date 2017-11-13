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
#        self.temp_data_dir = kwargs.get('temp_save_dir', 
#            '/home/ms30715/ssd_data/shumko/data/temp_data')
        self.saveDir = kwargs.get('saveDir', 
            '/home/ms30715/ssd_data/shumko/results/conjunctions/{}'.format(
                datetime.now().date().isoformat()))
        if not os.path.exists(self.saveDir):
            os.makedirs(self.saveDir)
            print('Created directory: {}'.format(self.saveDir))
        return
            
    def readConjunctiondata(self, cDataPath=None, cDataName=None, dL=1, dMLT=1): 
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
        
    def generatePlots(self, saveType='png', saveImg=True):
        """
        Runs a loop that generates and saves plots of RBSP
        and RB data.
        """
        
        # Loop over the conjunctions.   
        for t in range(len(self.cData['startTime'])):
            print('Processing {}-{}, {} to {} conjunction'.format(
                self.fb_id, self.rbsp_id,
                self.cData['startTime'][t].isoformat(), self.cData['endTime'][t].isoformat()))
            # Widen the time bounds a little.
            tBounds = [self.cData['startTime'][t] - self.tPad, 
                self.cData['endTime'][t] + self.tPad]
            
            # Plot RBSP
            self.plotMagEIS(self.rbsp_id, tBounds, ax[0:3])
            self.plotEMFISIS(self.sc_id, tBounds, ax[3])
            
            # Plot FIREBIRD
            
            # Plot Position
            self.plotPosition(ax[-2:])
            
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
                self.saveFig(saveName=saveName, saveType=saveType)
        return
        
    def plotMagEIS(self, sc_id, tBounds, axArr):
        """
        This function is a wrapper to plot rel03 MagEIS data
        for the 0, 90, and 180 degree pitch angles. The axArr
        is an array of three subplots to which plot the
        different pitch angle data.
        """
        rel03Obj = PlotMageis(rb_id, tBounds[0], dtype,
            tRange=tBounds)
        rel03Obj.loadMagEphem(Bmodel='T89D')
        self.RBSPmagEphem = rel03obj.magEphem
        rel03Obj.plotUnidirectionalFlux(0, ax=axArr[0])
        rel03Obj.plotUnidirectionalFlux(90, ax=axArr[1])
        rel03Obj.plotUnidirectionalFlux(180, ax=axArr[2])
        return
        
    def plotEMFISIS(self, sc_id, tBounds, ax):
        """
        Wrapper to plot EMFISIS data.
        """
        pObj = EMFISISspectra(sc_id, tBounds[0],
            tBounds=tBounds)
        pObj.loadWFRSpectra()
        pObj.loadMagEphem()
        pObj.plotSpectra(ax=ax)
        return 
        
    def plotPosition(self, axArr):
        # Plot RBSP magnetic ephemeris
        axArr[0].plot(self.RBSPmagEphem['DateTime'], 
            self.magEphem['Lstar'][:, 0], label='RBSP')
        axArr[1].plot(self.RBSPmagEphem['DateTime'], 
            self.magEphem['EDMAG_MLT'][:], label='RBSP')  
            
        # Plot FIREBIRD magnetic ephemeris  
        
        # Make plot readable
        axArr[0].set(ylim=(3, 10), ylabel='L')
        axArr[1].set(ylabel='MLT')
        return
        
    def saveFig(self, **kwargs):
        saveName = kwargs.get('saveName', None)
        saveType = kwargs.get('saveType', 'png')
        print(os.path.join(self.saveDir, saveName + '.' + saveType))
        plt.savefig(os.path.join(self.saveDir, saveName + '.' + saveType), dpi = 80)
        return
    
if __name__ == '__main__':
    for rbsp_id in ['A', 'B']:
        for fb_id in [3, 4]:
            cPlt = FIREBIRD_RBSP_Conjunction_Plots(
                rbsp_id, fb_id)
            cPlt.readConjunctiondata(dL='0dot1', dMLT=3)
            cPlt.generatePlots(saveImg=True)
