import numpy as np
import sys, os
from datetime import datetime, timedelta
import matplotlib.pylab as plt
import glob
from matplotlib import dates
import matplotlib.gridspec as gridspec
import multiprocessing
import dateutil.parser

import spacepy.datamodel

sys.path.insert(0, '/home/mike/research/mission-tools/rbsp/')
import plot_emfisis_spectra
import plot_mageis

class FIREBIRD_RBSP_Conjunction_Plots:
    def __init__(self, rbsp_id, fb_id, **kwargs):
        """
        This class loads in the conjunction data between FIREBIRD and RBSP,
        and generates EMFISIS & FIREBIRD data plots for those times, as well
        as 
        """
        self.rbsp_id = rbsp_id
        self.fb_id = fb_id
        self.plot_empty_data = kwargs.get('plot_empty_data', True)
        self.tPad = kwargs.get('tPad', timedelta(minutes = 2))
        self.fb_dir = ('/home/mike/research/firebird/'
            'Datafiles/FU_{}/hires/level2'.format(self.fb_id))
        self.saveDir = kwargs.get('saveDir', 
            '/home/mike/research/conjunction-tools/plots/{}'.format(
                datetime.now().date().isoformat()))
        if not os.path.exists(self.saveDir):
            os.makedirs(self.saveDir)
            print('Created directory: {}'.format(self.saveDir))
        return
            
    def readConjunctionData(self, path, dL=1, dMLT=1): 
        """
        This function reads in the conjunction data. 
        """
        # Read in the conjunction data
        self.cData = spacepy.datamodel.readJSONheadedASCII(path)
            
        # Convert conjunction times to datetime objects
        for key in ['startTime', 'endTime']:
            self.cData[key] = np.array(list(map(dateutil.parser.parse, self.cData[key])))
        return
        
    def generatePlots(self, saveType='png', saveImg=True):
        """
        Runs a loop that generates and saves plots of RBSP
        and RB data.
        """
        # Create subplots
        fig = plt.figure(figsize=(8, 11.5), dpi=80, facecolor = 'white')
        gs = gridspec.GridSpec(6,1)
        ax = [None]*6
        for i in range(len(ax)):
            ax[i] = fig.add_subplot(gs[i, 0])
        
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
            self.plotEMFISIS(self.rbsp_id, tBounds, ax[3])
            
            # Plot FIREBIRD
            self.plotFIREBIRD(self.fb_id, tBounds, ax[4])

            # Plot Position
            self.plotPosition(ax[-1])

            # If FIREBIRD data is empty and self.plot_empty_data
            # kwarg is false, skip the plot
            if ((not self.plot_empty_data) and 
                    (len(self.validFBIdt) == 0)):
                continue

            # Beautify plots
            datetimeStr = self.cData['startTime'][t].isoformat(
                ).replace('-', '').replace(':', '')[0:15]
            for a in ax[1:]:
                a.set_title('')
            ax[0].set_ylabel(r'$J(\alpha_L = 0^\circ)$')
            ax[1].set_ylabel(r'$J(\alpha_L = 90^\circ)$')
            ax[2].set_ylabel(r'$J(\alpha_L = 180^\circ)$')
            ax[0].set(title='FU{} RBSP{} conjunction {}'.format(
                self.fb_id, self.rbsp_id, self.cData['startTime'][t]))

            if saveImg: # Save plots
                saveName = '{}_FU{}_RBSP{}_conjunction'.format(
                    datetimeStr, self.fb_id, self.rbsp_id)
                self.saveFig(saveName=saveName, saveType=saveType)
            for a in ax: # Clear subplots
                a.cla()
        return
        
    def plotMagEIS(self, sc_id, tBounds, axArr):
        """
        This function is a wrapper to plot rel03 MagEIS data
        for the 0, 90, and 180 degree pitch angles. The axArr
        is an array of three subplots to which plot the
        different pitch angle data.
        """
        rel03Obj = plot_mageis.PlotMageis(rb_id, tBounds[0], 'rel03',
            tRange=tBounds)
        rel03Obj.loadMagEphem(Bmodel='T89D')
        self.RBSPmagEphem = rel03Obj.magEphem
        rel03Obj.plotUnidirectionalFlux(0, ax=axArr[0],
            pltLegendLoc=False)
        rel03Obj.plotUnidirectionalFlux(90, ax=axArr[1],
            pltLegendLoc=False)
        rel03Obj.plotUnidirectionalFlux(180, ax=axArr[2],
            pltLegendLoc=False)
        return
        
    def plotEMFISIS(self, sc_id, tBounds, zx):
        """
        Wrapper to plot EMFISIS data.
        """
        pObj = plot_emfisis_spectra.EMFISISspectra(sc_id, tBounds[0],
            tBounds=tBounds)
        pObj.loadWFRSpectra()
        pObj.loadMagEphem()
        pObj.plotSpectra(ax=zx, plotCb=False, grid=False, 
            legendLoc=False)
        return 

    def plotFIREBIRD(self, fb_id, tBounds, zx):
        fb_path = 'FU{}_Hires_{}_L2.txt'.format(
            self.fb_id, tBounds[0].date())
        self.hr = spacepy.datamodel.readJSONheadedASCII(
            os.path.join(self.fb_dir, fb_path))
        # Convert timestamps
        with multiprocessing.Pool() as p: # Convert HiRes times
            self.hr['Time'] = np.array(list(
                p.map(dateutil.parser.parse, self.hr['Time'])))
        # fix timestamps
        meanDt = np.mean(self.hr['Count_Time_Correction'])
        self.hr['Time'] = np.array([i + timedelta(seconds=meanDt) for i in self.hr['Time']])

        self.validFBIdt = np.where((self.hr['Time'] > tBounds[0]) & 
            (self.hr['Time'] < tBounds[1]))[0]
        zx.plot(self.hr['Time'][self.validFBIdt], 
            self.hr['Col_flux'][self.validFBIdt], 
            label=self.hr['Col_flux'].attrs['ELEMENT_LABELS'])
        return

    def plotPosition(self, zx):
        # Plot and annotate RBSP magnetic ephemeris
        zx.plot(self.RBSPmagEphem['DateTime'], 
            self.RBSPmagEphem['Lstar'][:, 0], 
            label='RBSP{}'.format(self.rbsp_id))
        zx.text(x=0.01, y=0.95, s='RBSP\nMLT={}\nMLAT={}'.format(
            round(np.mean(self.RBSPmagEphem['EDMAG_MLT'][:])),
            round(np.mean(self.RBSPmagEphem['EDMAG_MLAT'][:]))),
            transform=zx.transAxes, va='top')
            
        # Plot and annotate FIREBIRD magnetic ephemeris if data exists
        if len(self.validFBIdt) > 0:
            zx.plot(self.hr['McIlwainL'][self.validFBIdt], 
                label='FU{}'.format(self.fb_id))
            zx.text(x=.01, y=0, s='FB\nMLT={}\nLAT={}\nLON={}'.format(
                round(np.mean(self.hr['MLT'][self.validFBIdt])),
                round(np.mean(self.hr['Lat'][self.validFBIdt])),
                round(np.mean(self.hr['Lon'][self.validFBIdt]))
                ), transform=zx.transAxes, va='bottom')
        
        # Make plot readable
        zx.set(ylim=(3, 10), ylabel='L')
        zx.legend(loc=1)
        return
        
    def saveFig(self, **kwargs):
        saveName = kwargs.get('saveName', None)
        saveType = kwargs.get('saveType', 'png')
        print(os.path.join(self.saveDir, saveName + '.' + saveType))
        plt.savefig(os.path.join(self.saveDir, saveName + '.' + saveType), dpi = 80)
        return
    
if __name__ == '__main__':
    rb_id = 'A'
    fb_id = 3
    CONJUNCTION_DIR = ('/home/mike/research/conjunction-tools/data/'     'merged_conjunctions/2017-10-10_RBSP_FB_T89_conjunctions')
    paths = glob.glob('{}/FU{}_RBSP{}*'.format(CONJUNCTION_DIR, fb_id, rb_id))
    assert len(paths) == 1, 'None or multiple conjunction files found!'

    # Run summary plot generator.
    cPlt = FIREBIRD_RBSP_Conjunction_Plots(
        rb_id, fb_id)
    cPlt.readConjunctionData(paths[0])
    cPlt.generatePlots(saveImg=True)