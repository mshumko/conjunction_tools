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
            self.cData[key] = np.array(
                list(map(dateutil.parser.parse, self.cData[key])))
        return
        
    def generatePlots(self, saveType='png', saveImg=True, lowMagEISFlux=1E2):
        """
        Runs a loop that generates and saves plots of RBSP
        and RB data.
        """
        # Create subplots
        fig = plt.figure(figsize=(8, 11.5), dpi=80, facecolor = 'white')
        gs = gridspec.GridSpec(6,1, left=0.1, bottom=0.04, 
            right=0.99, top=0.97, wspace=0.1, hspace=0.06)
        ax = [None]*6
        ax[0] = fig.add_subplot(gs[0,0])
        for i in range(1, len(ax)):
            ax[i] = fig.add_subplot(gs[i, 0], sharex=ax[0])
        
        # Loop over the conjunctions.   
        for (st, et) in zip(self.cData['startTime'], self.cData['endTime']):
            print('Processing {}-{}, {} to {} conjunction'.format(
                self.fb_id, self.rbsp_id,
                st.isoformat(), et.isoformat()))
            # Widen the time bounds a little.
            tBounds = [st - self.tPad, 
                et + self.tPad]
            
            # Plot RBSP
            self.plotMagEIS(self.rbsp_id, tBounds, ax[0:3], lowFlux=lowMagEISFlux)
            self.plotEMFISIS(self.rbsp_id, tBounds, ax[3])
            
            # Plot FIREBIRD
            self.plotFIREBIRD(self.fb_id, tBounds, ax[4])

            # Plot Position
            self.plotPosition(ax[-1])

            #If FIREBIRD data is empty and self.plot_empty_data
            #kwarg is false, skip the plot
            if ((not self.plot_empty_data) and 
                    (len(self.validFBIdt) == 0)):
                for a in ax: # Clear subplots
                    a.cla()
                continue

            # Beautify plots
            datetimeStr = st.isoformat().replace('-', '').replace(':', '')[0:15]
            for a in ax[1:]:
                a.set_title('')
            for a in ax[:-1]:
                a.set_xlabel('')
            ax[-1].set_xlabel('UTC')
            
            ax[0].set_ylabel(r'MagEIS $J \ (\alpha_L \approx 0^\circ)$' + '\n(probably not LC)')
            ax[1].set_ylabel(r'MagEIS $J \ (\alpha_L = 90^\circ)$')
            ax[2].set_ylabel(r'MagEIS $J \ (\alpha_L \approx 180^\circ)$' + '\n(probably not LC)')
            ax[3].set_ylabel('EMFISIS WFR (Hz)')
            ax[0].set(title='FU{} RBSP{} conjunction {}'.format(
                self.fb_id, self.rbsp_id, st))
            ax[0].set_xlim(tBounds)

            ax[0].set_ylim((1E2, 1E5))
            ax[1].set_ylim((1E2, 1E6))
            ax[2].set_ylim((1E2, 1E5))
            #for a in ax[0:3]:
            #    a.set_ylim(bottom=lowMagEISFlux)

            for a in ax[:-1]:
                plt.setp(a.get_xticklabels(), visible=False)
            #gs.tight_layout(fig)

            if saveImg: # Save plots
                saveName = '{}_FU{}_RBSP{}_conjunction'.format(
                    datetimeStr, self.fb_id, self.rbsp_id)
                self.saveFig(saveName=saveName, saveType=saveType)
            for a in ax: # Clear subplots
                a.cla()
        return
        
    def plotMagEIS(self, sc_id, tBounds, axArr, lowFlux=False):
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
        # # Shrink legend fontsize
        # for label in rel03Obj.fluxLegend.get_texts():
        #     label.set_fontsize('5')
        rel03Obj.plotUnidirectionalFlux(90, ax=axArr[1],
            pltLegendLoc=False)
        rel03Obj.plotUnidirectionalFlux(180, ax=axArr[2],
            pltLegendLoc=False)
        
        # # Set lower limit flux level
        # if lowFlux:
        #     for aa in axArr:
        #         aa.set_ylim(bottom=lowFlux)
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
        try:
            self.hr = spacepy.datamodel.readJSONheadedASCII(
                os.path.join(self.fb_dir, fb_path))
        except FileNotFoundError as err: # If no file exists
            self.validFBIdt = []
            return

        # Convert timestamps
        with multiprocessing.Pool() as p: # Convert HiRes times
            self.hr['Time'] = np.array(list(
                p.map(dateutil.parser.parse, self.hr['Time'])))
        # fix timestamps
        meanDt = np.mean(self.hr['Count_Time_Correction'])
        self.hr['Time'] = np.array([i + timedelta(seconds=meanDt) 
            for i in self.hr['Time']])

        self.validFBIdt = np.where((self.hr['Time'] > tBounds[0]) & 
            (self.hr['Time'] < tBounds[1]))[0]
        for ee in range(6):
            zx.plot(self.hr['Time'][self.validFBIdt], 
                self.hr['Col_counts'][self.validFBIdt, ee], 
                label=self.hr['Col_counts'].attrs['ELEMENT_LABELS'][ee]
                )
        zx.legend(loc=1, fontsize=8)
        zx.set(yscale='log', ylabel='FIREBIRD\nCol counts')
        return

    def plotPosition(self, zx):
        # Plot and annotate RBSP magnetic ephemeris
        zx.plot(self.RBSPmagEphem['DateTime'], 
            self.RBSPmagEphem['Lstar'][:, 0], 
            label='RBSP{}'.format(self.rbsp_id))
        zx.text(x=0.01, y=0.95, s='RBSP\nMLT={}\nMLAT={}'.format(
            round(np.mean(self.RBSPmagEphem['EDMAG_MLT'][:])),
            round(np.mean(self.RBSPmagEphem['EDMAG_MLAT'][:]))),
            transform=zx.transAxes, va='top', fontsize=8)
            
        # Plot and annotate FIREBIRD magnetic ephemeris if data exists
        if len(self.validFBIdt) > 0:
            zx.plot(self.hr['Time'][self.validFBIdt],
                self.hr['McIlwainL'][self.validFBIdt], 
                label='FU{}'.format(self.fb_id))
            zx.text(x=0.1, y=0.95, 
                s='FB\nMLT={}\nLAT={}\nLON={}'.format(
                round(np.mean(self.hr['MLT'][self.validFBIdt])),
                round(np.mean(self.hr['Lat'][self.validFBIdt])),
                round(np.mean(self.hr['Lon'][self.validFBIdt]))
                ), transform=zx.transAxes, va='top', fontsize=8)
        
        # Make plot readable
        zx.set(ylim=(3, 10), ylabel='L')
        zx.legend(loc=1, fontsize=8)
        return
        
    def saveFig(self, **kwargs):
        saveName = kwargs.get('saveName', None)
        saveType = kwargs.get('saveType', 'png')
        print('Saving to:', os.path.join(
            self.saveDir, saveName + '.' + saveType))
        plt.savefig(os.path.join(self.saveDir, saveName + '.' + saveType), dpi = 80)
        return
    
if __name__ == '__main__':
    CONJUNCTION_DIR = ('/home/mike/research/conjunction-tools/data/'     
        'merged_conjunctions/camp13_RBSP_FB')
    for rb_id in ['A', 'B']:
        for fb_id in [3, 4]:
            paths = glob.glob('{}/FU{}_RBSP{}*'.format(CONJUNCTION_DIR, fb_id, rb_id))
            assert len(paths) == 1, 'None or multiple conjunction files found!'

            # Run summary plot generator.
            cPlt = FIREBIRD_RBSP_Conjunction_Plots(
                rb_id, fb_id, plot_empty_data=True)
            cPlt.readConjunctionData(paths[0])
            cPlt.generatePlots(saveImg=True)