### This is an exploratory script for the AC6-RBSP conjunction filters ###
import spacepy.datamodel as dm
import multiprocessing
import dateutil.parser
import numpy as np
import os, sys
import time
import matplotlib.pylab as plt
import matplotlib.dates
import matplotlib.gridspec as gridspec

from datetime import datetime, date, timedelta

sys.path.insert(0, '/home/ms30715/research_code/auxiliary')
import dates_in_filenames
import read_ac_data

# Intialize the random module with current system time seed.
import random
random.seed()

def plotExampleEphemeris(ac_id, rb_id, data, idx, plotNum = 100, **kwargs):
    """
    This script will generate, and optionally save plotNum number of plots
    from a AC6-RBSP conjunction dataset, and filtered indicies idx.
    """
    save_dir = kwargs.get('save_dir', None)
    verbose = kwargs.get('verbose', True)
    plotWidth = kwargs.get('plotWidth', 2) # Minutes

    N = len(idx)

    # Load in the list of all magnetic ephemeris files
    # AC6 MagEphem files
    acDir = '/home/ms30715/ssd_data/ac6/ac6{}/ascii/level2'.format(ac_id.lower())
    pathsAC, datesAC = dates_in_filenames.find_dates(acDir, 
        dateTimeConvert = True, wildCard = '*coords*.csv')

    # RBSP MagEphem files
    rbDir = ('/home/ms30715/ssd_data/shumko/data/RBSP/orbit/RBSP{}'
             '/magephem_interpolated/AC6{}'.format(rb_id.upper(), ac_id.upper()))
    pathsRB, datesRB =  dates_in_filenames.find_dates(rbDir, dateTimeConvert = True, 
        wildCard = '*TS04D*.txt')

    # Convert to serial time
    datesAC = matplotlib.dates.date2num(datesAC)
    datesRB = matplotlib.dates.date2num(datesRB)

    while plotNum != 0:
        # Pick a random filtered index, and get its date
        i = random.choice(idx)
        startDate = matplotlib.dates.date2num(data['startTime'][i].date())

        ### Find the files to plot ###
        acEphemIdx = np.where(datesAC == startDate)[0]
        rbEphemIdx = np.where(datesRB == startDate)[0]
        print(rbEphemIdx)
        assert (len(acEphemIdx) == 1 and len(rbEphemIdx) == 1), ('ERROR:'
            ' None or multiple matches found!')

        fDirAC = pathsAC[acEphemIdx[0]]
        fDirRB = pathsRB[rbEphemIdx[0]]
        if verbose: 
            print('Loading files: {} and {}.'.format(os.path.basename(fDirAC),
                os.path.basename(fDirRB)))

        ### LOAD IN DATA ###
        rbMagEphem = dm.readJSONheadedASCII(fDirRB)
        # Convert RBSP times
        tConvert = lambda x: dateutil.parser.parse(x).replace(tzinfo=None)
        rbMagEphem['DateTime'] = np.array(list(map(tConvert, 
            rbMagEphem['DateTime'])))
        acMagEphem = read_ac_data.read_ac_data(fDirAC)
        
        ### PLOT DATA ###
        fig = plt.figure(figsize = (10, 8), dpi = 80)
        gs = gridspec.GridSpec(2, 1)
        lPlt = plt.subplot(gs[0, 0])
        mltPlt = plt.subplot(gs[1, 0], sharex = lPlt)

        # Boolean time arrays to filter conjunction times
        acInd = (acMagEphem['dateTime'] > data['startTime'][i] - timedelta(
            minutes = plotWidth)) & (acMagEphem['dateTime'] < data['endTime'][i]
            + timedelta(minutes = plotWidth))
        rbInd = (rbMagEphem['DateTime'] > data['startTime'][i] - timedelta(
            minutes = plotWidth)) & (rbMagEphem['DateTime'] < data['endTime'][i]
            + timedelta(minutes = plotWidth))

        acTimeIdx = np.where(acInd)[0]
        rbTimeIdx = np.where(rbInd)[0]

        # Plot L shells
        lPlt.plot(acMagEphem['dateTime'][acInd], acMagEphem['Lm_OPQ'][acInd], 
            label = 'AC6-{} L'.format(ac_id.upper()))
        lPlt.plot(rbMagEphem['DateTime'][rbInd], rbMagEphem['Lstar'][rbInd, -2], 
            label = 'RBSP-{} L'.format(rb_id.upper()))
        lPlt.set(title = 'AC-{} - RBSP-{} Conjunction \n {}'.format(ac_id, rb_id,
            data['startTime'][i].date().isoformat()), ylabel = 'L')
        lPlt.legend(loc = 1)
        lPlt.set_ylim(top = 10, bottom = 1)

        # Plot MLT
        mltPlt.plot(acMagEphem['dateTime'][acInd], acMagEphem['MLT_OPQ'][acInd], 
            label = 'AC6-{} MLT'.format(ac_id.upper()))
        mltPlt.plot(rbMagEphem['DateTime'][rbInd], 
            rbMagEphem['EDMAG_MLT'][rbInd], 
            label = 'RBSP-{} MLT'.format(rb_id.upper()))
        mltPlt.set(ylabel = 'MLT', xlabel = 'UTC')
        mltPlt.legend(loc = 1)

        # Plot vertical lines at the start/end of conjunction.
        lPlt.axvline(data['startTime'][i]); lPlt.axvline(data['endTime'][i])
        mltPlt.axvline(data['startTime'][i]); mltPlt.axvline(data['endTime'][i])
        
        # Save the data if the flag is set.
        if save_dir is None:
            plt.show()
        else:   
            save_name = 'AC{}_RBSP{}_{}_conjunction_magepgem.png'.format(
                ac_id, rb_id, data['startTime'][i].strftime('%Y%m%d_%H%M%SUT'))  
            plt.savefig(os.path.join(save_dir, save_name), dpi = 200)   
        plotNum -= 1
    return

def dMLT(A, B):
    """
    This function finds the signed difference between MLT arrays A and B. This
    accounts for the periodicity of MLT over midnight.

    Negative output means A is west of B. 
    """
    # Find the straight difference.
    dMLT = A - B
    # Apply boolean mask if difference is greater than 12. 
    # This accounts for the periodicity acorss midnight... most cases
    dMLT[dMLT > 12] = -((B[dMLT > 12] - A[dMLT > 12]) % 24)
    # Address cases when A near mightnight and B right after midnight)
    dMLT[dMLT < -12] = (-(B[dMLT < -12] - A[dMLT < -12]) % 24)

##    ### OLD CODE THAT WORKS ###
##    dMLT[abs(dMLT) > 12] = -(B[abs(dMLT) > 12] - A[abs(dMLT) > 12]) % 24
##    # Address cases when A near mightnight and B right after midnight)
##    dMLT[abs(dMLT) > 12] = -((B[abs(dMLT) > 12] - A[abs(dMLT) > 12]) % 24)
    return dMLT
### TEST FOR dMLT function.
##MLTA = np.arange(24)
##MLTB = np.arange(24, 0, -1)
##dMLT = dMLT(MLTA, MLTB)
##for i in range(len(MLTA)):
##    print(MLTA[i], MLTB[i], dMLT[i])#, MLTA[i] -  MLTB[i])

startTime = time.time()
ac_id = 'A'
rb_id = 'A'
fDir = '/home/ms30715/ssd_data/shumko/0_code/conjunctiontoolkit/data/merged_conjunctions/2017_06_13_AC6_OPQ_RBSP_TS04_conjunctions'
fName = 'AC{}_RBSP{}_dmlt_3_dL_0dot25_conjunctions.txt'.format(ac_id, rb_id)
Lthresh = 0.01

# Load data and convert datetimes
data = dm.readJSONheadedASCII(os.path.join(fDir, fName))

for key in ['startTime', 'endTime']:
    with multiprocessing.Pool() as p:
        data[key] = np.array(list(p.map(dateutil.parser.parse, data[key])))

# Create a index array  
idx = np.arange(len(data['startTime']))
print('Original database length', len(idx))

# Indicies where dL < Lthresh
idxL = np.where(np.abs(data['meanL_AC6{}'.format(ac_id)] - 
           data['meanL_RBSP{}'.format(rb_id)]) < Lthresh)[0]

# Indicies where dMLT > 0 (AC6 is ahead in MLT)
idxMLTdrift = np.where(dMLT(data['meanMLT_AC6{}'.format(ac_id)], 
           data['meanMLT_RBSP{}'.format(rb_id)]) > 0)[0]

# Use set notation to get the intersection of the L
idx = list(set(idx) & set(idxL) & set(idxMLTdrift))
print('Filtered database length', len(idx))

print('Running random conjunction plotter')
plotExampleEphemeris(ac_id, rb_id, data, idx, save_dir = '/home/ms30715/ssd_data/shumko/temp/conjunction_validation')

print('Execution time: {}'.format(time.time() - startTime))


