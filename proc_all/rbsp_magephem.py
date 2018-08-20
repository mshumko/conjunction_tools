# This script maps RBSP ephem to magephem using IRBEM

import numpy as np
import csv
import glob
import os
import sys
import IRBEM

sys.path.insert(0, '/home/mike/research/firebird/data_processing/magnetic_ephemeris')
import make_magnetic_ephemeris

inDir = '/home/mike/research/conjunction-tools/proc_all/rbsp_camp_ephem'
outDir = '/home/mike/research/conjunction-tools/proc_all/rbsp_camp_magephem'
# If none, will look for existing kp in data. Change to a number when forward propagating!
singleKp = 20

overwrite = False # Flag to overwrite any magephem files already generated
inPaths = sorted(glob.glob(os.path.join(inDir, '*.txt')))

# Now get a list of outFiles (magephem files to save to).
inBasename = [os.path.basename(f) for f in inPaths]
outBasename = [x.replace('.', '_magephem.') for x in inBasename]
outPaths = [os.path.join(outDir, x) for x in outBasename]

# Now loop over every file, and generate the corresponding
# magephem file. If the file already exists, check is overwrite
# flag is set.
for inF, outF in zip(inPaths, outPaths):
    if not overwrite and not os.path.isfile(outF):
        # Process the magnetic ephemeris.
        print('Processing {} -> {}'.format(os.path.basename(inF), os.path.basename(outF)))
        magEphem = make_magnetic_ephemeris.Magnetic_Ephemeris(
            'rbsp', ephemDir=os.path.dirname(inF), ephemName=os.path.basename(inF), 
             outputDir=os.path.dirname(outF), outputName=os.path.basename(outF),
             singleKp=singleKp)
        magEphem.run_mag_model(multiprocess=True, nonPeriodicFlag=False,
            mirrorPointAlt=False)
        magEphem.outputName = os.path.basename(outF)
        magEphem.outputDir = os.path.dirname(outF)
        magEphem.saveToFile(daily=False)