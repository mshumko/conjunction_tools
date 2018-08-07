# This script processes FIREBIRD ephemeris into the future.
import sys
import os
import dateutil.parser
import numpy as np

import spacepy.datamodel

sys.path.insert(0, '/home/mike/research/firebird/data_processing/'         'magnetic_ephemeris')
import make_magnetic_ephemeris

for sc_id in ['A', 'B']:
    ephemDir = ('/home/mike/research/conjunction-tools/'
        '2018_04_predicted_ephem')
    ephemName = ('rbsp{}_ephem_def.txt'.format(
        sc_id.lower()))
    saveName = ephemName.split('_')
    saveName[-2] = 'magephem'
    saveName = '_'.join(saveName)
    print(saveName)

    # Now run the magnetic ephemeris code
    magEphem = make_magnetic_ephemeris.Magnetic_Ephemeris(
        'rbsp', singleKp=20, sc_id=str(sc_id),
        ephemName=ephemName, ephemDir=ephemDir,
        outputName=ephemName, outputDir=ephemDir)
    magEphem.run_mag_model(multiprocess=True, nonPeriodicFlag=False,
        mirrorPointAlt=False)
    magEphem.outputName = saveName
    magEphem.outputDir = ephemDir
    magEphem.saveToFile(daily=False)
