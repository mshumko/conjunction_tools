# This script processes FIREBIRD ephemeris into the future.
import sys
import os
import dateutil.parser
import numpy as np

import spacepy.datamodel

sys.path.insert(0, '/home/mike/research/firebird/data_processing/'         'magnetic_ephemeris')
import make_magnetic_ephemeris

for sc_id in [3, 4]:
    ephemDir = '/home/mike/research/conjunction-tools/2018_04_predicted_ephem'
    ephemName = ('FU{}_2018-04-19_2018-05-14_LLA_ephemeris_def.csv'.format(sc_id))
    saveName = ephemName.split('_')
    saveName[-2] = 'TLE_magephem.txt'
    saveName = '_'.join(saveName)
    print(saveName)

    # Read in the ephemeris
    #ephem = spacepy.datamodel.readJSONheadedASCII(
    #    os.path.join(ephemDir, ephemName))
    #ephem['dateTime'] = np.array([dateutil.parser.parse(i) for i in
    #    ephem['dateTime']])

    # Now run the magnetic ephemeris code
    magEphem = make_magnetic_ephemeris.Magnetic_Ephemeris(
        'FB', ephemDir=ephemDir, ephemName=ephemName, singleKp=20, sc_id=str(sc_id),
        outputName=saveName, outputDir=ephemDir)
    magEphem.run_mag_model(multiprocess=True, nonPeriodicFlag=False,
        mirrorPointAlt=False)
    # magEphem.outputName = saveName
    # magEphem.outputDir = ephemDir
    magEphem.saveToFile(daily=False)
