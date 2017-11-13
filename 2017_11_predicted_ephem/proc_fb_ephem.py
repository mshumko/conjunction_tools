# This script processes FIREBIRD ephemeris into the future.
import sys
import os
import dateutil.parser
import numpy as np

import spacepy.datamodel

sys.path.insert(0, '/home/mike/research/firebird/data_processing/'         'magnetic_ephemeris')
import make_magnetic_ephemeris

sc_id = 4
ephemDir = '/home/mike/research/mission-tools/orbit/data'
ephemName = ('FU{}_SGP4_LLA_2017-11-13_to_2017-12-22'                        '_gen_w_2017-11-13_TLE.txt'.format(sc_id))
saveName = ephemName.split('_')
saveName[-1] = 'TLE_magepehem.txt'
saveName = '_'.join(saveName)
print(saveName)

# Read in the ephemeris
ephem = spacepy.datamodel.readJSONheadedASCII(
    os.path.join(ephemDir, ephemName))
ephem['dateTime'] = np.array([dateutil.parser.parse(i) for i in
    ephem['dateTime']])

# Now run the magnetic ephemeris code
magEphem = make_magnetic_ephemeris.Magnetic_Ephemeris(
    'FB', ephem=ephem, singleKp=20, sc_id=str(sc_id),
    outputName=ephemName, outputDir=ephemDir)
magEphem.run_mag_model(multiprocess=True, nonPeriodicFlag=False,
    mirrorPointAlt=False)
# magEphem.outputName = saveName
# magEphem.outputDir = ephemDir
magEphem.saveToFile(daily=False)