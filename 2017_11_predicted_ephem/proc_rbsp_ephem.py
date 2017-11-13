# This script processes FIREBIRD ephemeris into the future.
import sys
import os
import dateutil.parser
import numpy as np

import spacepy.datamodel

sys.path.insert(0, '/home/mike/research/firebird/data_processing/'         'magnetic_ephemeris')
import make_magnetic_ephemeris

sc_id = 'B'
ephemDir = ('/home/mike/research/conjunction-tools/'
    '2017_11_predicted_ephem')
ephemName = ('RBSP{}_2017-11-13_2017-12-23_ephemeris.txt'.format(
    sc_id))
saveName = ephemName.split('_')
saveName[-1] = 'magephem.txt'
saveName = '_'.join(saveName)
print(saveName)

# Read in the ephemeris
# ephem = spacepy.datamodel.readJSONheadedASCII(
#     os.path.join(ephemDir, ephemName))
# ephem['dateTime'] = np.array([dateutil.parser.parse(i) for i in
#     ephem['dateTime']])

# Now run the magnetic ephemeris code
magEphem = make_magnetic_ephemeris.Magnetic_Ephemeris(
    'rbsp', singleKp=20, sc_id=str(sc_id),
    ephemName=ephemName, ephemDir=ephemDir,
    outputName=ephemName, outputDir=saveName)
magEphem.run_mag_model(multiprocess=True, nonPeriodicFlag=False,
    mirrorPointAlt=False)
# magEphem.outputName = saveName
# magEphem.outputDir = ephemDir
magEphem.saveToFile(daily=False)