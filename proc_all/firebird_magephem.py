# This script maps FIREBIRD ephem to magephem using IRBEM
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
singleKp = None 

overwrite = False # Flag to overwrite any magephem files already generated
inPaths = sorted(glob.glob(os.path.join(inDir, '*')))