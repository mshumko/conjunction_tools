# This script maps RBSP ephem to magephem using IRBEM

import numpy as np
import csv
import glob
import os

import IRBEM

inDir = '/home/mike/research/conjunction-tools/proc_all/rbsp_camp_ephem'
outDir = '/home/mike/research/conjunction-tools/proc_all/rbsp_camp_magephem'
overwrite = False # Flag to overwrite any magephem files already generated
inPaths = sorted(glob.glob(os.path.join(inDir, '*')))

# Now get a list of outFiles (magephem files to save to).
inBasename = [os.path.basename(f) for f in inPaths]
outBasename = [x.replace('.', '_magephem.') for x in inBasename]
outPaths = [os.path.join(outDir, x) for x in outBasename]

# Now loop over every file, and generate the corresponding
# magephem file. If the file already exists, check is overwrite
# flag is set.
for inF, outF in zip(inPaths, outPaths):
    if not overwrite and not os.path.isfile(outF):
        # Read the ephem file

        # Generate magephem paramaters

        # Write variables to outF.