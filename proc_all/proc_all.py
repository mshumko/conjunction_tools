# This script processes all of the RBSP-FIREBIRD conjunctions

import os
import sys
import glob

sys.path.insert(0, '/home/mike/research/conjunction-tools/')
import conjunctiontools_v2

fb_dir = lambda sc: '/home/mike/research/firebird/Datafiles/FU_{}/ephem/'.format(sc)
rbsp_dir = '/home/mike/research/conjunction-tools/proc_all/rbsp_camp_ephem/'

for fb_id in [3]:
    for rb_id in ['a']:
        for camp in range(1, 11):
            # Find matching files
            fbFiles = glob.glob(os.path.join(fb_dir(fb_id), 'FU{}*camp{:02d}*'.format(fb_id, camp)))
            assert len(fbFiles) == 1, '0 or multiple FIREBIRD ephem files found!'

            rbFiles = glob.glob(os.path.join(rbsp_dir, 'rbsp{}*camp{:02d}*'.format(rb_id, camp)))
            assert len(rbFiles) == 1, '0 or multiple RBSP ephem files found!'

            # Map file to magnetic coordinates.