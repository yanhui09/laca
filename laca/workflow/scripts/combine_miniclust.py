#!/usr/bin/env python

import os
import sys

with open(sys.argv[1], 'w') as outfile:
    outfile.write('bin_id,cluster,read_id,read_len,clust_read_score,frac_max_score\n')
    for file in os.listdir(sys.argv[2]):
        if file.endswith('.csv'):
            batch_id = file.split('.')[0]
            file_path = os.path.join(sys.argv[2], file)
            with open(file_path) as infile:
                next(infile)
                for line in infile:
                    # add batch_id before the second ","
                    line_new = ','.join(line.split(',', 2)[:2]) + batch_id + ',' + ','.join(line.split(',', 2)[2:])
                    outfile.write(line_new)
