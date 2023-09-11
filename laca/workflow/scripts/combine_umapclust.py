#!/usr/bin/env python

import os
import sys

with open(sys.argv[1], 'w') as outfile:
    outfile.write('read\tlength\tD1\tD2\tbin_id\tbatch_id\n')
    for file in os.listdir(sys.argv[2]):
        if file.endswith('.tsv'):
            batch_id = file.split('.')[0].split('_')[-1].replace('batch', 'b')
            file_path = os.path.join(sys.argv[2], file)
            with open(file_path) as infile:
                next(infile)
                for line in infile:
                    outfile.write(line.rstrip() + '\t' + batch_id + '\n')
