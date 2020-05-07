#! /usr/bin/env python
import screed
import sys
import os

for filename in sys.argv[1:]:
    print(filename)
    with open(filename + '.nostop.fa', 'wt') as fp:
        for record in screed.open(filename):
            record.sequence = record.sequence.replace('*', '')
            fp.write('>{}\n{}\n'.format(record.name, record.sequence))

        
