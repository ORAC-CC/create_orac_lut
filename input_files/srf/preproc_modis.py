#! /usr/bin/env python

import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'For nerf herding')

parser.add_argument(dest = 'input_file',
                    metavar = '<input filename>',
                    type = str,
                    help = 'Input filename')

parser.add_argument(dest = 'output_file',
                    metavar = '<output filename>',
                    type = str,
                    help = 'Output filename')

parser.add_argument('--nm_to_um',
                    action = 'store_true',
                    default = False,
                    dest = 'nm_to_um',
                    help = '')

parser.add_argument('--skip_rows',
                    default = 0,
                    dest = 'skip_rows',
                    metavar = '<rows to skip>',
                    type = int,
                    help = '')

args = parser.parse_args()

data = np.loadtxt(args.input_file, skiprows=args.skip_rows)

if args.nm_to_um:
    scale = 1.e-3
else:
    scale = 1.

f = open(args.output_file, 'w')

index = np.where((data[:,1] == 1) & (data[:,3] != -99.))
for i in index[0]:
    f.write("%.6e %.5e\n" % (data[i,2] * scale, data[i,3]))

f.close()
