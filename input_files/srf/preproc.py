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

parser.add_argument('--min_abs_value',
                    default = 0.,
                    dest = 'min_abs_value',
                    metavar = '<minumum absolute value>',
                    type = float,
                    help = '')

parser.add_argument('--missing',
                    default = None,
                    dest = 'missing',
                    metavar = '<minimum wavelength> <maximum wavelength>',
                    nargs = 3,
                    type = float,
                    help = '')

parser.add_argument('--nm_to_um',
                    action = 'store_true',
                    default = False,
                    dest = 'nm_to_um',
                    help = '')

parser.add_argument('--reverse',
                    action = 'store_true',
                    default = False,
                    dest = 'reverse',
                    help = '')

parser.add_argument('--skip_rows',
                    default = 0,
                    dest = 'skip_rows',
                    metavar = '<rows to skip>',
                    type = int,
                    help = '')

parser.add_argument('--use_columns',
                    default = (0,1),
                    dest = 'use_columns',
                    metavar = ('<wavelength column>', '<SR column>'),
                    nargs = 2,
                    type = int,
                    help = '')

args = parser.parse_args()

if args.missing == None:
    data = np.loadtxt(args.input_file, skiprows=args.skip_rows, usecols=args.use_columns)
else:
    data = np.empty((int(args.missing[0]), 2))
    data[0,0] = args.missing[1]
    data[1,0] = args.missing[2]
    data[:,1] = 1.

if args.nm_to_um:
    scale = 1.e-3
else:
    scale = 1.

if args.reverse:
    data = data[::-1,:]

f = open(args.output_file, 'w')

for i in range(0, data.shape[0]):
    if abs(data[i,1]) >= args.min_abs_value:
         f.write("%0.6e % .5e\n" % (data[i,0] * scale, data[i,1]))

f.close()
