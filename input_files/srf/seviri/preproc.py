#! /usr/bin/env python

import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'For nerf herding')

parser.add_argument('--min_abs_value',
                    default = 0.,
                    dest = 'min_abs_value',
                    metavar = '<minumum absolute value>',
                    type = float,
                    help = '')

args = parser.parse_args()


def write_srf_chan(i_chan):
    data = np.loadtxt('../seviri.org/seviri-chan_%d.srf' % i_chan)
    write_srf_msg_chan(1, i_chan, data)
    write_srf_msg_chan(2, i_chan, data)
    write_srf_msg_chan(3, i_chan, data)
    write_srf_msg_chan(4, i_chan, data)


def write_srf_msg_chan(i_msg, i_chan, data):
    f = open('seviri-msg%d-chan_%d.srf' % (i_msg, i_chan), 'w')

    for i in range(0, data.shape[0]):
        if abs(data[i,1]) >= args.min_abs_value:
             f.write("%0.6e % .5e\n" % (data[i,0], data[i,i_msg]))

    f.close()


def write_srf_chan2(i_chan):
    data = np.loadtxt('../seviri.org/seviri-chan_%d.srf' % i_chan)
    write_srf_msg_chan2(1, i_chan, 1, 95, data)
    write_srf_msg_chan2(1, i_chan, 2, 85, data)
    write_srf_msg_chan2(2, i_chan, 1, 95, data)
    write_srf_msg_chan2(2, i_chan, 2, 85, data)
    write_srf_msg_chan2(3, i_chan, 1, 95, data)
    write_srf_msg_chan2(3, i_chan, 2, 85, data)
    write_srf_msg_chan2(4, i_chan, 1, 95, data)
    write_srf_msg_chan2(4, i_chan, 2, 85, data)


def write_srf_msg_chan2(i_msg, i_chan, i_temp, temp, data):
    f = open('seviri-msg%d-chan_%d-temp_%d.srf' % (i_msg, i_chan, temp), 'w')

    for i in range(0, data.shape[0]):
        if abs(data[i,1]) >= args.min_abs_value:
             f.write("%0.6e % .5e\n" % (data[i,0], data[i,(i_msg-1) * 2 + i_temp]))

    f.close()


write_srf_chan(1)
write_srf_chan(2)
write_srf_chan(3)

write_srf_chan2(4)
write_srf_chan2(5)
write_srf_chan2(6)
write_srf_chan2(7)
write_srf_chan2(8)
write_srf_chan2(9)
write_srf_chan2(10)
write_srf_chan2(11)


def write_srf_chan3(i_chan):
    data = np.loadtxt('../seviri.org/seviri-chan_%d.srf' % i_chan)
    write_srf_msg_chan3(1, i_chan, '',          data, 1)
    write_srf_msg_chan3(1, i_chan, '-extended', data, 2)
    write_srf_msg_chan3(2, i_chan, '',          data, 3)
    write_srf_msg_chan3(2, i_chan, '-extended', data, 4)
    write_srf_msg_chan3(3, i_chan, '',          data, 5)
    write_srf_msg_chan3(4, i_chan, '',          data, 6)
    write_srf_msg_chan3(4, i_chan, '-extended', data, 7)


def write_srf_msg_chan3(i_msg, i_chan, tag, data, i_col):
    f = open('seviri-msg%d-chan_%d%s.srf' % (i_msg, i_chan, tag), 'w')

    for i in range(0, data.shape[0]):
        if abs(data[i,1]) >= args.min_abs_value:
             f.write("%0.6e % .5e\n" % (data[i,0], data[i,i_col]))

    f.close()

write_srf_chan3(12)
