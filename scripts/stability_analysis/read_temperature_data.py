'''
This script reads stored temperature data from the B-60 lab
'''

import numpy as np
import os
from matplotlib import pyplot as plt


datapath = r'D:\measuring\data\environment'

# dict of sensor ID and location
sensors = {'F300000783E0F228': 'Inside airco at M1 setup', '3C0000078361CD28': 'Short test lead', 'DA00000090419A26-H': 'Humidity, next to whiteboard & fuse box', '8D0000078403E328': 'On Table next to Montana M1 (starting 201607 08:09:14)', '0E000007834A8728': 'Taped to magnet Y-stage (starting 20160711 15:00:00)', '2D00000782775728': 'Above table rack hanging from cable gutter', 'C400000783BCC028': 'DS18B20 Inside air inlet above small table', '58000007839D8428': 'Inside airco at other half of room', 'E60000078284E828': 'Hanging from cable gutter in other half of room'}

for id in sensors:
    if id[-2:] != '-H':
        print id[-2:]
        print id
    # data = readTemperatureData(id.append('.txt'))

def readTemperatureData(datafile):
    fid = open(datapath + '\\' + datafile + '.txt', 'r')
    fid.readline()
    