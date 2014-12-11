'''script to analyze multiple ssro runs '''

import sys
sys.path.append(r'D:/measuring')
from analysis.lib.m2.ssro import ssro
reload(ssro)
from analysis.lib.tools import toolbox


timestamp = '20141122_235000'

while 1:
	timestamp, folder = toolbox.latest_data('AdwinSSRO', older_than = timestamp, return_timestamp = True)

	ssro.ssrocalib(folder = folder) 