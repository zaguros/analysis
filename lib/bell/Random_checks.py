import numpy as np


def get_list_block_sizes(Random_num):
	"""
	Returns a list with the block sizes in the random string.
	"""

	block_size = 0
	blocks = []

	for i, j in enumerate(Random_num):
		if i < (len(Random_num) - 1):
			if j != Random_num[i+1]:
				block_size = block_size + 1
				blocks.append(block_size)
				block_size = 0
			else:
				block_size = block_size + 1
		else:
			if j != Random_num[i-1]:
				blocks.append(1)
			else:
				block_size = block_size + 1
				blocks.append(block_size)

	return blocks

def autocorr(Random_num):
	"""
	Returns the autocorrelation of a string of Random numbers
	"""

	result = np.correlate(Random_num, Random_num, mode='full')
	print len(result)

	return result, result[result.size/2:]