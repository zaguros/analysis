with open('.git/info/attributes','w') as f:
	f.write('*.ipynb\tfilter=dropoutput_ipynb\n')

with open('.git/config','a') as f:
	f.write('[filter "dropoutput_ipynb"]\n')
	f.write('\tclean = python git_ipynb_output_filter.py\n')
	f.write('\tsmudge = cat\n')



	
	

