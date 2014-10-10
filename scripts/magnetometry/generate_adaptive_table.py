
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
reload(magnetometry)

def adaptive_table(tau0=1e-9):

	print '******* Generating ADAPTIVE TABLE *********'
	print
	for n in [5]:
		for m in [1]:
			
			print '##### N = '+str(n)+' --- M = '+str(m)
			
			t = magnetometry.AdaptiveTable (N=n,M=m)
			ttt = int(tau0*1e9)
			t.set_tau0(tau0=ttt*1e-9)
			t.save_folder = 'D:/measuring/measurement/scripts/Magnetometry/adaptive_tables_lt1/tau0='+str(ttt)+'ns_new/'
			t.verbose = True
			t.generate()
			#t.save_table()
	print 
	print '**** DONE!'


adaptive_table(tau0=1e-9)