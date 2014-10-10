
from analysis.lib.magnetometry import adaptive_magnetometry as magnetometry
reload(magnetometry)

def adaptive_table(tau0=1e-9):

	print '******* Generating ADAPTIVE TABLE *********'
	print
	for n in [6]:
		for m in [3]:
			
			print '##### N = '+str(n)+' --- M = '+str(m)
			
			t = magnetometry.AdaptiveTable (N=n,M=m)
			ttt = int(tau0*1e9)
			t.set_tau0(tau0=ttt*1e-9)
			t.save_folder = '/home/cristian/Work/Research/teamdiamond/measurement/scripts/Magnetometry/adaptive_tables_lt1/tau0='+str(ttt)+'ns/'
			#t.save_folder = 'D:/measuring/measurement/scripts/Magnetometry/adaptive_tables_lt1/tau0='+str(ttt)+'ns/'
			t.verbose = False
			t.generate()
			t.save_table()
	print 
	print '**** DONE!'


def test_table ():
	t = magnetometry.AdaptiveTable (N=6,M=2)
	t.set_tau0(tau0=1e-9)
	t.verbose = True
	t.msmnt_to_position (msmnt_results = [0,2,2,2,1,0])


adaptive_table(tau0=1e-9)
#test_table()
