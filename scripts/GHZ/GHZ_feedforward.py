import numpy as np

def GHZ_rho():
	GHZ_111 =np.array([
			[ 	(1./2.)	,0	,0	,0	,0	,0	,0	,(1./2.)],
			[	0		,0	,0	,0	,0	,0	,0	,0],
			[	0		,0	,0	,0	,0	,0	,0	,0],
			[	0		,0	,0	,0	,0	,0	,0	,0],
			[	0		,0	,0	,0	,0	,0	,0	,0],
			[	0		,0	,0	,0	,0	,0	,0	,0],
			[	0		,0	,0	,0	,0	,0	,0	,0],
			[	(1./2.)	,0	,0	,0	,0	,0	,0	,(1./2.)]])

	GHZ_101 =np.array([
			[ 	0	,0	,0		,0	,0	,0		,0	,0],
			[	0	,0	,0		,0	,0	,0		,0	,0],
			[	0	,0	,(1./2.)	,0	,0	,-(1./2.)	,0	,0],
			[	0	,0	,0		,0	,0	,0		,0	,0],
			[	0	,0	,0		,0	,0	,0		,0	,0],
			[	0	,0	,-(1./2.)	,0	,0	,(1./2.)	,0	,0],
			[	0	,0	,0		,0	,0	,0		,0	,0],
			[	0	,0	,0		,0	,0	,0		,0	,0]])

	GHZ_110 =np.array([	
			[	0	,0		,0	,0	,0	,0	,0		,0],
			[	0	,1./2.	,0	,0	,0	,0	,-(1./2.)	,0],
			[	0	,0		,0	,0	,0	,0	,0		,0],
			[	0	,0		,0	,0	,0	,0	,0		,0],
			[	0	,0		,0	,0	,0	,0	,0		,0],
			[	0	,0		,0	,0	,0	,0	,0		,0],
			[	0	,-(1./2.)	,0	,0	,0	,0	,1./2.	,0],
			[	0	,0		,0	,0	,0	,0	,0		,0]])

	GHZ_100 =np.array([
				[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,1./2.	,1./2.	,0	,0	,0],
				[0	,0	,0	,1./2.	,1./2.	,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0]])

	GHZ_011 =np.array([
				[0	,0	,0	,0		,0		,0	,0	,0],
			 	[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,1./2.	,-(1./2.)	,0	,0	,0],
				[0	,0	,0	,-(1./2.)	,1./2.	,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0],
				[0	,0	,0	,0		,0		,0	,0	,0]])

	GHZ_001 =np.array([
				[0	,0		,0	,0	,0	,0	,0		,0],
				[0	,(1./2.)	,0	,0	,0	,0	,(1./2.)	,0],
				[0	,0		,0	,0	,0	,0	,0		,0],
				[0	,0		,0	,0	,0	,0	,0		,0],
				[0	,0		,0	,0	,0	,0	,0		,0],
				[0	,0		,0	,0	,0	,0	,0		,0],
				[0	,(1./2.)	,0	,0	,0	,0	,(1./2.)	,0],
				[0	,0		,0	,0	,0	,0	,0		,0]])

	GHZ_010 =np.array([
			 	[0	,0	,0		,0	,0	,0		,0	,0],
				[0	,0	,0		,0	,0	,0		,0	,0],
				[0	,0	,(1./2.)	,0	,0		,(1./2.)	,0	,0],
				[0	,0	,0		,0	,0	,0		,0	,0],
				[0	,0	,0		,0	,0	,0		,0	,0],
				[0	,0	,(1./2.)	,0	,0		,(1./2.)	,0	,0],
				[0	,0	,0		,0	,0	,0		,0	,0],
				[0	,0	,0		,0	,0	,0		,0	,0]])

	GHZ_000 =np.array([
				[1./2.		,0	,0	,0	,0	,0	,0	,-(1./2.)],
				[	0		,0	,0	,0	,0	,0	,0	,0],
				[	0		,0	,0	,0	,0	,0	,0	,0],
				[	0		,0	,0	,0	,0	,0	,0	,0],
				[	0		,0	,0	,0	,0	,0	,0	,0],
				[	0		,0	,0	,0	,0	,0	,0	,0],
				[	0		,0	,0	,0	,0	,0	,0	,0],
				[	-(1./2.)	,0	,0	,0	,0	,0	,0	,1./2.]])


	return(GHZ_000,GHZ_001,GHZ_010,GHZ_011,GHZ_100,GHZ_101,GHZ_110,GHZ_111)

def GHZ_bra():
	GHZ_111 = 1./np.sqrt(2)*np.array([1,0,0,0,0,0,0,1])
	GHZ_000 = 1./np.sqrt(2)*np.array([1,0,0,0,0,0,0,-1])
	GHZ_001 = 1./np.sqrt(2)*np.array([0,1,0,0,0,0,1,0])
	GHZ_110 = 1./np.sqrt(2)*np.array([0,1,0,0,0,0,-1,0])
	GHZ_010 = 1./np.sqrt(2)*np.array([0,0,1,0,0,1,0,0])
	GHZ_101 = 1./np.sqrt(2)*np.array([0,0,1,0,0,-1,0,0])
	GHZ_100 = 1./np.sqrt(2)*np.array([0,0,0,1,1,0,0,0]) 
	GHZ_011 = 1./np.sqrt(2)*np.array([0,0,0,1,-1,0,0,0])
	return(GHZ_000,GHZ_001,GHZ_010,GHZ_011,GHZ_100,GHZ_101,GHZ_110,GHZ_111)

def GHZ_order():
	return(np.array(['000','001','010','011','100','101','110','111']))

def sigmax():
	return(np.array([[0, 1],[1, 0]]))
def sigmay():
	return(np.array([[0, -1j],[1j, 0]]))
def sigmaz():
	return(np.array([[1, 0],[0, -1]]))

def sx1():
	return(tensor(sigmax(), np.eye(2), np.eye(2)))
def sy1():
    return(tensor(sigmay(), np.eye(2), np.eye(2)))
def sz1():
    return(tensor(sigmaz(), np.eye(2), np.eye(2)))    
def sx2():
    return(tensor(np.eye(2), sigmax(), np.eye(2)))
def sy2():
	return(tensor(np.eye(2), sigmay(), np.eye(2)))
def sz2():
	return(tensor(np.eye(2), sigmaz(), np.eye(2)))
def sx3():
	return(tensor(np.eye(2), np.eye(2), sigmax()))
def sy3():
	return(tensor(np.eye(2), np.eye(2), sigmay()))
def sz3():
	return(tensor(np.eye(2), np.eye(2), sigmaz()))


    # Calculate products
def oper(q1,q2,q3):
	if q1 == 'x':
		m1 = sx1()
	elif q1 == 'y':
		m1 = sy1()
	elif q1 == 'z':
		m1 = sz1()
	elif q1 == 'i':
		m1 = np.eye(8)
	if q2 == 'x':
		m2 = sx2()
	elif q2 == 'y':
		m2 = sy2()
	elif q2 == 'z':
		m2 = sz2()
	elif q2 == 'i':
		m2 = np.eye(8)
	if q3 == 'x':
		m3 = sx3()
	elif q3 == 'y':
		m3 = sy3()
	elif q3 == 'z':
		m3 = sz3()
	elif q3 == 'i':
		m3 = np.eye(8)
	return(np.dot(m1,np.dot(m2,m3)))

def xyyop():
	return(np.dot(sx1(),np.dot(sy2(),sy3())))
def yxyop():
	return(np.dot(sy1(),np.dot(sx2(),sy3())))
def yyxop():
	return(np.dot(sy1(),np.dot(sy2(),sx3())))
def xxxop():
	return(np.dot(sx1(),np.dot(sx2(),sx3())))
def yyyop():
	return(np.dot(sy1(),np.dot(sy2(),sy3())))	

def tensor(in1,in2,in3):
	return(np.kron(in1,np.kron(in2,in3)))

def calc_exp_value(state,operator):
	e=np.dot(state,np.dot(operator,np.transpose(state)))
	return e

def GHZ_expectation_values(operator):
	GHZ_states = GHZ_bra()
	exp=np.zeros(GHZ_order().shape)
	for k,ghz in enumerate(GHZ_states):
		exp[k]=calc_exp_value(ghz,operator)
	return(exp)

def all_GHZ_expectations():
	bases = np.array(['x','y','z','i'])
	names = GHZ_order()
	table = np.append(np.array([['-']]),np.transpose(names))
	print(table)
	for i1,base1 in enumerate(bases):
		for i2,base2 in enumerate(bases):
			for i3,base3 in enumerate(bases):
				operator = oper(base1,base2,base3)
				exp_values = GHZ_expectation_values(operator)
				table_input = np.append(np.array([base1+base2+base3]),exp_values)
				table = np.vstack((table,table_input))
	for row in table:
		printArray([str(x) for x in row])

def nonzero_expectations():
	bases = np.array(['x','y','z','i'])
	names = GHZ_order()
	table = np.append(np.array([['-']]),np.transpose(names))
	for i1,base1 in enumerate(bases):
		for i2,base2 in enumerate(bases):
			for i3,base3 in enumerate(bases):
				if base1+base2+base3 in ['xxx','xyy','yxy','yyx','zzi','ziz','izz']:
					operator = oper(base1,base2,base3)
					exp_values = GHZ_expectation_values(operator)
					table_input = np.append(np.array([base1+base2+base3]),exp_values)
					table = np.vstack((table,table_input))
	for row in table:
		printArray([str(x) for x in row])
	# return(table)


def printArray(args):
    print "\t".join(args)
