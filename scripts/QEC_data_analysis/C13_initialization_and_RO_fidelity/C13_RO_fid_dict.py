import numpy as np

# single-qubit expectation values
Z_C1 = [0.83166, 0.00278]
Z_C2 = [0.8541, 0.00311]
Z_C5 = [0.85399, 0.00312]

# nitrogen initialization fidelity
MBI_N = [0.94125, 0.0283]

# define correction using init = RO
C_C1 = [0,0]
C_C1[0] = np.sqrt(Z_C1[0]/MBI_N[0])
C_C1[1] = C_C1[0]*np.sqrt((MBI_N[1]/(2*MBI_N[0]))**2+(Z_C1[1]/(2*Z_C1[0]))**2)
C_C2 = [0,0]
C_C2[0] = np.sqrt(Z_C2[0]/MBI_N[0])
C_C2[1] = C_C2[0]*np.sqrt((MBI_N[1]/(2*MBI_N[0]))**2+(Z_C2[1]/(2*Z_C2[0]))**2)
C_C5 = [0,0]
C_C5[0] = np.sqrt(Z_C5[0]/MBI_N[0])
C_C5[1] = C_C5[0]*np.sqrt((MBI_N[1]/(2*MBI_N[0]))**2+(Z_C5[1]/(2*Z_C5[0]))**2)

ZZZ = [0.70088, 0.02024]

ZC1_ZZZ = [0.83088, 0.017492]
ZC2_ZZZ = [0.83088, 0.017637]
ZC5_ZZZ = [0.83314, 0.017609]

C_C1C2C5 = [0,0]
C_C1C2C5[0] = ZZZ[0]*np.sqrt(MBI_N[0]*Z_C1[0]*Z_C2[0]*Z_C5[0])/(ZC1_ZZZ[0]*ZC2_ZZZ[0]*ZC5_ZZZ[0])
C_C1C2C5[1] = C_C1C2C5[0]*np.sqrt((ZZZ[1]/ZZZ[0])**2+(MBI_N[1]/(2*MBI_N[0]))**2+(Z_C1[1]/(2*Z_C1[0]))**2+(Z_C2[1]/(2*Z_C2[0]))**2+(Z_C5[1]/(2*Z_C5[0]))**2
							+(ZC1_ZZZ[1]/ZC1_ZZZ[0])**2+(ZC2_ZZZ[1]/ZC2_ZZZ[0])**2+(ZC5_ZZZ[1]/ZC5_ZZZ[0])**2)

#two-qubit in three-qubit tomo
Z1Z2_ZZZ = [0.76419, 0.018880]
Z1Z5_ZZZ = [0.75062, 0.019209]
Z2Z5_ZZZ = [0.71671, 0.019957]


C_C2C5 = [0,0]
C_C2C5[0] = Z2Z5_ZZZ[0]*np.sqrt(Z_C5[0]*Z_C2[0])/(ZC5_ZZZ[0]*ZC2_ZZZ[0])
C_C2C5[1] = C_C2C5[0]*np.sqrt((Z2Z5_ZZZ[1]/Z2Z5_ZZZ[0])**2+(Z_C5[1]/(2*Z_C5[0]))**2+(Z_C2[1]/(2*Z_C2[0]))**2
							+(ZC5_ZZZ[1]/ZC5_ZZZ[0])**2+(ZC2_ZZZ[1]/ZC2_ZZZ[0])**2)

#two-qubit in two-qubit tomo
Z1Z2 = [0.74631, 0.01375]
ZC1_Z1Z2 = [0.82793, 0.012527]
ZC2_Z1Z2 = [0.79886, 0.013060]

C_C1C2 = [0,0]
C_C1C2[0] = 1/2.*(Z1Z2_ZZZ[0]*np.sqrt(Z_C1[0]*Z_C2[0])/(ZC1_ZZZ[0]*ZC2_ZZZ[0])
				+Z1Z2[0]*np.sqrt(Z_C1[0]*Z_C2[0])/(ZC1_Z1Z2[0]*ZC2_Z1Z2[0]))
C_C1C2[1] = 1/2.*np.sqrt((Z_C1[1]/(2*Z_C1[0]))**2+(Z_C2[1]/(2*Z_C2[0]))**2
				+(Z1Z2[1]*np.sqrt(Z_C1[0]*Z_C2[0])/(ZC1_Z1Z2[0]*ZC2_Z1Z2[0]))**2
				+(Z1Z2_ZZZ[1]*np.sqrt(Z_C1[0]*Z_C2[0])/(ZC1_ZZZ[0]*ZC2_ZZZ[0]))**2
				+(ZC1_Z1Z2[1]*Z1Z2[0]*np.sqrt(Z_C1[0]*Z_C2[0])/(ZC1_Z1Z2[0]**2*ZC2_Z1Z2[0]))**2
				+(ZC2_Z1Z2[1]*Z1Z2[0]*np.sqrt(Z_C1[0]*Z_C2[0])/(ZC1_Z1Z2[0]*ZC2_Z1Z2[0]**2))**2
				+(ZC1_ZZZ[1]*Z1Z2_ZZZ[0]*np.sqrt(Z_C1[0]*Z_C2[0])/(ZC1_ZZZ[0]**2*ZC2_ZZZ[0]))**2
				+(ZC2_ZZZ[1]*Z1Z2_ZZZ[0]*np.sqrt(Z_C1[0]*Z_C2[0])/(ZC1_ZZZ[0]*ZC2_ZZZ[0]**2))**2				
				)

Z1Z5 = [0.76308, 0.01352]
ZC1_Z1Z5 = [0.83576, 0.012493]
ZC5_Z1Z5 = [0.86147, 0.012493]

C_C1C5 = [0,0]
C_C1C5[0] = 1/2.*(Z1Z5_ZZZ[0]*np.sqrt(Z_C1[0]*Z_C5[0])/(ZC1_ZZZ[0]*ZC5_ZZZ[0])
				+Z1Z5[0]*np.sqrt(Z_C1[0]*Z_C5[0])/(ZC1_Z1Z5[0]*ZC5_Z1Z5[0]))
C_C1C5[1] = 1/2.*np.sqrt((Z_C1[1]/(2*Z_C1[0]))**2+(Z_C5[1]/(2*Z_C5[0]))**2
				+(Z1Z5[1]*np.sqrt(Z_C1[0]*Z_C5[0])/(ZC1_Z1Z5[0]*ZC5_Z1Z5[0]))**2
				+(Z1Z5_ZZZ[1]*np.sqrt(Z_C1[0]*Z_C5[0])/(ZC1_ZZZ[0]*ZC5_ZZZ[0]))**2
				+(ZC1_Z1Z5[1]*Z1Z5[0]*np.sqrt(Z_C1[0]*Z_C5[0])/(ZC1_Z1Z5[0]**2*ZC5_Z1Z5[0]))**2
				+(ZC5_Z1Z5[1]*Z1Z5[0]*np.sqrt(Z_C1[0]*Z_C5[0])/(ZC1_Z1Z5[0]*ZC5_Z1Z5[0]**2))**2
				+(ZC1_ZZZ[1]*Z1Z5_ZZZ[0]*np.sqrt(Z_C1[0]*Z_C5[0])/(ZC1_ZZZ[0]**2*ZC5_ZZZ[0]))**2
				+(ZC5_ZZZ[1]*Z1Z5_ZZZ[0]*np.sqrt(Z_C1[0]*Z_C5[0])/(ZC1_ZZZ[0]*ZC5_ZZZ[0]**2))**2				
				)


def get_C13_correction(order = [1,5,2], errorbar = None):
	if order == [1,5,2]: #all experiments and 000 tomo
	# RO order = [1,5,2]

		Tomo_correction_list = [
		C_C1[0],C_C1[0],C_C1[0],
		C_C5[0],C_C5[0],C_C5[0],
		C_C2[0],C_C2[0],C_C2[0],

		C_C1C5[0],C_C1C5[0],C_C1C5[0],
		C_C1C5[0],C_C1C5[0],C_C1C5[0],
		C_C1C5[0],C_C1C5[0],C_C1C5[0],

		C_C1C2[0],C_C1C2[0],C_C1C2[0],
		C_C1C2[0],C_C1C2[0],C_C1C2[0],
		C_C1C2[0],C_C1C2[0],C_C1C2[0],

		C_C2C5[0],C_C2C5[0],C_C2C5[0],
		C_C2C5[0],C_C2C5[0],C_C2C5[0],
		C_C2C5[0],C_C2C5[0],C_C2C5[0],

		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],

		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],

		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0]]

		if errorbar == None:
			error_Tomo_correction_list = np.zeros(len(Tomo_correction_list))
		else: 
			error_Tomo_correction_list = [
			C_C1[1],C_C1[1],C_C1[1],
			C_C5[1],C_C5[1],C_C5[1],
			C_C2[1],C_C2[1],C_C2[1],

			C_C1C5[1],C_C1C5[1],C_C1C5[1],
			C_C1C5[1],C_C1C5[1],C_C1C5[1],
			C_C1C5[1],C_C1C5[1],C_C1C5[1],

			C_C1C2[1],C_C1C2[1],C_C1C2[1],
			C_C1C2[1],C_C1C2[1],C_C1C2[1],
			C_C1C2[1],C_C1C2[1],C_C1C2[1],

			C_C2C5[1],C_C2C5[1],C_C2C5[1],
			C_C2C5[1],C_C2C5[1],C_C2C5[1],
			C_C2C5[1],C_C2C5[1],C_C2C5[1],

			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],

			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],

			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1]]

		return Tomo_correction_list, error_Tomo_correction_list

	elif order == [2,5,1]: #all other tomo's
	# RO order = [1,5,2]
		Tomo_correction_list = [
		C_C2[0],C_C2[0],C_C2[0],
		C_C5[0],C_C5[0],C_C5[0],
		C_C1[0],C_C1[0],C_C1[0],

		C_C2C5[0],C_C2C5[0],C_C2C5[0],
		C_C2C5[0],C_C2C5[0],C_C2C5[0],
		C_C2C5[0],C_C2C5[0],C_C2C5[0],

		C_C1C2[0],C_C1C2[0],C_C1C2[0],
		C_C1C2[0],C_C1C2[0],C_C1C2[0],
		C_C1C2[0],C_C1C2[0],C_C1C2[0],

		C_C1C5[0],C_C1C5[0],C_C1C5[0],
		C_C1C5[0],C_C1C5[0],C_C1C5[0],
		C_C1C5[0],C_C1C5[0],C_C1C5[0],

		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],

		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],

		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0],
		C_C1C2C5[0],C_C1C2C5[0],C_C1C2C5[0]]

		if errorbar == None:
			error_Tomo_correction_list = np.zeros(len(Tomo_correction_list))
		else: 
			error_Tomo_correction_list = [
			C_C2[1],C_C2[1],C_C2[1],
			C_C5[1],C_C5[1],C_C5[1],
			C_C1[1],C_C1[1],C_C1[1],

			C_C2C5[1],C_C2C5[1],C_C2C5[1],
			C_C2C5[1],C_C2C5[1],C_C2C5[1],
			C_C2C5[1],C_C2C5[1],C_C2C5[1],

			C_C1C2[1],C_C1C2[1],C_C1C2[1],
			C_C1C2[1],C_C1C2[1],C_C1C2[1],
			C_C1C2[1],C_C1C2[1],C_C1C2[1],

			C_C1C5[1],C_C1C5[1],C_C1C5[1],
			C_C1C5[1],C_C1C5[1],C_C1C5[1],
			C_C1C5[1],C_C1C5[1],C_C1C5[1],

			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],

			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],

			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1],
			C_C1C2C5[1],C_C1C2C5[1],C_C1C2C5[1]]

		return Tomo_correction_list, error_Tomo_correction_list

def get_C13_correction_state(order = [1,5,2], state = 'Z', errorbar = None):
	
	Tomo_correction_list, error_Tomo_correction_list = get_C13_correction(order = order)

	Tomo = {}
	Tomo['X_list'] = [9,18,27,50,52,58,62]
	Tomo['Y_list'] = [9,18,27,49,53,59,61]
	Tomo['Z_list'] = [0,3,6,9,18,27,36]
	Tomo['mX_list'] = [9,18,27,50,52,58,62]
	Tomo['mY_list'] = [9,18,27,49,53,59,61]
	Tomo['mZ_list'] = [0,3,6,9,18,27,36]

	state_correction_list = np.zeros(7)
	error_state_correction_list = np.zeros(7)

	for ii, jj in enumerate(Tomo[state + '_list']):
		state_correction_list[ii] = Tomo_correction_list[jj]
		if errorbar == None:
			error_state_correction_list = np.zeros(len(state_correction_list))
		else: 		
			error_state_correction_list[ii] = error_Tomo_correction_list[jj]

	return state_correction_list, error_state_correction_list


def get_C13_correction_C1C2( errorbar = None):
	Tomo_correction_list = [
		C_C1[0],C_C1[0],C_C1[0],
		C_C2[0],C_C2[0],C_C2[0],

		C_C1C2[0],C_C1C2[0],C_C1C2[0],
		C_C1C2[0],C_C1C2[0],C_C1C2[0],
		C_C1C2[0],C_C1C2[0],C_C1C2[0]]

	if errorbar == None:
		error_Tomo_correction_list = np.zeros(len(Tomo_correction_list))
	else: 
		error_Tomo_correction_list = [
			C_C1[1],C_C1[1],C_C1[1],
			C_C2[1],C_C2[1],C_C2[1],

			C_C1C2[1],C_C1C2[1],C_C1C2[1],
			C_C1C2[1],C_C1C2[1],C_C1C2[1],
			C_C1C2[1],C_C1C2[1],C_C1C2[1]]

	return Tomo_correction_list, error_Tomo_correction_list

def get_C13_correction_C5C1( errorbar = None):
	Tomo_correction_list = [
		C_C5[0],C_C5[0],C_C5[0],
		C_C1[0],C_C1[0],C_C1[0],
		
		C_C1C5[0],C_C1C5[0],C_C1C5[0],
		C_C1C5[0],C_C1C5[0],C_C1C5[0],
		C_C1C5[0],C_C1C5[0],C_C1C5[0]]

	if errorbar == None:
		error_Tomo_correction_list = np.zeros(len(Tomo_correction_list))
	else: 
		error_Tomo_correction_list = [
			C_C5[1],C_C5[1],C_C5[1],
			C_C1[1],C_C1[1],C_C1[1],

			C_C1C5[1],C_C1C5[1],C_C1C5[1],
			C_C1C5[1],C_C1C5[1],C_C1C5[1],
			C_C1C5[1],C_C1C5[1],C_C1C5[1]]

	return Tomo_correction_list, error_Tomo_correction_list