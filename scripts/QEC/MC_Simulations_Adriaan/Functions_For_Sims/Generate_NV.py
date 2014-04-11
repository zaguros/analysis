
def Generate_NV(Carb_Conc=0.0011,N_NV=2,N=25):
    """
    Takes Carbon Concentration, number of centres and gridsize as inputs
    Returns a list of 2D numpy arrays with hyperfine strengths in Hz (real freq)
    """
    import numpy as np
    #Constants
    #lattice parameter
    a0 = 3.57 * 10**(-10)
    # Hyperfine related constants
    gam_el = 1.760859 *10**11 #Gyromagnetic ratio rad s-1 T-1
    gam_n = 67.262 *10**6 #rad s-1 T-1
    hbar = 1.05457173*10**(-34)
    pi = np.pi
    mu0 = 4*pi*10**(-7)

    ##Carbon Lattice Definition
    #Rotation matrix to get b along z-axis
    Rz=np.array([[np.cos(pi/4),-np.sin(pi/4),0],[np.sin(pi/4),np.cos(pi/4),0],[0,0,1]])
    Rx=np.array([[1,0,0],[0,np.cos(np.arctan(np.sqrt(2))),-np.sin(np.arctan(np.sqrt(2)))],[0,np.sin(np.arctan(np.sqrt(2))),np.cos(np.arctan(np.sqrt(2)))]])
    # Basis vectors
    a = np.array([0,0,0])
    b =a0/4*np.array([1,1,1])
    b = Rx.dot(Rz).dot(b)
    # Basisvectors of Bravais lattice
    i =a0/2*np.array([0,1,1])
    i = Rx.dot(Rz).dot(i)
    j = a0/2*np.array([1,0,1])
    j=Rx.dot(Rz).dot(j)
    k =a0/2*np.array([1,1,0])
    k = Rx.dot(Rz).dot(k)
    NVPos = round(N/2) *i +round(N/2)*j+round(N/2)*k

    prefactor = mu0*gam_el*gam_n/(4*pi)*hbar**2 /hbar/(2*pi) #Last /hbar/2pi is to convert from Joule to Hz

    #Initialise
    L_size = 2*(N)**3
    A = np.empty([L_size,2])
    o=0

    #Calculate Hyperfine strength for all gridpoints
    for n in range(N):
        for m in range(N):
            for l in range(N):
                if (n== round(N/2) and m==round(N/2) and l == round(N/2)) :#Omit the Nitrogen and the Vacancy centre in the calculations
                    A[o,0]=0
                    A[o,1]= 0
                    o+=1
                    A[o,0]=0
                    A[o,1]= 0
                    o+=1
                else:
                    pos1 = n*i + m*j+l*k -NVPos
                    pos2 = pos1 + b
                    r = np.sqrt(pos1.dot(pos1))
                    A[o,0] =prefactor*np.power(r,-3)*(3*np.power(pos1[2],2)*np.power(r,-2)-1)
                    A[o,1] = prefactor*np.power(r,-3)*3*(np.sqrt(np.power(pos1[0],2)+np.power(pos1[1],2))*pos1[2]*np.power(r,-2))
                    o +=1
                    r = np.sqrt(pos2.dot(pos2))
                    A[o,0] =prefactor*np.power(r,-3)*(3*np.power(pos2[2],2)*np.power(r,-2)-1)
                    A[o,1] = prefactor*np.power(r,-3)*3*(np.sqrt(np.power(pos2[0],2)+np.power(pos2[1],2))*pos2[2]*np.power(r,-2))
                    o+=1
    # Generate different NV-Objects by randomly selecting which gridpoints contain a carbon.
    NV_list = []
    for p in range(N_NV):
        Sel = np.where(np.random.rand(L_size)<Carb_Conc)
        A_NV =[ A[u,:] for u in Sel]
        NV_list.append(A_NV[0]) #index 0 is to get rid of outher brackets in A_NV
    return NV_list


