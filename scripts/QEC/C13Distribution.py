from numpy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

##Constants
#Grid size
N = 25
NumSims = 1
#lattice parameter
a0 = 3.57 * 10**(-10)
#Probability that a carbon atom is C13, for natural carbon 1.1%
P = 1.
# Hyperfine related constants
gam_el = 1.760859 *10**11
gam_n = 67.262 *10**6
hbar = 1.05457173*10**(-34)
mu0 = 4*pi*10**(-7)

#Filename where NV- centres are saved in
filename = r'D:\measuring\data\QEC_data\NV_sims\NV_Sim_1p1_'

##Carbon Lattice Definition
#Rotation matrix to get b along z-axis
Rz=array([[cos(np.pi/4),-sin(np.pi/4),0],[sin(np.pi/4),cos(np.pi/4),0],[0,0,1]])
Rx=array([[1,0,0],[0,cos(arctan(sqrt(2))),-sin(arctan(sqrt(2)))],[0,sin(arctan(sqrt(2))),cos(arctan(sqrt(2)))]])
# Primitive basisvectors
a = array([0,0,0])
b =a0/4*array([1,1,1])
b = Rx.dot(Rz).dot(b)
for m in [0,1,2]:
    b[m] = round(b[m],24) # for rounding errors in getting b along z-axis 
print 'b '+str(b)
# Basisvectors of Bravais lattice
i =a0/2*array([0,1,1])
i = Rx.dot(Rz).dot(i)
j = a0/2*array([1,0,1])
j=Rx.dot(Rz).dot(j)
k =a0/2*array([1,1,0])
k = Rx.dot(Rz).dot(k)
for m in [0,1,2]:
    b[m] = round(b[m],24) # for rounding errors in getting b along z-axis 
    i[m] = round(i[m],24) # for rounding errors in getting b along z-axis 
    j[m] = round(j[m],24) # for rounding errors in getting b along z-axis 
    k[m] = round(b[m],24) # for rounding errors in getting b along z-axis     
##Lattice generation
#Location of the NV- centre is in the middle of the grid
NVPos = round(N/2) *i +round(N/2)*j+round(N/2)*k
NVx = NVPos[0]
NVy = NVPos[1]
NVz = NVPos[2]
# print 'NV pos'
# print NVPos
# forloop generates the location elements
x = zeros(2*(N)**3-1)
y= zeros(2*(N)**3-1)
z=zeros(2*(N)**3-1)
o = 0
for n in range(N):
    for m in range(N):
        for l in range(N):
            # print 'n ' +str(n)+ 'm ' +str(m)+'l ' +str(l)
            pos = n*i + m*j+l*k
            x[o] = pos[0]
            y[o] = pos[1]
            z[o] = pos[2]
            # print 'pos ' +str(o)
            if x[o]== NVx and y[o]==NVy and z[o]==NVz:
                o = o
                # print 'NV pos!'
                # print o
            else:
                o=o+1
            x[o] = pos2[0]
            y[o] = pos2[1]
            z[o] = pos2[2]
            # print 'pos2 ' +str(o)
            # if x[o]== NVx and y[o]==NVy and z[o]==NVz:
            #     o = o
                # print 'NV pos!'
            # else:
            o=o+1
# print o 
# print 'len x ' +str(len(x))
# print x
# #Remove the NV- Center from the carbon lattice
# indNV = [2*(N/2)**3,2*(N/2)**3-1]

# x = delete(x,indNV)
# y= delete(y,indNV)
# z= delete(z,indNV)
#Offsetting all positions relative to the NV
x = x -NVx
y = y -NVy
z = z -NVz


##Determining which lattice points are C13's
for p in xrange(NumSims):
    Sel = random.rand(size(x))<P
    iSel = where(Sel==1)
    xs = [x[u] for u in iSel]
    ys = [y[u] for u in iSel]
    zs = [z[u] for u in iSel]

    print 'Number of C 13s found in grid: ' + str(size(xs))
    print str(p)
     ## Calculating the Hyperfine strength
    r = power(power(xs,2)+power(ys,2)+power(zs,2),.5)
    print ' max'
    print max(r[0])
    Aparr = mu0*gam_el*gam_n/(4*pi)*hbar**2*power(r,-3)*(3*power(zs,2)*power(r,-2)-1) #Note 3 [0-1] - 1 can be negative and it should be able to be negative!
    Aorth = mu0*gam_el*gam_n/(4*pi)*hbar**2*power(r,-3)*3*(sqrt(power(xs,2)+power(ys,2))*zs*power(r,-2)) # also adapted this!
    #Save Hyperfine strengths to file
    print min(np.sqrt((Aparr/hbar/(2*pi))**2+(Aorth/hbar/(2*pi))**2)[0])
    q= transpose(vstack((xs,ys,zs,Aparr/hbar/(2*pi),Aorth/hbar/(2*pi))))
    savetxt(filename +str(p)+'.dat', q,delimiter='\t') #, header ='x[m]\t y[m] \t z[m] \t Aparr[Hz] \t Aorth[Hz]')


print 'done'
