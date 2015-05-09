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
P = 0.011
# Hyperfine related constants
gam_el = 1.760859 *10**11
gam_n = 67.262 *10**6
hbar = 1.05457173*10**(-34)
mu0 = 4*pi*10**(-7)

#Filename where NV- centres are saved in
filename = r'D:\measuring\data\QEC_data\NV_sims\NV_Sim_1p1_'

##Carbon Lattice Definition
#Rotation matrix to get b along z-axis
Rz=array([[cos(pi/4),-sin(pi/4),0],[sin(pi/4),cos(pi/4),0],[0,0,1]])
Rx=array([[1,0,0],[0,cos(arctan(sqrt(2))),-sin(arctan(sqrt(2)))],[0,sin(arctan(sqrt(2))),cos(arctan(sqrt(2)))]])
# Primitive basisvectors
a = array([0,0,0])
b =a0/4*array([1,1,1])
b = Rx.dot(Rz).dot(b)
# Basisvectors of Bravais lattice
i =a0/2*array([0,1,1])
i = Rx.dot(Rz).dot(i)
j = a0/2*array([1,0,1])
j=Rx.dot(Rz).dot(j)
k =a0/2*array([1,1,0])
k = Rx.dot(Rz).dot(k)

##Lattice generation
#Location of the NV- centre is in the middle of the grid
NVPos = round(N/2) *i +round(N/2)*j+round(N/2)*k
NVx = NVPos[0]
NVy = NVPos[1]
NVz = NVPos[2]
# forloop generates the location elements
x = zeros(2*(N)**3)
y= zeros(2*(N)**3)
z=zeros(2*(N)**3)
o = 0
for n in range(N):
    for m in range(N):
        for l in range(N):
            pos = n*i + m*j+l*k
            pos2 = pos + b
            x[o] = pos[0]
            y[o] = pos[1]
            z[o] = pos[2]
            o=o+1
            x[o] = pos2[0]
            y[o] = pos2[1]
            z[o] = pos2[2]
            o = o+1
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
    Aparr = mu0*gam_el*gam_n/(4*pi)*hbar**2*power(r,-3)*(3*power(zs,2)*power(r,-2)-1) #Note 3 [0-1] - 1 can be negative and it should be able to be negative!
    Aorth = mu0*gam_el*gam_n/(4*pi)*hbar**2*power(r,-3)*3*(sqrt(power(xs,2)+power(ys,2))*zs*power(r,-2)) # also adapted this!
    #Save Hyperfine strengths to file
    q= transpose(vstack((xs,ys,zs,Aparr/hbar/(2*pi),Aorth/hbar/(2*pi))))
    savetxt(filename +str(p)+'.dat', q,delimiter='\t') #, header ='x[m]\t y[m] \t z[m] \t Aparr[Hz] \t Aorth[Hz]')

# ##Plotting
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(xs,ys,zs,c='r')
# ax.scatter(0,0,0,c='b')
# ax.grid(True)

# #Histogram of distances
# plt.figure()
# hist, bins = histogram(r, bins=50)
# width = 0.7 * (bins[1] - bins[0])
# center = (bins[:-1] + bins[1:]) / 2
# plt.bar(center, hist, align='center', width=width)

# # #Histograms of parallel Hyperfine interaction
# plt.figure()
# hist2, bins2 = histogram(Aparr/hbar/(2*pi), bins=logspace(3,14,200))
# width = 0.7 * (bins2[1] - bins2[0])
# center = (bins2[:-1] + bins2[1:]) / 2
# plt.bar(center, hist2, align='center', width=width)
# plt.gca().set_xscale("log")
# plt.gca().grid()



# plt.show()

print 'done'
