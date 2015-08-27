'''
MAB 12-5-15
Script to map nuclear spins around an inpurity
'''
from sympy.solvers import nsolve
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
from math import floor, log10
import os

'''Physical paramaters'''
mu0 =  4.*np.pi * 1e-7 #H/m
gamma_c = 0.67262e8 #rad s^-1 T^-1
gamma_h = 2.67515255e8 #rad s^-1 T^-1
gamma_e = 1.760859708e11 #rad s^-1 T^-1
#r = 1e-9 #nm
hbar = 1.05457173e-34 #Js

fontsize = 15
'''
b: blue
g: green
r: red
c: cyan
m: magenta
y: yellow
k: black
'''


carbon_params = {}
carbon_params['C1']  = {'label':'C2', 'par' : -36.0e3, 'perp':25.0e3, 'cl':'b'}
carbon_params['C2']  = {'label':'C1', 'par' : 21.2e3, 'perp':43.0e3, 'cl':'c'}
carbon_params['C3']  = {'label':'C4', 'par' : -12.0e3, 'perp':50.0e3, 'cl':'g'}
carbon_params['C5']  = {'label':'C3', 'par' : 24.7e3, 'perp':26.0e3, 'cl':'orange'}
carbon_params['C6']  = {'label':'C5', 'par' : -48.7e3, 'perp':12.0e3, 'cl':'m'}

carbon_radthet = {}
carbon_radthet['C3'] = {'label':'C4', 'r':8.40497018151e-10, 'r_err':0.0575929146301e-10, 't':0.832629794493, 't_err':0.0105199736889}
carbon_radthet['C2'] = {'label':'C1', 'r':7.94700291357e-10, 'r_err':0.0578080254284e-10, 't':1.16616965881, 't_err':0.008824378407}
carbon_radthet['C1'] = {'label':'C2', 'r':9.48815017638e-10, 'r_err':0.0736277928139e-10, 't':0.399107964421, 't_err':0.0146329798043}
carbon_radthet['C6'] = {'label':'C5', 'r':9.22504231636e-10, 'r_err':0.0615681577991e-10, 't':0.160717889439, 't_err':0.0132003653301}
carbon_radthet['C5'] = {'label':'C3', 'r':8.50891184787e-10, 'r_err':0.0815073054208e-10, 't':1.28727042467, 't_err':0.0106360511769}


'''
[0, 1, 0, 1]
[ 0.06029272  0.          0.3939056   0.0602927   0.39859622]
[-1.35526964  0.          3.21790382 -0.07917014  2.52893894]

[1, 1, 0, 1]
2.40096326  4.53280171  1.2933426   3.84273371
'''

carbon_radthet['C1']['phi'] = -1.35526964
carbon_radthet['C1']['phi_err'] = 0.06029272
carbon_radthet['C2']['phi'] = 0
carbon_radthet['C2']['phi_err'] = 0
carbon_radthet['C5']['phi'] = -0.07917014
carbon_radthet['C5']['phi_err'] = 0.0602927
carbon_radthet['C3']['phi'] = 3.21790382
carbon_radthet['C3']['phi_err'] = 0.3939056
carbon_radthet['C6']['phi'] = 2.52893894
carbon_radthet['C6']['phi_err'] = 0.39859622

carbon_radthet['C1']['phi'] = 0
carbon_radthet['C2']['phi'] = 2.40096326
carbon_radthet['C5']['phi'] = 4.53280171
carbon_radthet['C3']['phi'] = 1.2933426
carbon_radthet['C6']['phi'] = 3.84273371


carbon_radthet['C1']['posneg'] = 'pos'
carbon_radthet['C2']['posneg'] = 'neg'
carbon_radthet['C3']['posneg'] = 'neg'
carbon_radthet['C5']['posneg'] = 'pos'
carbon_radthet['C6']['posneg'] = 'neg'


'''
26 = 1.892987 +- 0.024541 #1
13 = 2.308485 +- 0.065680 #2
16 = 1.148795 +- 0.009923 #3
36 = 2.351193 +- 0.052753 #4
12 = 0.825097 +- 0.010388 #5
25 = 0.454692 +- 0.038310 #6
15 = 1.229432 +- 0.025608 #7
35 = 0.566582 +- 0.024878 #8
23 None #9
56 None #10
'''
from scipy.optimize import newton_krylov, minimize, brentq, fmin_l_bfgs_b#from mpldatacursor import datacursor
import scipy.optimize as optimize

def get_hyperfine(radius,theta):
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * radius**3)
    Apar = prefactor * (1. - 3. * np.cos(theta) **2.)
    Aper = abs(prefactor * 3. * np.cos(theta)*np.sin(theta))
    return Apar, Aper



def un2str(x, xe, precision=1):
    """pretty print nominal value and uncertainty

    x  - nominal value
    xe - uncertainty
    precision - number of significant digits in uncertainty

    returns shortest string representation of `x +- xe` either as
        x.xx(ee)e+xx
    or as
        xxx.xx(ee)"""
    # base 10 exponents
    x_exp = int(floor(log10(x)))
    xe_exp = int(floor(log10(xe)))

    # uncertainty
    un_exp = xe_exp-precision+1
    un_int = round(xe*10**(-un_exp))

    # nominal value
    no_exp = un_exp
    no_int = round(x*10**(-no_exp))

    # format - nom(unc)exp
    fieldw = x_exp - no_exp
    fmt = '%%.%df' % fieldw
    result1 = (fmt + '(%.0f)e%d') % (no_int*10**(-fieldw), un_int, x_exp)

    # format - nom(unc)
    fieldw = max(0, -no_exp)
    fmt = '%%.%df' % fieldw
    result2 = (fmt + '(%.0f)') % (no_int*10**no_exp, un_int*10**max(0, un_exp))

    # return shortest representation
    if len(result2) <= len(result1):
        return result2
    else:
        return result1

print carbon_radthet


def hyperfine_pars(x, *args):
    #For use with fsolve
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
    Apar = prefactor * (1. - 3. * np.cos(x[1]) **2.)
    Aper = abs(prefactor * 3. * np.cos(x[1])*np.sin(x[1]))
    if x[1] > np.pi/2 or x[1] < 0.:
        return 10e6, 10e6
    else:
        return (Apar - args[0]),(Aper - args[1])

for key in carbon_params.keys():
    xFinal = optimize.fsolve(hyperfine_pars,(1e-9,0.45*np.pi),args=(carbon_params[key]['par'],carbon_params[key]['perp']),full_output=True)
    carbon_params[key]['r'] = xFinal[0][0]
    carbon_params[key]['theta'] = xFinal[0][1]

for keys in carbon_radthet.keys():
    print keys, get_hyperfine(carbon_radthet[keys]['r'],carbon_radthet[keys]['t'])




def spher2cart(r,theta,phi):
	x = r*np.sin(theta) * np.cos(phi)
	y = r*np.sin(theta) * np.sin(phi)
	z = r*np.cos(theta)
	return x, y, z

def cart2spher(x,y,z):
	r = (x**2 + y**2 + z**2) ** 0.5
	theta = np.arccos(z / r)
	phi = np.arctan(y / x)
	return r, theta, phi

def distance_in_spher(r1,r2):
    return (r1[0]**2+r2[0]**2 -2.*r1[0]*r2[0]*(np.sin(r1[1])*np.sin(r2[1])*np.cos(r1[2]-r2[2])+np.cos(r1[1])*np.cos(r2[1])))**0.5



def r2coupling(gamma_i,gamma_j, radius):
	mu0 =  4.*np.pi * 1e-7
	return abs(((mu0 * gamma_i*gamma_j * hbar) / (8.*np.pi**2. * r**3) ))

def coupling2r(gamma_i,gamma_j, coupling, theta):
	mu0 =  4.*np.pi * 1e-7
	return ((mu0 * gamma_i * gamma_j * hbar) * (3. * np.cos(theta)**2. - 1.) / (8.*np.pi**2. * coupling)) **(1./3)

def get_spher_from_hf(carbon_nr):
	par = abs(carbon_params['C'+str(carbon_nr)]['par'])
	per = abs(carbon_params['C'+str(carbon_nr)]['perp'])
	theta = np.arctan(per / par)
	print 'theta', theta * 180 / np.pi
	radius = coupling2r(gamma_e,gamma_c,par,theta)
	print 'radius', radius*1e9
	return radius, theta

def plot_hf_circle(ax,carbon_nr,posneg='both'):
    C_str= 'C' + str(carbon_nr)
    label = carbon_params[C_str]['label']
    theta = carbon_radthet[C_str]['t']
    theta_err = carbon_radthet[C_str]['t_err']
    r = carbon_radthet[C_str]['r']*1e9
    r_err = carbon_radthet[C_str]['r_err']*1e9
    color = carbon_params[C_str]['cl']
    phi = np.linspace(-1.*np.pi,1.*np.pi,101)
    if posneg in ['both', 'pos']:
        for r2 in np.linspace(r-r_err*2,r+r_err*2,6):
            for t2 in np.linspace(theta-theta_err*2,theta+theta_err*2,18):
                x,y,z = spher2cart(r2, t2, phi)
                # opacity = np.exp(-(r2-r)**2./(2*(r_err/1.)**2.))*np.exp(-(t2-theta)**2./(2*(theta_err/1.)**2.))
                ax.plot(x, y, z, color= color,alpha=0.025)
        # x,y,z = spher2cart(r, theta, phi)
        templine1 = ax.plot(x, y, 100, color=color,label=label,alpha=1.)
    if posneg in ['both', 'neg']:
        theta = np.pi - theta
        for r2 in np.linspace(r-r_err*2,r+r_err*2,6):
            for t2 in np.linspace(theta-theta_err*2,theta+theta_err*2,18):
                x,y,z = spher2cart(r2, t2, phi)
                # opacity = np.exp(-(r2-r)**2./(2*(r_err/1.)**2.))*np.exp(-(t2-theta)**2./(2*(theta_err/1.)**2.))
                ax.plot(x, y, z, color= color,alpha=0.025)
        # x,y,z = spher2cart(r, theta, phi)
        templine2 = ax.plot(x, y, 100, color=color,alpha=1.)
    return [templine1,templine2]

def plot_hf_circle2(ax,carbon_nr,posneg='both'):
    C_str= 'C' + str(carbon_nr)
    theta = carbon_radthet[C_str]['t']
    r = carbon_radthet[C_str]['r']*1e9
    color = carbon_params[C_str]['cl']
    phi = np.linspace(-1.*np.pi,1.*np.pi,101)
    if posneg in ['both', 'pos']:
        x,y,z = spher2cart(r, theta, phi)
        ax.plot(x, y, z, color= color, label=C_str, alpha=0.)
    if posneg in ['both', 'neg']:
        theta = np.pi - theta
        x,y,z = spher2cart(r, theta, phi)
        ax.plot(x, y, z, color= color, label=C_str, alpha=0.)

def plot_carbon(ax,carbon_nr,phi,posneg='both'):
    C_str = 'C' + str(carbon_nr)
    print C_str

    theta = carbon_radthet[C_str]['t']
    r = carbon_radthet[C_str]['r']*1e9
    color = carbon_params[C_str]['cl']
    if posneg in ['pos']:
		x,y,z = spher2cart(r, theta, phi)
		ax.scatter([x], [y], [z], c=color, marker='o', s=np.pi * (6 * np.ones((3,)))**2)
        # ax.scatter([x], [y], [z], c='w', marker='o', label=C_str, s=np.pi * (4 * np.ones((3,)))**2, linewidths=0)
    if posneg in ['neg']:
		theta = np.pi - theta
		x,y,z = spher2cart(r, theta, phi)
		ax.scatter([x], [y], [z], c=color, marker='o', s=np.pi * (6 * np.ones((3,)))**2)
        # ax.scatter([x], [y], [z], c='w', marker='o', label=C_str, s=np.pi * (4 * np.ones((3,)))**2, linewidths=0)

plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X (nm)', fontsize = fontsize)
ax.set_ylabel('Y (nm)', fontsize = fontsize)
ax.set_zlabel('Z (nm)', fontsize = fontsize)
ax.scatter([0.], [0.], [0.], c='r', marker=u'o', label='NV',s=np.pi * (6 * np.ones((3,)))**2,linewidths=0.01)
# for carbon in [1,3,5]:
# print 'Hello'
lines = []
for carbon in [2,1,5,3,6]:
    C_str = 'C' + str(carbon)
    templine = plot_hf_circle(ax,carbon)
    lines.extend(templine)
    plot_carbon(ax,carbon,carbon_radthet[C_str]['phi'],posneg=carbon_radthet[C_str]['posneg'])
# for carbon in [3,5]:
#     plot_hf_circle2(ax,carbon)
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.view_init(elev=20)
ax.set_xlim3d(-1.0, 1.0)
ax.set_ylim3d(-1.0,1.0)
ax.set_zlim3d(-1.0,1.0)

# plot_carbon(ax,1,0,posneg='pos')
# plot_carbon(ax,2,1.35526965, posneg='pos')
# plot_carbon(ax,3,4.57317311, posneg='neg')
# plot_carbon(ax,5,1.2760995, posneg='pos')
# plot_carbon(ax,6,3.88420822, posneg='neg')
ax.legend(loc=(-0.1,0.3),fontsize=15,numpoints = 1,scatterpoints=1,frameon = False,columnspacing=0.5,handletextpad=1)
plt.savefig(r'/Users/'+os.getlogin()+r'/Dropbox/Thesis/Images/hyperfinecircles.png')
plt.show()





### fitted to the ms = -1 data




# def randrange(n, vmin, vmax):
#     return (vmax-vmin)*np.random.rand(n) + vmin


# n = 100



# # for c, m, zl, zh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
# #     xs = randrange(n, 23, 32)
# #     ys = randrange(n, 0, 100)
# #     zs = randrange(n, zl, zh)
# #     ax.scatter(xs, ys, zs, c=c, marker=m)



# plt.show()


# fig = plt.figure()
# ax = fig.gca(projection='3d')
# theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
# z = np.linspace(-2, 2, 100)
# r = z**2 + 1
# x = r * np.sin(theta)
# y = r * np.cos(theta)
# ax.plot(x, y, z, label='parametric curve')
