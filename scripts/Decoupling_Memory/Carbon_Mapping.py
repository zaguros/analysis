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


'''Physical paramaters'''
mu0 =  4.*np.pi * 1e-7 #H/m
gamma_c = 0.67262e8 #rad s^-1 T^-1
gamma_h = 2.67515255e8 #rad s^-1 T^-1
gamma_e = 1.760859708e11 #rad s^-1 T^-1
#r = 1e-9 #nm
hbar = 1.05457173e-34 #Js

fontsize = 21
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
carbon_params['C1']  = {'par' : -36.0e3, 'perp':25.0e3, 'cl':'b'}
carbon_params['C2']  = {'par' : 21.2e3, 'perp':43.0e3, 'cl':'c'}
carbon_params['C3']  = {'par' : -12.0e3, 'perp':50.0e3, 'cl':'g'}
carbon_params['C5']  = {'par' : 24.7e3, 'perp':26.0e3, 'cl':'y'}
carbon_params['C6']  = {'par' : -48.7e3, 'perp':12.0e3, 'cl':'m'}


C2C_coupling = {}
C2C_coupling['12'] = {'C' : 0.825097, 'C_u' : 0.010388}
C2C_coupling['13'] = {'C' : 2.308485, 'C_u' : 0.065680}
C2C_coupling['15'] = {'C' : 1.229432, 'C_u' : 0.025608}
C2C_coupling['16'] = {'C' : 1.148795, 'C_u' : 0.009923}
C2C_coupling['23'] = {'C' : None, 'C_u' : None}
C2C_coupling['25'] = {'C' : 0.454692, 'C_u' : 0.038310}
C2C_coupling['26'] = {'C' : 1.892987, 'C_u' : 0.024541}
C2C_coupling['35'] = {'C' : 0.566582, 'C_u' : 0.024878}
C2C_coupling['36'] = {'C' : 2.351193, 'C_u' : 0.052753}
C2C_coupling['56'] = {'C' : None, 'C_u' : None}

for key in C2C_coupling.keys():
    C2C_coupling[key[::-1]]=C2C_coupling[key]

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


def thetasolve(theta, Apar, Aper):
    return  (3. * np.sin(theta) * np.cos(theta)) / (3. * np.cos(theta) **2. - 1.) - 1.* (Apar/Aper)

def radiussolve(radius,Apar,theta):
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    return (-1.*(mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * radius**3) ) * (3. * np.cos(theta) **2. - 1.) + Apar

def get_hyperfine(radius,theta):
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * radius**3)
    Apar = prefactor * (1. - 3. * np.cos(theta) **2.)
    Aper = abs(prefactor * 3. * np.cos(theta)*np.sin(theta))
    return Apar, Aper

def hyperfine_fit(x,*args):
    print args[0]
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
    Apar = prefactor * (1. - 3. * np.cos(x[1]) **2.)
    Aper = prefactor * 3. * np.cos(x[1])*np.sin(x[1])
    return (Apar - args[0] + Aper - args[1])**2.



# def hyperfine_gradient(x):
#     mu0 =  4.*np.pi * 1e-7
#     gamma_c = 0.67262e8
#     gamma_e = 1.760859708e11
#     hbar = 1.05457173e-34
#     prefactor3 = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
#     prefactor4 =
#     dfdtheta = prefactor3 () 


    # return 

# def hyperfine_par(x):
#     mu0 =  4.*np.pi * 1e-7
#     gamma_c = 0.67262e8
#     gamma_e = 1.760859708e11
#     hbar = 1.05457173e-34
#     prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
#     Apar = prefactor * (1. - 3. * np.cos(x[1]) **2.)
#     return Apar + 36.0e3

def hyperfine_pars(x, *args):
    #For use with fsolve
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
    Apar = prefactor * (1. - 3. * np.cos(x[1]) **2.)
    Aper = abs(prefactor * 3. * np.cos(x[1])*np.sin(x[1]))
    return (Apar - args[0]),(Aper - args[1])

def calc_C2C_coupling(r1,r2):
    mu0 =  4.*np.pi * 1e-7 #H/m
    gamma_c = 0.67262e8 #rad s^-1 T^-1
    hbar = 1.05457173e-34 #Js

    x1,y1,z1 = spher2cart(r1[0],r1[1],r1[2])
    x2,y2,z2 = spher2cart(r2[0],r2[1],r2[2])
    distance = ((x1-x2)**2. + (y1-y2)**2. + (z1-z2)**2.)**0.5
    cos2 = ((z1-z2) / distance)**2.
    # print 'distance', distance,'cos2', cos2
    return abs( (mu0 * gamma_c**2. * hbar) / (8.*np.pi**2. * distance**3.) * (1. - 3. * cos2) )




for key in carbon_params.keys():
    xFinal= optimize.fsolve(hyperfine_pars,(1e-9,0.45*np.pi),args=(carbon_params[key]['par'],carbon_params[key]['perp']),full_output=True)
    # x = optimize.fmin_l_bfgs_b(hyperfine_pars_single_f, np.array([xFinal[0][0],xFinal[0][1]]),
    #     fprime=hyperfine_pars_single_fprime,args=(carbon_params[key]['par'],carbon_params[key]['perp']))
    print xFinal[0]
    print sum(xFinal[1]['fvec'])
    print (carbon_params[key]['par'],carbon_params[key]['perp'])
    xFinal[0][0]
    carbon_params[key]['r'] = xFinal[0][0]
    carbon_params[key]['theta']= xFinal[0][1]
    print get_hyperfine(xFinal[0][0],xFinal[0][1])

def calcal_fsolve_func(x, *args):
    coor1 = (args[1],args[6],args[0])
    coor2 = (args[2],args[7],x[0])
    coor3 = (args[3],args[8],x[1])
    coor5 = (args[4],args[9],x[2])
    coor6 = (args[5],args[10],x[3])

    C12 = calc_C2C_coupling(coor1,coor2)-0.825097
    C13 = calc_C2C_coupling(coor1,coor3)-2.308485
    C15 = calc_C2C_coupling(coor1,coor5)-1.229432
    C16 = calc_C2C_coupling(coor1,coor6)-1.148795
    C25 = calc_C2C_coupling(coor2,coor5)-0.454692
    C26 = calc_C2C_coupling(coor2,coor6)-1.892987
    C35 = calc_C2C_coupling(coor3,coor5)-0.566582
    C36 = calc_C2C_coupling(coor3,coor6)-2.351193
    return C12**2,C13**2,C35**2,C25**2

def Apar_pars(x):
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
    Apar = prefactor * (1. - 3. * np.cos(x[1]) **2.)
    return Apar/10e3

def Aper_pars(x):
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
    Aper = abs(prefactor * 3. * np.cos(x[1])*np.sin(x[1]))
    return Aper/10e3

def spher2cart(r,theta,phi):
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x,y,z

def carcar_resid_func2(p, y):
    r = p[4:9]
    theta = p[9:15]
    phi = p[0:4]

    coor1 = (p[4],p[9],0.)
    coor2 = (p[5],p[10],p[0])
    coor3 = (p[6],p[11],p[1])
    coor5 = (p[7],p[12],p[2])
    coor6 = (p[8],p[13],p[3])

    C12 = y[0]-calc_C2C_coupling(coor1,coor2)
    C13 = y[1]-calc_C2C_coupling(coor1,coor3)
    C15 = y[2]-calc_C2C_coupling(coor1,coor5)
    C16 = y[3]-calc_C2C_coupling(coor1,coor6)
    C25 = y[4]-calc_C2C_coupling(coor2,coor5)
    C26 = y[5]-calc_C2C_coupling(coor2,coor6)
    C35 = y[6]-calc_C2C_coupling(coor3,coor5)
    C36 = y[7]-calc_C2C_coupling(coor3,coor6)

    Par1 = y[8]-Apar_pars((coor1[0],coor1[1]))
    Par2 = y[9]-Apar_pars((coor2[0],coor2[1]))
    Par3 = y[10]-Apar_pars((coor3[0],coor3[1]))
    Par5 = y[11]-Apar_pars((coor5[0],coor5[1]))
    Par6 = y[12]-Apar_pars((coor6[0],coor6[1]))

    Per1 = y[13]-Aper_pars((coor1[0],coor1[1]))
    Per2 = y[14]-Aper_pars((coor2[0],coor2[1]))
    Per3 = y[15]-Aper_pars((coor3[0],coor3[1]))
    Per5 = y[16]-Aper_pars((coor5[0],coor5[1]))
    Per6 = y[17]-Aper_pars((coor6[0],coor6[1]))

    if np.any(phi> 2.*np.pi) or np.any(phi < 0) or phi[0]>np.pi or np.any(r<0.) or np.any(r>5e-9) or np.any(theta < 0) or np.any(theta > np.pi) or theta[0] > np.pi/2.:
        return np.ones((18,))*1e6
    else:
        return np.array([C12,C13,C15,C16,C25,C26,C35,C36,Par1,Par2,Par3,Par5,Par6,Per1,Per2,Per3,Per5,Per6])

def carcar_resid_func(p, y, args2):
    coor1 = (args2[1],args2[6],args2[0])
    coor2 = (args2[2],args2[7],p[0])
    coor3 = (args2[3],args2[8],p[1])
    coor5 = (args2[4],args2[9],p[2])
    coor6 = (args2[5],args2[10],p[3])

    C12 = y[0]-calc_C2C_coupling(coor1,coor2)
    C13 = y[1]-calc_C2C_coupling(coor1,coor3)
    C15 = y[2]-calc_C2C_coupling(coor1,coor5)
    C16 = y[3]-calc_C2C_coupling(coor1,coor6)
    C25 = y[4]-calc_C2C_coupling(coor2,coor5)
    C26 = y[5]-calc_C2C_coupling(coor2,coor6)
    C35 = y[6]-calc_C2C_coupling(coor3,coor5)
    C36 = y[7]-calc_C2C_coupling(coor3,coor6)

    if np.any(p > 2.*np.pi) or np.any(p < 0) or p[0]>np.pi:
            return np.ones((8,))*1e6
    else:
        return np.array([C12,C13,C15,C16,C25,C26,C35,C36])



def head_func(args2):
    def carcar_curvefit_func(x, phi2, phi3, phi5,phi6):
        coor1 = (x[1],x[6],x[0])
        coor2 = (x[2],x[7],phi2)
        coor3 = (x[3],args2[8],phi3)
        coor5 = (x[4],args2[9],phi5)
        coor6 = (x[5],args2[10],phi6)

        C12 = calc_C2C_coupling(coor1,coor2)
        C13 = calc_C2C_coupling(coor1,coor3)
        C15 = calc_C2C_coupling(coor1,coor5)
        C16 = calc_C2C_coupling(coor1,coor6)
        C25 = calc_C2C_coupling(coor2,coor5)
        C26 = calc_C2C_coupling(coor2,coor6)
        C35 = calc_C2C_coupling(coor3,coor5)
        C36 = calc_C2C_coupling(coor3,coor6)
        p = np.array([phi2,phi3,phi5,phi6])
        if np.any(p > 2.*np.pi) or np.any(p < 0) or p[0]>np.pi:
            return np.ones((8,))*1e6
        else:
            return np.array([C12,C13,C15,C16,C25,C26,C35,C36])
    return carcar_curvefit_func

random.seed()
print 'Start'
p1 = 0.
r1 = carbon_params['C1']['r']
r2 = carbon_params['C2']['r']
r3 = carbon_params['C3']['r']
r5 = carbon_params['C5']['r']
r6 = carbon_params['C6']['r']
t1 = carbon_params['C1']['theta']

diglist= [
[0, 0, 0, 0],
[0, 0, 0, 1],
[0, 0, 1, 0],
[0, 0, 1, 1],
[0, 1, 0, 0],
[0, 1, 0, 1],
[0, 1, 1, 0],
[0, 1, 1, 1],
[1, 0, 0, 0],
[1, 0, 0, 1],
[1, 0, 1, 0],
[1, 0, 1, 1],
[1, 1, 0, 0],
[1, 1, 0, 1],
[1, 1, 1, 0],
[1, 1, 1, 1]]


coupling_array =np.array([0.825097,2.308485,1.229432,1.148795,0.454692,1.892987,0.566582,2.351193,\
    -3.60,2.12,-1.20,2.47,-4.87,2.50,4.30,5.00,2.60,1.20])
if False:
    for diglist1 in diglist:
        if diglist1[0] == 0:
            t2 = carbon_params['C2']['theta']
        else:
            t2 = np.pi-carbon_params['C2']['theta']
            
        if diglist1[1] == 0:
            t3 = carbon_params['C3']['theta']
        else:
            t3 = np.pi-carbon_params['C3']['theta']
                
        if diglist1[2] == 0:
            t5 = carbon_params['C5']['theta']
        else:
            t5 = np.pi-carbon_params['C5']['theta']
                    
        if diglist1[3] == 0:
            t6 = carbon_params['C6']['theta']
        else:
            t6 = np.pi-carbon_params['C6']['theta']
        print diglist1
        args2 = np.array([p1,r1,r2,r3,r5,r6,t1,t2,t3,t5,t6])
        # xdata = args2[:8]
        Antwoord = 10e4
        timeschange = 0
        result = None
        for aa in range(10000):
            try:
                x0 = (random.uniform(0.,np.pi),random.uniform(0.,2.*np.pi),random.uniform(0.,2.*np.pi),random.uniform(0.,2.*np.pi))
                # print x0
                # res= optimize.curve_fit(head_func(args2),xdata,coupling_array,p0=x0,full_output=True,maxfev=1000)
                res= optimize.leastsq(carcar_resid_func,x0,args=(coupling_array,args2),full_output=True)

                # print res
                # print res[0]
                if np.sum(np.abs(res[2]['fvec'])) < Antwoord:
                    print aa
                    print res[0]
                    timeschange += 1
                    Antwoord = np.sum(np.abs(res[2]['fvec']))
                    result = res
            except:
                pass
            
        print '============'
        print Antwoord
        print timeschange
        print result
        print '============'
                        # print sum(xFinal[1]['fvec'])

if True: 
    # xdata = args2[:8]
    Antwoord = 10e4
    timeschange = 0
    result = None
    for aa in range(200000):
        try:
            x0 = (random.uniform(0.,np.pi),random.uniform(0.,2.*np.pi),random.uniform(0.,2.*np.pi),random.uniform(0.,2.*np.pi), \
                1e-9,1e-9,1e-9,1e-9,1e-9,random.uniform(0.,np.pi/2.),random.uniform(0.,np.pi),random.uniform(0.,np.pi),random.uniform(0.,np.pi),random.uniform(0.,np.pi))
            # print x0
            # res= optimize.curve_fit(head_func(args2),xdata,coupling_array,p0=x0,full_output=True,maxfev=1000)
            res= optimize.leastsq(carcar_resid_func2,x0,args=(coupling_array),full_output=True,maxfev=10000)

            # print res
            # print res[0]
            if np.sum(np.abs(res[2]['fvec'])) < Antwoord:
                print '-------- next -------'
                print 'complete res', res
                print 'index + timeschange', aa, timeschange
                print 'res0',res[0]
                print 'error', np.sum(np.abs(res[2]['fvec']))
                print 'error_all', res[2]['fvec']
                print 
                timeschange += 1
                Antwoord = np.sum(np.abs(res[2]['fvec']))

        except:
            pass
    print '============'
    print Antwoord
    print timeschange
    print result
    print res
    print '============'
                        # print sum(xFinal[1]['fvec'])

# dfafesfa
def hyperfine_pars_single_f(x, *args):
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
    Apar = prefactor * (1. - 3. * np.cos(x[1]) **2.)
    Aper = prefactor * 3. * np.cos(x[1])*np.sin(x[1])
    return Apar - args[0] + Aper - args[1]

def hyperfine_pars_single_fprime(x,*args):
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2.)
    fprime_r = prefactor * -3./(x[0]**4.) * ((1.-3.*np.cos(x[1])**2.) + (3. * np.cos(x[1])*np.sin(x[1])))
    fprime_t = prefactor * 1./(x[0]**3.) * (3. * np.cos(2*x[1]) + 6.* np.cos(x[1])*np.sin(x[1]) )
    return np.array([fprime_r,fprime_t])

def hyperfine_pars_curve_fit(x):
    mu0 =  4.*np.pi * 1e-7
    gamma_c = 0.67262e8
    gamma_e = 1.760859708e11
    hbar = 1.05457173e-34
    prefactor = (mu0 * gamma_e*gamma_c * hbar) / (8.*np.pi**2. * x[0]**3)
    Apar = prefactor * (1. - 3. * np.cos(x[1]) **2.)
    Aper = prefactor * 3. * np.cos(x[1])*np.sin(x[1])
    return np.array([Apar,Aper])

# key = 'C1'
# x,f,d = fmin_l_bfgs_b(hyperfine_fit,x0=np.array([1e-10,0.45*np.pi]),
#     args=(carbon_params[key]['par'],carbon_params[key]['perp']),bounds=[(0.,10e-9),(0,np.pi/2.)],approx_grad=True)
# print x
# print f
# print d

# adsfa



# AparP, AperP = get_hyperfine(9.99213002471048e-10,0.5275582194172149)
# print 'Par',AparP
# print 'Per',AperP


def get_radius_theta(key):
    t_x=np.linspace(0,np.pi/2,1001)
    A_par = carbon_params[key]['par']
    A_per = carbon_params[key]['perp']
    plt.plot(t_x,thetasolve(t_x,A_par,A_per))
    xmin= t_x[np.min(np.where(thetasolve(t_x,A_par,A_per)<0))]
    # print thetasolve(xmin,A_par,A_per)
    xmax= t_x[np.max(np.where(thetasolve(t_x,A_par,A_per)>0))]
    # print thetasolve(xmin,A_par,A_per)
    theta = brentq(thetasolve,xmin,xmax,args=(A_par,A_per))
    print 
    print key
    print 'A_par\t', A_par/1e3
    print 'A_per\t', A_per/1e3
    try:
        radius= brentq(radiussolve,0.01e-9,10e-9,args=(A_par,theta))
        print 'radius\t', radius*1e9
    except:
        pass

    print 'cos term\t', (3. * np.cos(theta) **2. - 1.)
    print 'theta \t',theta * 180/np.pi
    return radius, theta


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
	theta = carbon_params[C_str]['theta']
	r = carbon_params[C_str]['r']*1e9
	color = carbon_params[C_str]['cl']
	phi = np.linspace(-1.*np.pi,1.*np.pi,101)
	if posneg in ['both', 'pos']:
		x,y,z = spher2cart(r, theta, phi)
		ax.plot(x, y, z, color= color, label=C_str)
	if posneg in ['both', 'neg']:
		theta = np.pi - theta
		x,y,z = spher2cart(r, theta, phi)
		ax.plot(x, y, z, color= color, label=C_str)

def plot_hf_circle2(ax,carbon_nr,posneg='both'):
    C_str= 'C' + str(carbon_nr)
    theta = carbon_params[C_str]['theta']
    r = carbon_params[C_str]['r']*1e9
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
	C_str= 'C' + str(carbon_nr)
	theta = carbon_params[C_str]['theta']
	r = carbon_params[C_str]['r']*1e9
	color = carbon_params[C_str]['cl']
	if posneg in ['pos']:
		x,y,z = spher2cart(r, theta, phi)
		ax.scatter([x], [y], [z], c=color, marker='o', label=C_str, s=np.pi * (6 * np.ones((3,)))**2)
        # ax.scatter([x], [y], [z], c='w', marker='o', label=C_str, s=np.pi * (4 * np.ones((3,)))**2, linewidths=0)
	if posneg in ['neg']:
		theta = np.pi - theta
		x,y,z = spher2cart(r, theta, phi)
		ax.scatter([x], [y], [z], c=color, marker='o', label=C_str, s=np.pi * (6 * np.ones((3,)))**2)
        # ax.scatter([x], [y], [z], c='w', marker='o', label=C_str, s=np.pi * (4 * np.ones((3,)))**2, linewidths=0)

	# C_str= 'C' + str(carbon_nr)
	# theta = carbon_params[C_str]['theta']
	# r = carbon_params[C_str]['r']*1e9
	# color = carbon_params[C_str]['cl']
	# phi = np.linspace(-1.*np.pi,1.*np.pi,101)
	# if posneg in ['both', 'pos']:
	# 	x,y,z = spher2cart(r, theta, phi)
	# 	ax.plot(x, y, z, color= color, label=C_str)
	# if posneg in ['both', 'pos']:
	# 	theta = np.pi - theta
	# 	x,y,z = spher2cart(r, theta, phi)
	# 	ax.plot(x, y, z, color= color, label=C_str)

def func_for_fit(x):
    C1 = (9.99213002471048e-10,0.5275582194172149,0)
    C2 = (x[0],x[1],x[2])
    C3 = (x[3],x[4],x[5])

    return C1


# for key in carbon_params.keys():
#     # carbon_nr = int(key[1])
#     #get_radius_theta(key)
#     carbon_params[key]['r'], carbon_params[key]['theta'] = get_radius_theta(key)
#     AparP, AperP = get_hyperfine(carbon_params[key]['r'], carbon_params[key]['theta'] )
#     print key, 'Apar', AparP, 'Aper', AperP
# print carbon_params[key]

# aadfae

# key1 = 'C1'
# carbon_params[key1]['phi'] = 0
# carbon_params['C5']['phi'] = 70.5
# r1 = (carbon_params[key1]['r'],carbon_params[key1]['theta'],0)
# thetas = np.linspace(0,2*np.pi,201)
# fig = plt.figure()
# ax2 = fig.add_subplot(111)

# for key2 in carbon_params.keys():
#     if key2 != key1 and C2C_coupling[key1[1] + key2[1]]['C'] != None:
#         fig = plt.figure()
#         ax2 = fig.add_subplot(111)
#         coupling = []
#         coupling2 = []
#         for ii, theta in enumerate(thetas):
#             r2 = (carbon_params[key2]['r'],carbon_params[key2]['theta'],theta)
#             coupling.append(calc_C2C_coupling(r1,r2)/2.)
#             r2 = (carbon_params[key2]['r'],np.pi - carbon_params[key2]['theta'],theta)
#             coupling2.append(calc_C2C_coupling(r1,r2)/2.)
            
#         ax2.plot(thetas*180/np.pi,coupling,label='same')
#         ax2.axhline(C2C_coupling[key1[1] + key2[1]]['C'],0,360, color = 'r',label = 'msmt_pos')
#         ax2.axhline(-1.*C2C_coupling[key1[1] + key2[1]]['C'],0,360, color = 'r',label = 'msmt_neg')
#         ax2.plot(thetas*180/np.pi,coupling2,label='diff')
#         ax2.set_xlabel('phi',fontsize =fontsize)
#         ax2.set_ylabel('Coupling (Hz)', fontsize =fontsize)
#         plt.title('set_'+key1+'_vary_' +key2)
#         ax2.set_ylim(-4,4)
#         ax2.set_xlim(0,360)
#         # plt.legend()

#     key1 = 'C1'
# carbon_params[key1]['phi'] = 0
# carbon_params['C5']['phi'] = 70.5
# carbon_params['C3']['phi'] = 310.5

# carbon_params['C3']['theta'] = carbon_params['C3']['theta']
# r1 = (carbon_params[key1]['r'],carbon_params[key1]['theta'],0)
# key3 = 'C5'
# r3 = (carbon_params[key3]['r'],carbon_params[key3]['theta'],carbon_params[key3]['phi'])

# thetas = np.linspace(0,2*np.pi,201)
# fig = plt.figure()
# ax2 = fig.add_subplot(111)


'''
key1 = 'C1'
carbon_params[key1]['phi'] = 0
carbon_params['C5']['phi'] = 70.5
carbon_params['C3']['phi'] = 310.8
r1 = (carbon_params[key1]['r'],carbon_params[key1]['theta'],0)
key3 = 'C5'
r3 = (carbon_params[key3]['r'],carbon_params[key3]['theta'],70.5)

thetas = np.linspace(0,2*np.pi,201)
fig = plt.figure()
ax2 = fig.add_subplot(111)

for key2 in carbon_params.keys():
    if key2 != key1 and C2C_coupling[key1[1] + key2[1]]['C'] != None:
        if key2 != key3 and C2C_coupling[key3[1] + key2[1]]['C'] != None:
            fig = plt.figure()
            ax2 = fig.add_subplot(111)
            coupling = []
            coupling2 = []
            coupling3 = []
            coupling4 = []
            for ii, theta in enumerate(thetas):
                r2 = (carbon_params[key2]['r'],carbon_params[key2]['theta'],theta)
                coupling.append(calc_C2C_coupling(r1,r2)/2.)
                r2 = (carbon_params[key2]['r'],carbon_params[key2]['theta'],theta)
                coupling3.append(calc_C2C_coupling(r3,r2)/2.)
                r2 = (carbon_params[key2]['r'],np.pi - carbon_params[key2]['theta'],theta)
                coupling2.append(calc_C2C_coupling(r1,r2)/2.)
                r2 = (carbon_params[key2]['r'],np.pi - carbon_params[key2]['theta'],theta)
                coupling4.append(calc_C2C_coupling(r3,r2)/2.)
            ax2.plot(thetas*180/np.pi,coupling,label='sameC1',color='g')
            ax2.axhline(C2C_coupling[key1[1] + key2[1]]['C'],0,360, color = 'r',label = 'msmt_pos_C1')
            ax2.axhline(-1.*C2C_coupling[key1[1] + key2[1]]['C'],0,360, color = 'r',label = 'msmt_neg_C1')
            ax2.axhline(C2C_coupling[key3[1] + key2[1]]['C'],0,360, color = 'm',label = 'msmt_pos_C2')
            ax2.axhline(-1.*C2C_coupling[key3[1] + key2[1]]['C'],0,360, color = 'm',label = 'msmt_neg_C2')
            ax2.plot(thetas*180/np.pi,coupling2,label='diffC1',color='k')
            ax2.plot(thetas*180/np.pi,coupling3,label='sameC2',color='g',linestyle='dashed')
            ax2.plot(thetas*180/np.pi,coupling4,label='diffC2',color='k',linestyle='dashed')
            # ax2.axvline(x=311, ymin=-10, ymax=10)
            ax2.set_xlabel('phi')
            ax2.set_ylabel('Coupling (Hz)')
            plt.title('set_'+key1+'_vary_' +key2)
            ax2.set_ylim(-3,0)
            ax2.set_xlim(0,360)
            plt.legend(loc = 'center left',bbox_to_anchor = (1.0, 0.5))
'''

# carbon_params['C1']['phi'] = 0
# carbon_params['C5']['phi'] = 70.5
# carbon_params['C3']['phi'] = 310.8

# known_carbons = []
# known_keys = []
# for carbons in [1,3,5]:
#     key = 'C' + str(carbons)
#     known_keys.append(key)
#     known_carbons.append((carbon_params[key]['r'],carbon_params[key]['theta'],carbon_params[key]['phi']))

# thetas = np.linspace(0,2*np.pi,501)

# for key in carbon_params.keys():
#     if key not in known_keys:
#         fig = plt.figure()
#         ax2 = fig.add_subplot(111)
#         for rr, r2 in enumerate(known_carbons):
#             coupling1 = []
#             coupling2 = []
#             for theta in thetas:
#                 r1 = (carbon_params[key]['r'],carbon_params[key]['theta'],theta)
#                 coupling1.append(calc_C2C_coupling(r1,r2)/2.)
#                 r2 = (carbon_params[key]['r'],-1.*carbon_params[key]['theta'],theta)
#                 coupling2.append(calc_C2C_coupling(r1,r2)/2.)
#             ax2.plot(thetas*180/np.pi,coupling1,label='up' + known_keys[rr],color='g')
#             ax2.plot(thetas*180/np.pi,coupling2,label='down' + known_keys[rr],color='g')
#             ax2.axhline(C2C_coupling[key[1] + known_keys[rr][1]]['C'],0,360, color = 'r',label = 'msmt_pos_'+ known_keys[rr])
#             ax2.axhline(-1.*C2C_coupling[key[1] + known_keys[rr][1]]['C'],0,360, color = 'r',label = 'msmt_neg_' + known_keys[rr])
#         ax2.set_xlabel('Phi')
#         ax2.set_ylabel('Coupling (Hz)')
#         ax2.set_ylim(-4,4)
#         ax2.set_xlim(0,360)
#         plt.title('_vary_' +key1)
#         plt.legend(loc = 'center left',bbox_to_anchor = (1.0, 0.5))
#         plt.show()


plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15) 
fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X (nm)', fontsize = fontsize)
ax.set_ylabel('Y (nm)', fontsize = fontsize)
ax.set_zlabel('Z (nm)', fontsize = fontsize)
ax.scatter([0.], [0.], [0.], c='r', marker='o', label='NV',s=np.pi * (6 * np.ones((3,)))**2,linewidths=1)
# for carbon in [1,3,5]:
for carbon in [1,2,3,5,6]:
	plot_hf_circle(ax,carbon)

# for carbon in [3,5]:
#     plot_hf_circle2(ax,carbon)
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
plot_carbon(ax,1,0,posneg='pos')
plot_carbon(ax,2,1.35526965, posneg='pos')
plot_carbon(ax,3,4.57317311, posneg='neg')
plot_carbon(ax,5,1.2760995, posneg='pos')
plot_carbon(ax,6,3.88420822, posneg='neg')
ax.legend()

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
