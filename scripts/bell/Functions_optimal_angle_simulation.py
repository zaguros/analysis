## IMPORT MODULES ##
import numpy as np
import scipy.optimize as opt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import BoundaryNorm
import time

## Define Rotation Operator ##
def rotation (x, y, Rho):
	a = x*np.pi/180.
	b = y*np.pi/180.
	ra = np.array ([[np.cos(a), -np.sin(a)],[np.sin(a), np.cos(a)]])
	rb = np.array ([[np.cos(b), -np.sin(b)],[np.sin(b), np.cos(b)]])
	rot = np.kron(ra, rb)
	rotD = rot.conj().T
	rho1 = np.dot(rot, np.dot(Rho, rotD))
	P = np.array ([rho1[0,0], rho1[1,1], rho1[2,2], rho1[3,3]])
	P = P.T
	return P

## Define readout Matrix ##
def Readout(P,F,G):
    F = np.array([[F[0], 1-F[1]],[1-F[0], F[1]]])
    G = np.array([[G[0], 1-G[1]],[1-G[0], G[1]]])
    RO = np.kron(F,G)
    X = np.dot(RO,P)
    return X

 ## Define E(x,y) Function ##
def E(N):
    E_a_b =  (N[0]-N[1]-N[2]+N[3])/(0.+N[0]+N[1]+N[2]+N[3])
    return E_a_b
######################################################################

## Define CHSH inequality ##
def CHSH(Ang, Rho,F,G):
    p = Readout(rotation(0.,Ang[0],Rho),F,G)
    q = Readout(rotation(Ang[1],Ang[0],Rho),F,G)
    r = Readout(rotation(0.,Ang[2],Rho),F,G)
    s = Readout(rotation(Ang[1],Ang[2],Rho),F,G)


    CHSH = abs(E(p) + E(q) + E(r) - E(s))
    return CHSH

######################################################################

## Invert CHSH ##
def Invert(Ang, Rho,F,G):
    a = - CHSH(Ang,Rho,F,G)
    return a

def Include_errors(Rho_pure,Error_dist,Error_dark,D,dc):
    a = (1-dc) *((1-D)*Rho_pure + D * Error_dist) + dc * Error_dark
    return a