#
#  peritecticPyCode.py
#  
#
#  Created by Claas Hueter on 5/9/12.
#  Copyright (c) 2012 __tinycodes__. All rights reserved.
#
'''
import numpy as np
from scipy.optimize import newton_krylov
from numpy import cosh, zeros_like, mgrid, zeros

# parameters
nx, ny = 75, 75
hx, hy = 1./(nx-1), 1./(ny-1)

P_left, P_right = 0, 0
P_top, P_bottom = 1, 0

def residual(P):
    d2x = zeros_like(P)
    d2y = zeros_like(P)

    d2x[1:-1] = (P[2:]   - 2*P[1:-1] + P[:-2]) / hx/hx
    d2x[0]    = (P[1]    - 2*P[0]    + P_left)/hx/hx
    d2x[-1]   = (P_right - 2*P[-1]   + P[-2])/hx/hx

    d2y[:,1:-1] = (P[:,2:] - 2*P[:,1:-1] + P[:,:-2])/hy/hy
    d2y[:,0]    = (P[:,1]  - 2*P[:,0]    + P_bottom)/hy/hy
    d2y[:,-1]   = (P_top   - 2*P[:,-1]   + P[:,-2])/hy/hy

    return d2x + d2y + 5*cosh(P).mean()**2

# solve
guess = zeros((nx, ny), float)
sol = newton_krylov(residual, guess, verbose=1)
#sol = broyden2(residual, guess, max_rank=50, verbose=1)
#sol = anderson(residual, guess, M=10, verbose=1)
print 'Residual', abs(residual(sol)).max()

# visualize
import matplotlib.pyplot as plt
x, y = mgrid[0:1:(nx*1j), 0:1:(ny*1j)]
plt.pcolor(x, y, sol)
plt.colorbar()
plt.show()
'''
import numpy as np
from scipy.optimize import newton_krylov
from numpy import cosh, zeros_like, mgrid, zeros
import pylab
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quadrature # returns array where at ...[0] the value of the integral is saved
from scipy.integrate import romberg
from scipy.special import *

def integrand_consLaw_LorR(xpVal,xObsVal,yObsVal,Ps_i,Peclet):
	ypVal = Ps_i.__call__(xpVal)
	#ypVal = -0.5*xpVal**2
	#yObsVal = Ps_i.__call__(xObsVal)
	etaVal = eta(xpVal,xObsVal,ypVal,yObsVal)
	
	if (  abs(xpVal - xObsVal) < 1.e-7 ):
		return 0.
	else:
		return 2.*np.exp(-Peclet*( yObsVal - ypVal )) * kv(0,Peclet*etaVal)	
	
def integrate_consLaw__L(xObsVal,yObsVal,Ps_L,Peclet):
	#return romberg(integrand_consLaw_LorR, xProl_leftArray[0], xProl_leftArray[-1-nAdd], args=(xObsVal,yObsVal,Ps_L,Peclet),rtol=1.48e-03, divmax=17)
	return quadrature(integrand_consLaw_LorR, xProl_leftArray[0], xProl_leftArray[-1-nAdd], args=(xObsVal,yObsVal,Ps_L,Peclet))[0]
	
def integrate_consLaw__R(xObsVal,yObsVal,Ps_R,Peclet):
	#return romberg(integrand_consLaw_LorR, xProl_rightArray[nAdd], xProl_rightArray[-1], args=(xObsVal,yObsVal,Ps_R,Peclet),rtol=1.48e-03, divmax=17)
	return quadrature(integrand_consLaw_LorR, xProl_rightArray[0], xProl_rightArray[-1-nAdd], args=(xObsVal,yObsVal,Ps_R,Peclet))[0]
	
def eta(xp,xObs,yp,yObs):
	return ( (xObs-xp)**2 + (yObs-yp)**2   )**(0.5)
	
	
#################################################################################
# parameters
nx_L = 41
nx_R = 41

nAdd = 5
hx_L = 10./(nx_L-1)
hx_R = 10./(nx_R-1) 
interpolationFactor = 5
ip = interpolationFactor

kappa = np.zeros(nx_L+nx_R)
integral = np.zeros(nx_L+nx_R)

Peclet = 0.1
DeltaVal = ( Peclet * np.pi )**(0.5) * np.exp(Peclet) * erfc( Peclet**(0.5) )
sigma = 1.
angle_delta = 0.25*np.pi
slope_L = -1./np.tan(-angle_delta)
slope_R = -1./np.tan(angle_delta)
#angle_two_phi = 2./3. * np.pi

xLeftArrayList = [ -(nx_L-1-ix)*hx_L for ix in range(nx_L)]
xRightArrayList = [ (ix)*hx_R for ix in range(nx_R)]
xArrayList = xLeftArrayList + xRightArrayList
xSymm_array = np.asarray(xArrayList)

xAddList_L = [ 1.0, 2.0, 3.0, 4.0 , 5.0 ]
xAddList_R = [ -5.0, -4.0, -3.0, -2.0, -1.0]
yAddList_L = ([ el*slope_L for el in xAddList_L ])
yAddList_R = ([ el*slope_R for el in xAddList_R ])

xProlLeftArrayList = [ -(nx_L-1-ix)*hx_L for ix in range(nx_L+nAdd)] 
xProlRightArrayList = [ -nAdd*hx_R+(ix)*hx_R for ix in range(nx_R+nAdd)]

xProlLeftArrayListInterp = [ -((nx_L-1)*ip-ix)*hx_L/ip for ix in range((nx_L+nAdd-1)*ip+1)] 
xProlRightArrayListInterp = [ -nAdd*hx_R+(ix)*hx_R/ip for ix in range((nx_R+nAdd-1)*ip+1)]

xProl_leftArray = np.asarray(xProlLeftArrayList)
xProl_rightArray = np.asarray(xProlRightArrayList)

xProl_leftArrayInterp = np.asarray(xProlLeftArrayListInterp)
xProl_rightArrayInterp = np.asarray(xProlRightArrayListInterp)

def residual(P):
	global integral, Peclet, DeltaVal, sigma, angle_delta, slope_L, slope_R#, angle_two_phi
	global xProl_leftArray, xProl_rightArray, xProl_leftArrayInterp, xProl_rightArrayInterp, yAddList_L, yAddList_R
	
	eqsSystem = np.zeros(nx_L+nx_R)
	
	P_L = P[0:nx_L]
	P_R = P[nx_L::]
	
	sigma = P_L[-1]
	angle_delta = 0.25*np.pi#P_R[0]
	slope_L = -1./np.tan(-angle_delta)
	slope_R = -1./np.tan(angle_delta)

	P_L[-1] = 0.
	P_R[0] = 0.
	
	P_Llist = list(P_L)
	extP_L =np.asarray(P_Llist+yAddList_L)
	P_Rlist = list(P_R)
	extP_R =np.asarray(yAddList_R+P_Rlist)
	
	Ps_L = UnivariateSpline(xProl_leftArray, extP_L, s=0)
	Ps_R = UnivariateSpline(xProl_rightArray, extP_R, s=0)
	P_Linterp = Ps_L(xProl_leftArrayInterp)
	P_Rinterp = Ps_R(xProl_rightArrayInterp)
	
	for iL in range(nx_L-1):
		xObsVal = xProl_leftArray[iL]
		yObsVal = Ps_L.__call__(xProl_leftArray[iL],0)
		dydxObsVal = Ps_L.__call__(xProl_leftArray[iL],1)
		d2ydx2ObsVal = Ps_L.__call__(xProl_leftArray[iL],2)
		dsObsVal = np.sqrt( (1.+dydxObsVal*dydxObsVal) )
		kappaObsVal = d2ydx2ObsVal/(dsObsVal*dsObsVal*dsObsVal) # this way the error on the curvature close to the tip is heavily dependent on the 
																#discretization close to the tip
		integral_L = integrate_consLaw__L(xObsVal, yObsVal,Ps_L,Peclet)
		integral_R = integrate_consLaw__R(xObsVal, yObsVal,Ps_R,Peclet)
		integral[iL] = Peclet/(2.*np.pi) * (integral_L + integral_R)
		
		eqsSystem[iL] = DeltaVal + sigma*kappaObsVal - integral[iL]
		
	
	y_m1 = Ps_L.__call__(xProl_leftArray[-1-(nAdd+1)])
	y_p1 = Ps_R.__call__(xProl_rightArray[nAdd+1])
	
		
	eqsSystem[nx_L-1] = y_m1-slope_L*xProl_leftArray[-1-(nAdd+1)]
	eqsSystem[nx_L] = y_p1-slope_R*xProl_rightArray[nAdd+1]
	
	for iR in range(nx_R-1):
		xObsVal = xProl_rightArray[nAdd+1+iR]
		yObsVal = Ps_R.__call__(xProl_rightArray[nAdd+1+iR])
		dydxObsVal = Ps_R.__call__(xProl_rightArray[nAdd+1+iR],1)
		d2ydx2ObsVal = Ps_R.__call__(xProl_rightArray[nAdd+1+iR],2)
		dsObsVal = np.sqrt( (1.+dydxObsVal*dydxObsVal) )
		kappaObsVal = d2ydx2ObsVal/(dsObsVal*dsObsVal*dsObsVal) # this way the error on the curvature close to the tip is heavily dependent on the 
																#discretization close to the tip
		integral_L = integrate_consLaw__L(xObsVal, yObsVal,Ps_L,Peclet)
		integral_R = integrate_consLaw__R(xObsVal, yObsVal,Ps_R,Peclet)
		integral[iR] = Peclet/(2.*np.pi) * (integral_L + integral_R)
		eqsSystem[nx_L+1+iR] = DeltaVal + sigma*kappaObsVal - integral[iR]
	
	return eqsSystem
	
guess = np.asarray([ -0.5*el**2 for el in xArrayList ])
guess[nx_L-1] = 1.
guess[nx_L] = angle_delta

sol = newton_krylov(residual, guess, verbose=1)
#sol = broyden2(residual, guess, max_rank=50, verbose=1)
#sol = anderson(residual, guess, M=10, verbose=1)
print 'Residual', abs(residual(sol)).max()

solList = [el for el in list(sol)]
sol_array = np.asarray( solList )

pylab.plot(x_array, sol_array )
pylab.xlabel(' ')
pylab.ylabel(' ')
pylab.title(' ')
pylab.grid(True)
pylab.savefig('testMeiron')

pylab.show()

exit(1)













