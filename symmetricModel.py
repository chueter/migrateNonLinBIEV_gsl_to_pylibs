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
from scipy.integrate import quad # returns array where at ...[0] the value of the integral is saved
from scipy.integrate import romberg
from scipy.special import *

def integrand_consLaw_LorR(xpVal,xObsVal,yObsVal,Ps_i,Peclet):
	ypVal = Ps_i.__call__(xpVal)
	etaVal = eta(xpVal,xObsVal,ypVal,yObsVal)
	#return 2.*np.exp(-Peclet*( yObsVal - ypVal )) * kv(0,Peclet*etaVal)
	
	if (  abs(xpVal - xObsVal) < 1.e-7 ):
		return 0.
	else:
		return 2.*np.exp(-Peclet*( yObsVal - ypVal )) * kv(0,Peclet*etaVal)	
	
	
def integrate_consLaw__L(xObsVal,yObsVal,Ps_L,Peclet):
	#return romberg(integrand_consLaw_LorR, xProl_leftArray[0], xProl_leftArray[-1-nAdd], args=(xObsVal,yObsVal,Ps_L,Peclet),rtol=1.48e-04, divmax=10)
	return quad(integrand_consLaw_LorR, xProl_leftArray[0], xProl_leftArray[-1-nAdd], args=(xObsVal,yObsVal,Ps_L,Peclet))[0]
	
def integrate_consLaw__R(xObsVal,yObsVal,Ps_R,Peclet):
	#return romberg(integrand_consLaw_LorR, xProl_rightArray[nAdd], xProl_rightArray[-1], args=(xObsVal,yObsVal,Ps_R,Peclet),rtol=1.48e-04, divmax=10)
	return quad(integrand_consLaw_LorR, xProl_rightArray[nAdd], xProl_rightArray[-1], args=(xObsVal,yObsVal,Ps_R,Peclet))[0]
	
def eta(xp,xObs,yp,yObs):
	return ( (xObs-xp)**2 + (yObs-yp)**2   )**(0.5)
'''	
def quad_routine(func, a, b, x_list, w_list):
    c_1 = (b-a)/2.0
    c_2 = (b+a)/2.0
    eval_points = map(lambda x: c_1*x+c_2, x_list)
    func_evals = map(func, eval_points)
    return c_1 * sum(array(func_evals) * array(w_list))

def quad_gauss_7(func, a, b):
    x_gauss = [-0.949107912342759, -0.741531185599394, -0.405845151377397, 0, 0.405845151377397, 0.741531185599394, 0.949107912342759]
    w_gauss = array([0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277,0.129484966168870])
    return quad_routine(func,a,b,x_gauss, w_gauss)

def quad_kronrod_15(func, a, b):
    x_kr = [-0.991455371120813,-0.949107912342759, -0.864864423359769, -0.741531185599394, -0.586087235467691,-0.405845151377397, -0.207784955007898, 0.0, 0.207784955007898,0.405845151377397, 0.586087235467691, 0.741531185599394, 0.864864423359769, 0.949107912342759, 0.991455371120813]
    w_kr = [0.022935322010529, 0.063092092629979, 0.104790010322250, 0.140653259715525, 0.169004726639267, 0.190350578064785, 0.204432940075298, 0.209482141084728, 0.204432940075298, 0.190350578064785, 0.169004726639267, 0.140653259715525,  0.104790010322250, 0.063092092629979, 0.022935322010529]
    return quad_routine(func,a,b,x_kr, w_kr)

class Memorize(object):
    def __init__(self, func):
        self.func = func
        self.eval_points = {}
    def __call__(self, *args):
        if args not in self.eval_points:
            self.eval_points[args] = self.func(*args)
        return self.eval_points[args]

def quad(func,a,b):
    ### Output is the 15 point estimate; and the estimated error 
    func = Memorize(func) #  Memoize function to skip repeated function calls.
    g7 = quad_gauss_7(func,a,b)
    k15 = quad_kronrod_15(func,a,b)
    # I don't have much faith in this error estimate taken from wikipedia
    # without incorporating how it should scale with changing limits
    return [k15, (200*scipy.absolute(g7-k15))**1.5]
'''	
#################################################################################
# parameters
nx_L = 21
nx_R = 21

nAdd = 5
hx_L = 10./(nx_L-1)
hx_R = 10./(nx_R-1) 
interpolationFactor = 2
ip = interpolationFactor

kappa = np.zeros(nx_L+nx_R)
integral = np.zeros(nx_L+nx_R)

Peclet = 0.1
DeltaVal = ( Peclet * np.pi )**(0.5) * np.exp(Peclet) * erfc( Peclet**(0.5) )
sigma = 1.
angle_delta = 0.5*np.pi
slope_L = -1./np.tan(-angle_delta)
slope_R = -1./np.tan(angle_delta)
#angle_two_phi = 2./3. * np.pi

xLeftArrayList = [ -(nx_L-1-ix)*hx_L for ix in range(nx_L)]
xRightArrayList = [ (ix)*hx_R for ix in range(nx_R)]
xArrayList = xLeftArrayList + xRightArrayList
xSymm_array = np.asarray(xArrayList)

xAddList_L = [ 1.0, 2.0, 3.0, 4.0 , 5.0 ]
xAddList_R = [ -5.0, -4.0, -3.0, -2.0, -1.0]
#yAddList_L = ([ el*slope_L for el in xAddList_L ])
#yAddList_R = ([ el*slope_R for el in xAddList_R ])

xProlLeftArrayList = [ -(nx_L-1-ix)*hx_L for ix in range(nx_L+nAdd)] 
xProlRightArrayList = [ -nAdd*hx_R+(ix)*hx_R for ix in range(nx_R+nAdd)]

xProlLeftArrayListInterp = [ -((nx_L-1)*ip-ix)*hx_L/ip for ix in range((nx_L+nAdd-1)*ip+1)] 
xProlRightArrayListInterp = [ -nAdd*hx_R+(ix)*hx_R/ip for ix in range((nx_R+nAdd-1)*ip+1)]

xProl_leftArray = np.asarray(xProlLeftArrayList)
xProl_rightArray = np.asarray(xProlRightArrayList)

xProl_leftArrayInterp = np.asarray(xProlLeftArrayListInterp)
xProl_rightArrayInterp = np.asarray(xProlRightArrayListInterp)

cntCalls = 0

def residual(P):
	global integral, Peclet, DeltaVal, sigma, angle_delta, slope_L, slope_R#, angle_two_phi
	global xProl_leftArray, xProl_rightArray, xProl_leftArrayInterp, xProl_rightArrayInterp#, yAddList_L, yAddList_R
	global cntCalls
	cntCalls += 1
	eqsSystem = np.zeros(nx_L+nx_R)
	
	#print cntCalls, 'cntCalls'
	#print P
	#if cntCalls >=3:
			#exit(1)
	
	P_L = P[0:nx_L]
	P_R = P[nx_L::]
	
	sigma = P_L[-1]
	angle_delta = P_R[0]
	slope_L = -1./np.tan(-angle_delta)
	slope_R = -1./np.tan(angle_delta)
	print slope_L, slope_R
	
	P_L[-1] = 0.
	P_R[0] = 0.
	
	yAddList_L = ([ el*slope_L for el in xAddList_L ])
	yAddList_R = ([ el*slope_R for el in xAddList_R ])
	
	P_Llist = list(P_L)
	extP_L =np.asarray(P_Llist+yAddList_L)
	P_Rlist = list(P_R)
	extP_R =np.asarray(yAddList_R+P_Rlist)
	
	#print extP_R
	#print extP_L
	
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
		print 'iL, Delta, DeltaCalc', iL, DeltaVal, integral[iL]
		eqsSystem[iL] = DeltaVal + sigma*kappaObsVal - integral[iL]
		
	
	y_m1 = Ps_L.__call__(xProl_leftArray[-1-(nAdd+1)])
	y_p1 = Ps_R.__call__(xProl_rightArray[nAdd+1])
	
		
	eqsSystem[nx_L-1] = y_m1-slope_L*xProl_leftArray[-1-(nAdd+1)]
	eqsSystem[nx_L] = y_p1-slope_R*xProl_rightArray[nAdd+1]
	#exit(1)
	#print range(nx_R-1)
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
		print 'iR, Delta, DeltaCalc', iR, DeltaVal, integral[iR]
		eqsSystem[nx_L+1+iR] = DeltaVal + sigma*kappaObsVal - integral[iR]
		
	print eqsSystem 
	exit(1)
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













