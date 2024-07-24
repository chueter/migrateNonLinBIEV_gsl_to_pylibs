#
#  exampleOfIntegroDiffEquation.py
#  
#
#  Created by Claas Hueter on 5/9/12.
#  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
#
'''
np.trapz([1,2,3])
4.0
'''



import numpy as np
from scipy.optimize import newton_krylov
from numpy import cosh, zeros_like, mgrid, zeros
import pylab
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
from scipy.special import *



'''
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
########################################################################
'''
>>> from scipy.integrate import quad
>>> def integrand(t,n,x):
...     return exp(-x*t) / t**n
>>> def expint(n,x):
...     return quad(integrand, 1, Inf, args=(n, x))[0]
>>> vec_expint = vectorize(expint)
>>> vec_expint(3,arange(1.0,4.0,0.5))
array([ 0.1097,  0.0567,  0.0301,  0.0163,  0.0089,  0.0049])
'''


def integrandTest1(t,n,x):
	return np.exp(-x*t) / t**n
	
def expint(n,x):
	return quad(integrandTest1, 1, np.Inf, args=(n, x))[0] # 1st element of returned tuple is val of integral, second is upper bound for error

vec_expint = np.vectorize(expint)

print '#########################'
print '#########################'	
print 'test1: exponential integral for n = 3, x = np.arange(1.0,4.0,0.5)', vec_expint(3,np.arange(1.0,4.0,0.5))	
print '#########################'
print '#########################'
########################################################################

#################################################################################
def integrand(xc_val,P_val):
	return P_val-(xc_val) 

def integrand_xObs(xc_val,xcObs,P_val):
	return P_val-xc_val+xcObs

def integrand_xObs_interp(xc_val,xcObs,P_val,Ps):
	return ( (P_val-Ps.__call__(xcObs))**2+(xc_val-xcObs)**2 )**(0.5)

def intVal(P_val):
	return quad(integrand,xc[0], xc[-1], args=(P_val))[0]

def intObsVal(xcObs,P_val):
	return quad(integrand_xObs,xc[0], xc[-1], args=(xcObs,P_val))[0]

def intObsInterpVal(xcObs,P_val,Ps):
	return quad(integrand_xObs_interp,xc[0], xc[-1], args=(xcObs,P_val,Ps))[0]

vec_intVal = np.vectorize(intVal)
vec_intObsVal = np.vectorize(intObsVal)
#################################################################################


# parameters
nx = 111
hx = 20./(nx)
interpolationFactor = 2

ip = interpolationFactor

P_left, P_right = 0, 0

#x_array = np.arange(0.0,nx*hx,hx)
xSymm_array = np.arange(-0.5*nx*hx,0.5*nx*hx,hx)

#xnew = np.arange(0.0,nx*hx,hx/interpolationFactor)
xSymm_new = np.arange(-0.5*nx*hx,0.5*nx*hx,hx/interpolationFactor)

#xc = x_array
xcSymm = xSymm_array

DeltaVal = 0.01
Peclet = (DeltaVal)**2 * 1./np.pi


def integrand_ivantsov(xpVal,xObsVal,Ps,peclet):
	#ypVal = Ps.__call__(xpVal)
	ypVal = -0.5*xpVal**2
	yObsVal = Ps.__call__(xObsVal)
	etaVal = eta(xpVal,xObsVal,ypVal,yObsVal)
	
	#### integrand is definitely correct, checked at different places all contributions 
	#if (abs(xObsVal) <= 3.):
		#print xpVal, xObsVal, ypVal, yObsVal, etaVal, peclet , kv(0,peclet*etaVal), np.exp(-peclet*( yObsVal - ypVal )),'in integrand_ivantsov'
		#exit(1)
	#else:
		#pass
	
	if (  abs(xpVal - xObsVal) < 1.e-7 ):
		return 0.
	else:
		return np.exp(-peclet*( yObsVal - ypVal )) * kv(0,peclet*etaVal)
	
	

def integrand_ivantsovVal(xObsVal,Ps,peclet):
	#print xcSymm[0], xcSymm[-1], peclet, xObsVal ,quad(integrand_ivantsov,xcSymm[0], xcSymm[-1], args=(xObsVal,Ps,peclet))[0], 'in integrand_ivantsovVal'
	#exit(1)
	return quad(integrand_ivantsov,xcSymm[0], xcSymm[-1], args=(xObsVal,Ps,peclet))[0]

def eta(xp,xObs,yp,yObs):
	return ( (xObs-xp)**2 + (yObs-yp)**2   )**(0.5)


def residual(P):
	d2x = zeros_like(P)
	integral = zeros_like(P)
	
	
	Peclet = 0.1
	DeltaVal = ( Peclet * np.pi )**(0.5) * np.exp(Peclet) * erfc( Peclet**(0.5) )
									
	Ps = UnivariateSpline(xSymm_array, P, s=0)
	P_interp = Ps(xSymm_new)
	
	Pnew = Ps
	
	#d2x[1:-1] = (Pnew.__call__(xc[2:])		- 2*Pnew.__call__(xc[1:-1])		+ Pnew.__call__(xc[:-2]) )	/ hx/hx
	#d2x[0]    = (Pnew.__call__(xc[1])		- 2*Pnew.__call__(xc[0])		+ P_left)					/ hx/hx
	#d2x[-1]   = (P_right					- 2*Pnew.__call__(xc[-1])		+ Pnew.__call__(xc[-2]))	/ hx/hx
	
	
	#for i in range(len(xc)-2):
		#integral[i+1]	= integrand_ivantsovVal(xc[i],Ps,Peclet)
	#integral[0]		= integrand_ivantsovVal(xc[1],Ps,Peclet)
	#integral[-1]	= integrand_ivantsovVal(xc[-1],Ps,Peclet)
	
	for i in range(len(xcSymm)-2):
		integral[i+1]	= integrand_ivantsovVal(xcSymm[i+1],Ps,Peclet)
		
	integral[0]		= integrand_ivantsovVal(xcSymm[1],Ps,Peclet)
	integral[-1]	= integrand_ivantsovVal(xcSymm[-1],Ps,Peclet)

	
	
	
	print DeltaVal, Peclet/(np.pi) * integral[5]
	
	return DeltaVal - Peclet/(np.pi) * integral#vec_intVal(P) + vec_intObsVal(xc,P) + vec_intObsInterpVal(xc,P,Pnew)


	
# solve
#guess = zeros(nx, float)	
#print type(guess), ' = type guess'
guess = np.asarray([ -0.5*el**2 for el in xcSymm ])
#print guess


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
pylab.savefig('testIntEq')

pylab.show()

exit(1)

########################################################################


x = np.arange(0,2*np.pi+np.pi/4,2*np.pi/8)
y = np.sin(x)
tck = interpolate.splrep(x,y,s=0)
xnew = np.arange(0,2*np.pi,np.pi/50)
ynew = interpolate.splev(xnew,tck,der=0)
'''
plt.figure()
plt.plot(x,y,'x',xnew,ynew,xnew,np.sin(xnew),x,y,'b')
plt.legend(['Linear','Cubic Spline', 'True'])
plt.axis([-0.05,6.33,-1.05,1.05])
plt.title('Cubic-spline interpolation')
plt.show()
'''
def integ(x,tck,constant=-1):
	x = np.atleast_1d(x)
	out = np.zeros(x.shape, dtype=x.dtype)
	for n in xrange(len(out)):
		out[n] = interpolate.splint(0,x[n],tck)
	out += constant
	return out

yint = integ(xnew,tck)
plt.figure()
plt.plot(xnew,yint,xnew,-np.cos(xnew),'--')
plt.legend(['Cubic Spline', 'True'])
plt.axis([-0.05,6.33,-1.05,1.05])
plt.title('Integral estimation from spline')
plt.show()

'''
def func2(x):
    out = [x[0]*np.cos(x[1]) - 4]
    out.append(x[1]*x[0] - x[1] - 5)
    return out

from scipy.optimize import fsolve
x02 = fsolve(func2, [0, 0.5])
print x02
'''

'''
#########################################################################

#for (int i = 1; i< numberGridPoints;i++)
xObs = xGridG1[i+addGridPoints];	
yObs = gsl_spline_eval(current_splineG1,xObs,current_accelG1);
dydx = gsl_spline_eval_deriv(current_splineG1,xObs,current_accelG1);
d2ydx2 = gsl_spline_eval_deriv2(current_splineG1,xObs,current_accelG1);
kappa = gsl_spline_eval(kappaRHP_splineG1,xObs,kappaRHP_accelG1);

gsl_integration_qng(&integrateExactDiffKernel,0.00/*xGridG1[addGridPoints+1]*/, tailFactor*xGridG1[addGridPoints+numberGridPoints-2],0,1e-7,&Integral,&error,&n_evaluations);
rightIntegral_conservationLaw = Integral;
conservationLawIntegral = rightIntegral_conservationLaw
gsl_integration_qng(&integrateExactDiffKernel,/*xGridG1[addGridPoints+1]*/0.00, tailFactor*xGridG1[addGridPoints+numberGridPoints-2],0,1e-7,&Integral,&error,&n_evaluations);
rightIntegral_locEq = Integral;
locEqIntegral = rightIntegral_locEq;
gsl_integration_qng(&integrateExactDiffKernel, 0., yGridG12[addGridPoints+numberGridPoints_y-2],0,1e-7,&l1l2Integral,&error,&n_evaluations)
Integral = 1.*conservationLawIntegral + 1.*locEqIntegral + l1l2Integral

equation[i] = 0.5*(deltaG1-kappa*sigma)-(pG1/(2.*PI))*Integral


#########################################################################

equation[numberGridPoints-1] = B_N_R_end;
equation[numberGridPoints] = x_p1-slopeG12*yGridG12[addGridPoints+1];


#########################################################################

#for (int i=1; i<numberGridPoints_y-1;i++)
yObs = yGridG12[i+addGridPoints];	
xObs = gsl_spline_eval(current_splineG12,yObs,current_accelG12);
dxdy = gsl_spline_eval_deriv(current_splineG12,yObs,current_accelG12);
d2xdy2 = gsl_spline_eval_deriv2(current_splineG12,yObs,current_accelG12);
kappa = -d2xdy2/pow((1+dxdy*dxdy),1.5);
kappa = gsl_spline_eval(kappa_splineG12,yObs,kappa_accelG12);

gsl_integration_qng(&integrateExactDiffKernel,/*xGridG1[addGridPoints+1]*/0.00, tailFactor*xGridG1[addGridPoints+numberGridPoints-2],0,1e-7,&Integral,&error,&n_evaluations);
rightIntegral_conservationLaw = Integral;
conservationLawIntegral  =  rightIntegral_conservationLaw;
gsl_integration_qng(&integrateExactDiffKernel,0.00/*xGridG1[addGridPoints+1]*/, tailFactor*xGridG1[addGridPoints+numberGridPoints-2],0,1e-7,&Integral,&error,&n_evaluations);
rightIntegral_locEq = Integral;
locEqIntegral =  rightIntegral_locEq;
gsl_integration_qng(&integrateExactDiffKernel, 0., yGridG12[addGridPoints+numberGridPoints_y-2],0,1e-7,&l1l2Integral,&error,&n_evaluations);
Integral = 1.*conservationLawIntegral + 1.*locEqIntegral + l1l2Integral;

equation[i+numberGridPoints] = -0.5*DcLDelta_DcLGamma*kappa*sigma-(pG1/(2.*PI))*Integral; // 1.0 here comes from assumption that d_12 = d_1S !


#########################################################################
'''







#########################################################################






#########################################################################






#########################################################################







#########################################################################


















