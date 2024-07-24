#
#  fabiansThesisCalculations
#  
#
#  Created by Claas Hueter on 5/9/12.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#
from sys import argv, exit
import numpy as np
import pylab

################################################################################
################################################################################
################################################################################


"""
Example: simple line plot.
Show how to make and save a simple line plot with labels, title and grid
"""
'''
import numpy
import pylab

t = numpy.arange(0.0, 1.0+0.01, 0.01)
s = numpy.cos(2*2*numpy.pi*t)
pylab.plot(t, s)

pylab.xlabel('time (s)')
pylab.ylabel('voltage (mV)')
pylab.title('About as simple as it gets, folks')
pylab.grid(True)
pylab.savefig('simple_plot')

pylab.show()
'''

'''
from numpy import *
def f(i,j):
	return i**2 + j**2

fromfunction(f, (3,3)) # evaluate functiom for all combinations of indices [0,1,2]x[0,1,2]
array([[0, 1, 4],
       [1, 2, 5],
       [4, 5, 8]])

'''





################################################################################
################################################################################
################################################################################



"""
x^2 + ax + b = 0  (or ax^2 + bx + c = 0)
By substituting x = y-t and t = a/2, the equation reduces to 
    y^2 + (b-t^2) = 0 
which has easy solution
    y = +/- sqrt(t^2-b)
"""
def quadratic(a, b, c=None):
    import math, cmath
    if c:		# (ax^2 + bx + c = 0)
	a, b = b / float(a), c / float(a)
    t = a / 2.0
    r = t**2 - b
    if r >= 0:		# real roots
	y1 = math.sqrt(r)
    else:		# complex roots
	y1 = cmath.sqrt(r)
    y2 = -y1
    return y1 - t, y2 - t
	
	
	
"""
x^3 + ax^2 + bx + c = 0  (or ax^3 + bx^2 + cx + d = 0)
With substitution x = y-t and t = a/3, the cubic equation reduces to    
    y^3 + py + q = 0,
where p = b-3t^2 and q = c-bt+2t^3.  Then, one real root y1 = u+v can
be determined by solving 
    w^2 + qw - (p/3)^3 = 0
where w = u^3, v^3.  From Vieta's theorem,
    y1 + y2 + y3 = 0
    y1 y2 + y1 y3 + y2 y3 = p
    y1 y2 y3 = -q,
the other two (real or complex) roots can be obtained by solving
    y^2 + (y1)y + (p+y1^2) = 0
"""
def cubic(a, b, c, d=None):
    from math import cos
    if d:			# (ax^3 + bx^2 + cx + d = 0)
	a, b, c = b / float(a), c / float(a), d / float(a)
    t = a / 3.0
    p, q = b - 3 * t**2, c - b * t + 2 * t**3
    u, v = quadratic(q, -(p/3.0)**3)
    if type(u) == type(0j):	# complex cubic root
	r, w = polar(u.real, u.imag)
	y1 = 2 * cbrt(r) * cos(w / 3.0)
    else:			# real root
        y1 = cbrt(u) + cbrt(v)
    y2, y3 = quadratic(y1, p + y1**2)
    return y1 - t, y2 - t, y3 - t
	
	
"""
Ubiquitous Newton-Raphson algorithm for solving
    f(x) = 0
where a root is repeatedly estimated by
    x = x - f(x)/f'(x)
until |dx|/(1+|x|) < TOL is achieved.  This termination condition is a
compromise between 
    |dx| < TOL,  if x is small
    |dx|/|x| < TOL,  if x is large
"""
def newton(func, funcd, x, TOL=1e-6):	# f(x)=func(x), f'(x)=funcd(x)
    f, fd = func(x), funcd(x)
    count = 0
    while 1:
	dx = f / float(fd)
	if abs(dx) < TOL * (1 + abs(x)): return x - dx
	x = x - dx
	f, fd = func(x), funcd(x)
	count = count + 1
	print "newton(%d): x=%s, f(x)=%s" % (count, x, f)
	
	
	
"""
Similar to Newton's method, but the derivative is estimated by divided
difference using only function calls.  A root is estimated by
    x = x - f(x) (x - oldx)/(f(x) - f(oldx))
where oldx = x[i-1] and x = x[i].
"""
def secant(func, oldx, x, TOL=1e-6):	# f(x)=func(x)
    oldf, f = func(oldx), func(x)
    if (abs(f) > abs(oldf)):		# swap so that f(x) is closer to 0
	oldx, x = x, oldx
	oldf, f = f, oldf
    count = 0
    while 1:
	dx = f * (x - oldx) / float(f - oldf)
	if abs(dx) < TOL * (1 + abs(x)): return x - dx
	oldx, x = x, x - dx
	oldf, f = f, func(x)
	count = count + 1
	print "secant(%d): x=%s, f(x)=%s" % (count, x, f)
	
	
"""
Closed Simpson's rule for 
    \int_a^b f(x) dx
Divide [a,b] iteratively into h, h/2, h/4, h/8, ... step sizes; and,
for each step size, evaluate f(x) at a+h, a+3h, a+5h, a+7h, ..., b-3h,
b-h, noting that other points have already been sampled.

At each iteration step, data are sampled only where necessary so that
the total data is represented by adding sampled points from all
previous steps:
    step 1:	h	a---------------b
    step 2:	h/2 	a-------^-------b
    step 3:	h/4	a---^-------^---b
    step 4:	h/8	a-^---^---^---^-b
    total:		a-^-^-^-^-^-^-^-b
So, for step size of h/n, there are n intervals, and the data are
sampled at the boundaries including the 2 end points.

If old = Trapezoid formula for an old step size 2h, then Trapezoid
formula for the new step size h is obtained by 
    new = old/2 + h{f(a+h) + f(a+3h) + f(a+5h) + f(a+7h) +...+ f(b-3h)
	+ f(b-h)}
Also, Simpson formula for the new step size h is given by
    simpson = (4 new - old)/3
"""
def closedpoints(func, a, b, TOL=1e-6):		# f(x)=func(x)
    h = b - a
    old2 = old = h * (func(a) + func(b)) / 2.0
    count = 0
    while 1:
	h = h / 2.0
	x, sum = a + h, 0
	while x < b:
	    sum = sum + func(x)
	    x = x + 2 * h
	new = old / 2.0 + h * sum
	new2 = (4 * new - old) / 3.0
	if abs(new2 - old2) < TOL * (1 + abs(old2)): return new2
	old = new	# Trapezoid
	old2 = new2	# Simpson
	count = count + 1
	print 'closedpoints(%d): Trapezoid=%s, Simpson=%s' % (count, new, new2)
	
		
## Get the first eigenvalues with Power iterations.
def eig_power(A, tol=1.e-6):
	"""
	Computes the eigenvector with largest eigenvalue
	of A. Returns the pair (eigval, eigvec).
	"""
	assert A.ndim == 2
	b = NP.ones((A.shape[1], 1), NP.float32)
	while True:
		Ab = NP.dot(A, b)
		b_old = b.copy()
		b = Ab/norm(Ab)
		if (norm(b-b_old) / norm(b_old)) < tol:
			break
	eigval = NP.dot(NP.transpose(b), NP.dot(A, b))
	eigval /= norm2(b)
	return eigval, b
	
	
# -*- coding: utf-8 -*-
# Talbot suggested that the Bromwich line be deformed into a contour that begins
# and ends in the left half plane, i.e., z to -inf at both ends.
# Due to the exponential factor the integrand decays rapidly
# on such a contour. In such situations the trapezoidal rule converge
# extraordinarily rapidly.
# For example here we compute the inverse transform of F(s) = 1/(s+1) at t = 1
#
# >>> error = Talbot(1,24)-exp(-1)
# >>> error
#   (3.3306690738754696e-015+0j)
#
# Talbot method is very powerful here we see an error of 3.3e-015
# with only 24 function evaluations
#
# Created by Fernando Damian Nieuwveldt      
# email:fdnieuwveldt@gmail.com
# Date : 25 October 2009
#
# Reference
# L.N.Trefethen, J.A.C.Weideman, and T.Schmelzer. Talbot quadratures
# and rational approximations. BIT. Numerical Mathematics,
# 46(3):653 670, 2006.
from cmath import *
def Talbot(t):
    N = 24
#   Initiate the stepsize
    h = 2*pi/N;
#   Shift contour to the right in case there is a pole on the positive real axis : Note the contour will
#   not be optimal since it was originally devoloped for function with
#   singularities on the negative real axis
#   For example take F(s) = 1/(s-1), it has a pole at s = 1, the contour needs to be shifted with one
#   unit, i.e shift  = 1. But in the test example no shifting is necessary
    shift = 0.0;
    ans =   0.0;
    
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
        
#   The for loop is evaluating the Laplace inversion at each point theta which is based on the trapezoidal   rule
    for k in range(0,N):
        theta = -pi + (k+1./2)*h;
        z = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans = ans + exp(z*t)*F(z)*dz;
        
    return ((h/(2j*pi))*ans).real        
      

#   Here is the Laplace function to be inverted, should be changed manually        
def F(s):
    return 1.0/(s+1.0)

#   Define the trig functions cot(phi) and csc(phi)
def cot(phi):
    return 1.0/tan(phi)

def csc(phi):
    return 1.0/sin(phi)
       	
################################################################################
################################################################################
################################################################################



def beta(d0_div_delta, phi, Delta):
	return 2*d0_div_delta*phi/Delta

def alpha(DeltaPM,Delta):
	return DeltaPM/Delta

def discriminantOfMuThirdRoot(alphaVal, betaVal ):
	return (0.25*np.pi**(2./3.) * (betaVal-alphaVal-1)**2 + np.pi**(2./3.)*betaVal )**(1./2.)

def muThirdRoot(discriminantVal,alphaVal, betaVal):
	return discriminantVal - 0.5 * np.pi**(1./3.) * (betaVal-alphaVal-1)

def s_star(mu):
	return (2./(np.pi*mu**2))**(1./3.)

def transformOf_dxdy(alphaVal,betaVal,muVal,s):
	#return (muVal-1./s-(alphaVal/(s+betaVal/muVal)))/(muVal*s-(2./np.pi)**(0.5) * s**(-0.5))
	return (muVal*s -1-(alphaVal*s/(s+betaVal/muVal)))/(muVal*s**2-(2./np.pi)**(0.5) * s**(0.5))
	'''
	if abs(muVal*s**2-(2./np.pi)**(0.5) * s**(0.5)) <= e-6:
		shifted_s = s+0.0001
		X = (muVal*s -1-(alphaVal*s/(s+betaVal/muVal)))/ (muVal*shifted_s**2-(2./np.pi)**(0.5) * shifted_s**(0.5))  
	else:
		X = (muVal*s -1-(alphaVal*s/(s+betaVal/muVal)))/(muVal*s**2-(2./np.pi)**(0.5) * s**(0.5))
	return X
	'''
################################################################################	

def Talbot_for_X(mXi,N,alphaVal,betaVal,muVal):
#   Initiate the stepsize
    h = 2*pi/N;
#   Shift contour to the right in case there is a pole on the positive real axis : Note the contour will
#   not be optimal since it was originally devoloped for function with
#   singularities on the negative real axis
#   For example take F(s) = 1/(s-1), it has a pole at s = 1, the contour needs to be shifted with one
#   unit, i.e shift  = 1. But in the test example no shifting is necessary
    shift = muVal+1.
    ans =   0.0;
    
    if mXi < 0:
        print "ERROR:   Inverse transform can not be calculated for -xi < 0"
        return ("Error");
        
#   The for loop is evaluating the Laplace inversion at each point theta which is based on the trapezoidal   rule
    for k in range(0,N):
        theta = -pi + (k+1./2)*h;
        z = shift + N/mXi*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz = N/mXi*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans = ans + exp(z*mXi)*transformOf_dxdy(alphaVal,betaVal,muVal,z)*dz;
        
    return ((h/(2j*pi))*ans).real        
	
################################################################################
################################################################################
################################################################################


			
'''			
betaVal = beta(100., 0.0001, 0.01)
alphaVal = alpha(0.01,0.01)
discVal = discriminantOfMuThirdRoot(alphaVal, betaVal )
mu1o3Val = muThirdRoot(discVal,alphaVal, betaVal)
print mu1o3Val**3		
'''
################################################################################

'''
Delta_array = np.arange(0.005, 0.25, 0.005)
mu1o3Val_array_DeltaPM0p01 = muThirdRoot(discriminantOfMuThirdRoot(alpha(0.01,Delta_array), beta(100., 0.1, Delta_array) ),alpha(0.01,Delta_array), beta(100., 0.1, Delta_array))
mu1o3Val_array_DeltaPM0p05 = muThirdRoot(discriminantOfMuThirdRoot(alpha(0.01,Delta_array), beta(100., 0.0001, Delta_array) ),alpha(0.05,Delta_array), beta(100., 0.0001, Delta_array))
#pylab.plot(Delta_array,mu1o3Val_array_DeltaPM0p01**3,Delta_array,mu1o3Val_array_DeltaPM0p05**3)
pylab.plot(Delta_array,mu1o3Val_array_DeltaPM0p01**3)
pylab.xlabel('Delta')
pylab.ylabel('mu')
pylab.title('mu(Delta)|phi=0.1,d0/delta=100,DeltaPM=0.01')
pylab.grid(True)
pylab.savefig('muOfDeltaPlot')
pylab.show()
'''
################################################################################

shift_mu = 0.89275*e-4 ## trial and error currently for smoothing the singularit at the pole of X(s) - this shifts the singularit away from s_star

s_array = np.arange(0.0,1.,0.001)

mu1o3Val = muThirdRoot(discriminantOfMuThirdRoot(alpha(0.01,0.01), beta(100., 0.1, 0.01) ),alpha(0.01,0.01), beta(100., 0.1, 0.01))

transformOf_dxdy_array = transformOf_dxdy(alpha(0.01,0.01),beta(100., 0.1, 0.01),mu1o3Val**3+shift_mu,s_array)

#print s_star(mu1o3Val**3)
'''
pylab.plot(s_array,transformOf_dxdy_array)
pylab.xlabel('s')
pylab.ylabel('X(s)')
pylab.title('X(s)|phi=0.1,d0/delta=100,DeltaPM=0.01,Delta=0.01')
pylab.grid(True)
pylab.savefig('X_of_s_Plot')
pylab.show()
'''


################################################################################

'''
t = numpy.arange(0.0, 1.0+0.01, 0.01)
s = numpy.cos(2*2*numpy.pi*t)
pylab.plot(t, s)

pylab.xlabel('time (s)')
pylab.ylabel('voltage (mV)')
pylab.title('About as simple as it gets, folks')
pylab.grid(True)
pylab.savefig('simple_plot')

pylab.show()
'''
################################################################################
'''
t_array = np.arange(0.05,1.,0.05)
invTrafoList = [Talbot(float(el)) for el in list(t_array)]
invTrafo_array = np.asarray( invTrafoList )
pylab.plot(t_array, invTrafo_array)
pylab.xlabel('time t')
pylab.ylabel('L^(-1) of F(s)')
pylab.title('inv. transform: exp(-t)')
pylab.grid(True)
pylab.savefig('invTrafo')

pylab.show()
'''

'''
def beta(d0_div_delta, phi, Delta):
	return 2*d0_div_delta*phi/Delta

def alpha(DeltaPM,Delta):
	return DeltaPM/Delta


'''

t_array = np.arange(0.025,10.,0.1)
ivantsov_list = [ (1./(sqrt(2.*float(el)))) for el in list(t_array)] 
ivantsov_array = np.asarray(ivantsov_list)


invTrafoList = [Talbot_for_X( float(el),24,alpha(0.01,0.01),beta(100., 0.1, 0.01),muThirdRoot(discriminantOfMuThirdRoot(alpha(0.01,0.01), beta(100., 0.1, 0.01) ),alpha(0.01,0.01), beta(100., 0.1, 0.01))**3+shift_mu) for el in list(t_array)]
invTrafo_array = np.asarray( invTrafoList )
pylab.plot(t_array, invTrafo_array,t_array,ivantsov_array)
pylab.xlabel('dist. to tip -xi')
pylab.ylabel('blue: L^(-1) of X(s), green: Ivantsov')
pylab.title('inv. transform: dx/dy vs Ivantsov solution')
pylab.grid(True)
pylab.savefig('invTrafo_dxdy_VS_IvantsovSolution')

pylab.show()























	
	
