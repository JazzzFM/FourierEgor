#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 23:52:13 2021

@author: fernando
"""

# -This program attemps to calculate the coefficients of a cuadratic periodic 
#  spline given a set of points.

from pylab import *
from numpy.linalg import *

# -coef_quad_spline returns the coefficients of our spline.
# -'u' and 'v' are arrays withc together represents the points nedeed.
def coef_quad_spline2(u, v):
    n = min(len(u),len(v))       # -In case of having len(u) != len(v).
    m = n-1                      # -Number of coeficients to return.
    A = zeros((2*m, 2*m), float)        # -Initializing our linear sistem with
    b = zeros(2*m, float)               #  zeros.
    
    # -Creating our linear sistem.
    # -Using the b_i(u_{i+1} - u_i) + c_i(u_{i+1} - u_i)^2 = v_{i+1} - v_i
    #  equations.
    for i in range(m-1):
        # -The A[i][i] numbers represent the b_i's coefficients.
        A[i][i] = u[i+1]-u[i]
        # -The A[i][i+m] numbers represent the c_i's coefficients.
        A[i][i+m] = (u[i+1]-u[i])*(u[i+1]-u[i])
        b[i] = v[i+1]-v[i]
    
    # -Introducing the b_m(u_{m+1}-u_m) + c_m(u_{m+1}-u_m)^2 = v_1-v_m
    #  equation.
    A[m-1][m-1] = u[m]-u[m-1]
    A[m-1][2*m-1] = (u[m]-u[m-1])*(u[m]-u[m-1])
    b[m-1] = v[m]-v[m-1]

    # -Still creating our linear sistem.
    # -Using the b_i - b_{i+1} + 2c_i(u_{i+1}-u_i) = 0 equations.
    for i in range(m-1):
        A[i+m][i] = 1.0
        A[i+m][i+1] = -1.0
        A[i+m][i+m] = 2.0*(u[i+1]-u[i])
        
    # -Introducing the b_m -b_1 +2c_m(u_{m+1}-u_m) = 0 equation.    
    A[2*m-1][m-1] = 1.0
    A[2*m-1][2*m-1] = 2.0*(u[m]-u[m-1])
    
    # -Playing with the derivatives.
    A[2*m-1][0] = -1.0  
    #A[2*m-1][0] = 1.0
    if det(A) == 0:
        A[2*m-1][0] = 1.0
    
    #Debugging porpuses.
    #print(A); print(b)
    
    # -Solving our linear sistem.
    x = solve(A,b)
    
    # -This is just to show our results properly.
    y = zeros((m,3), float)
    for i in range(m):
        y[i][0] = v[i]
        y[i][1] = x[i]
        y[i][2] = x[i+m]

    return y


def sumalt(z, r, star, stop):
    result = 0
    for i in range(star, stop+1):
        result += (-1)**(i+r)*z[i]
        
    return result


def coef_quad_spline(u, v):
    m = len(u)-1
    w = zeros(m,float)
    y = zeros(m,float)
    
    for i in range(m):
        w[i] = u[i+1]-u[i]
        y[i] = v[i+1]-v[i]
    
    z = y/w
    coef = zeros((m,3), float)
    const = sumalt(z, 0, 0, m-2)
    for i in range(m-1):
        coef[i,0] = v[i]
        coef[i,1] = (-1)**i*z[m-1] + 2.*sumalt(z, -i, i, m-2)\
            + (-1)**(i+1)*const
    
    coef[m-1,0] = v[m-1]
    coef[m-1,1] = z[m-1]-const
    
    for i in range(m-2):
        coef[i,2] = ((-1)**(i+1)*z[m-1]-z[i] + 2.*sumalt(z,-i-1,i+1,m-2)\
                     + (-1)**(i)*const)/w[i]

    coef[m-2,2] = (z[m-1]-z[m-2] - const)/w[m-2]
    coef[m-1,2] = const/w[m-1]

    return coef        


# -eval_quad_spline function evaluates a quadratic spline over a vector x.
# -u isa vector with the partition of the domain.
# -coef is a nx3 matrix witch contains the coeficients of our spline.
def eval_quad_spline(u, coef, x):
    n = len(x)
    m = len(u)
    y = zeros(n)
    
    for i in range(n):
        
        # -In case of having elements outside the domain.
        if x[i] < u[0] or x[i] > u[m-1]: 
            continue
        
        # -Finding the [u_{j-1}, u_j) interval where x[i] belongs.
        for j in range(1,m):
            if u[j] > x[i]:
                z = x[i]-u[j-1]
                y[i] = coef[j-1][0]+coef[j-1][1]*z+coef[j-1][2]*z**2
                break
        if u[m-1] == x[i]:
            z = x[i]-u[m-2]
            y[i] = coef[m-2][0]+coef[m-2][1]*z+coef[m-2][2]*z**2
            
    return y 
    
        

#u = array([0.0, pi/2.0, 3.*pi/2., 2.0*pi])
#v = array([-1.0, 4.0, 4.0, -1.0])
#u = array([0., pi/3.0, 2.0*pi/3., 4.*pi/3., 5.*pi/3., 2*pi])
#v = array([-1., 4., 6, 6, 4, -1.])
u = array([0., pi/4., pi/2., 3.*pi/4., 5.*pi/4., 3.*pi/2., 7.*pi/4., 2.*pi])
v = array([-1., 4., 6., 8., 8., 6., 4., -1.])
#u = array([0., pi/5., 2.*pi/5., 3.*pi/5., 4.*pi/5., 6.*pi/5., 7.*pi/5.,\
#          8.*pi/5., 9.*pi/5, 2*pi])
#v = array([-1., 4., 6., 8.,10,10, 8., 6., 4., -1.])
"""
coef = coef_quad_spline(u, v)
#coef2 = coef_quad_spline2(u, v)

#print(coef)
#print(coef2)

#%%

x = linspace(0., 2.*pi, 100)
y = eval_quad_spline(u, coef, x)

for i in range(len(u)):
    print('f(',u[i],') =',eval_quad_spline(u, coef, list([u[i]]) )[0] )

#print('f(', u[0],') =', y[0])
#print('f(', u[len(u)-1], ') = ', y[len(y) -1])

plot(u,v,'.',lw=2,color='r')
plot(x,y,'-',lw=1,color='b')
grid(True)
show()
"""

#%%

def single_fourier_coef(u,coef,k):
    m = len(u)
    a = squeeze(coef[:,0])
    b = squeeze(coef[:,1])
    c = squeeze(coef[:,2])
    
    u1 = u[0:m-1]; u2 = u[1:m]
    w = u2-u1
    w2 = w*w
    w3 = w2*w
    if k==0:
        Reh_k = sum(a*w + b*w2/2. + c*w3/3.)
        Reh_k /= 2.*pi
        return complex(Reh_k, 0.)
    
    sinu1u2 = sin(k*u2) - sin(k*u1)
    cosu1u2 = cos(k*u1) - cos(k*u2)
    
    Reh_k = sum(a*sinu1u2 + b*(w*sin(k*u2) - cosu1u2/k)\
                + c*(w2*sin(k*u2) + 2.*w*cos(k*u2)/k - 2.*sinu1u2/(k*k)))
    Reh_k /= 2.*pi*k
    
    Imh_k = sum(-a*cosu1u2 + b*(w*cos(k*u2)-sinu1u2/k)\
                + c*(w2*cos(k*u2) - 2.*w*sin(k*u2)/k + 2.*cosu1u2/(4.*k*k)))
    Imh_k /= 2.*pi*k
    
    return complex(Reh_k, Imh_k)


def fourier_coef_cuadr_spline(u,v,n):
    al = zeros(n, dtype=complex)
    coef = coef_quad_spline(u, v)
    
    for k in range(n):
        al[k] = single_fourier_coef(u, coef, k)
    
    return al
    
    
#FC = fourier_coef_cuadr_spline(u, v, 4)
#print("Coeficientes de fourier:")
#print(FC)    
    
    
def eval_trig_sum(al, x):
    y = al[0]*ones(len(x),dtype=complex)
    
    n = len(al)
    for k in range(1,n):
        coefcos = complex(2.*real(al[k]), 0.)
        coefsin = complex(0., 2.*imag(al[k]))
        y = y+coefcos*cos(k*x) + coefsin*sin(k*x)
        
    return y

def test_quad_spline_fourier_sum(u,v,n):
    coef = coef_quad_spline(u, v)
    xg = linspace(0., 2.*pi, 601)
    yg = eval_quad_spline(u, coef, xg)
    al = fourier_coef_cuadr_spline(u, v, n)
    tg = eval_trig_sum(al, xg)
    err = max(abs(tg-yg))
    print(err)
    
test_quad_spline_fourier_sum(u, v, 200)