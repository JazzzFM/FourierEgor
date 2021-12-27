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
def coef_quad_spline(u, v):
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
    
        

u = list([0.0, pi/2.0, 3.*pi/2., 2.0*pi])
v = list([-1.0, 4.0, 4.0, -1.0])
#u = list([0., pi/2.0, 2.0*pi])
#v = list([-1., 4., -1.])

coef = coef_quad_spline(u, v)
x = linspace(0.,2.*pi,100)
#coef = array([[  0.0, -10./pi    , -40/(3.*pi*pi)],\
#              [pi/2., 10./(3.*pi), 0.0        ]])
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