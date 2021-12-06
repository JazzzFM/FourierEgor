#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 23:52:13 2021

@author: fernando
"""

# -This program attemps to calculate the coefficients of a cuadratic periodic 
#  spline given a set of points.

from numpy import *
from numpy.linalg import *

# -coef_quad_spline returns the coefficients of our spline.
# -'u' and 'v' are arrays withc together represents the points nedeed.
def coef_quad_spline(u, v):
    n = min(len(u),len(v))       # -In case of having len(u) != len(v).
    m = n-1                      # -Number of coeficients to return.
    A = zeros((2*m, 2*m))        # -Initializing our linear sistem with
    b = zeros(2*m)               #  zeros.
    
    # -Creating our linear sistem.
    # -Using the b_i(u_{i+1} - u_i) + c_i(u_{i+1} - u_i)^2 = v_{i+1} - v_i
    #  equations.
    for i in range(m-1):
        # -The A[i][i] numbers represent the b_i's coefficients.
        A[i][i] = u[i+1]-u[i]
        # -The A[i][i+m] numbers represent the c_i's coefficients.
        A[i][i+m] = (u[i+1]-u[i])**2
        b[i] = v[i+1]-v[i]
    
    # -Introducing the b_m(u_{m+1}-u_m) + c_m(u_{m+1}-u_m)^2 = v_1-v_m
    #  equation.
    A[m-1][m-1] = u[m]-u[m-1]
    A[m-1][2*m-1] = (u[m]-u[m-1])**2
    b[m-1] = v[0]-v[m]

    # -Still creating our linear sistem.
    # -Using the b_i - b_{i+1} + 2c_i(u_{i+1}-u_i) = 0 equations.
    for i in range(m-1):
        A[i+m][i] = 1.0
        A[i+m][i+1] = -1.0
        A[i+m][i+m] = 2.0*(u[i+1]-u[i])
        
    # -Introducing the b_m -b_1 +2c_m(u_{m+1}-u_m) = 0 equation.    
    A[2*m-1][m-1] = 1.0
    A[2*m-1][0] = -1.0
    A[2*m-1][2*m-1] = 2.0*(u[m]-u[m-1])
    
    #Debugging porpuses.
    print(A); print(b)
    
    # -Solving our linear sistem.
    x = solve(A,b)

    return x

u = list([0.0, pi/2.0, 2.0*pi])
v = list([-1.0, 4.0, -1.0])


x = coef_quad_spline(u, v)
print(x)
         