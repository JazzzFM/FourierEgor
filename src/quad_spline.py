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
    n = len(min(u,v))
    m = n-1
    A = zeros((2*m, 2*m))
    b = zeros(2*m)

    for i in range(m-1):
        A[i][i] = u[i+1]-u[i]
        A[i][i+m] = (u[i+1]-u[i])**2
        b[i] = v[i+1]-v[i]
    
    A[m-1][m-1] = u[m]-u[m-1]
    A[m-1][2*m-1] = (u[m]-u[m-1])**2
    b[m-1] = v[0]-v[m]

    for i in range(m-1):
        A[i+m][i] = 1.0
        A[i+m][i+1] = -1.0
        A[i+m][i+m] = 2.0*(u[i+1]-u[i])
        
        
    A[2*m-1][m-1] = 1.0
    A[2*m-1][0] = -1.0
    A[2*m-1][2*m-1] = 2.0*(u[m]-u[m-1])
    
    print(A); print(b)
    
    x = solve(A,b)

    return x

u = list([0.0, pi/2.0, 2.0*pi])
v = list([-1.0, 4.0, -1.0])


x = coef_quad_spline(u, v)
print(x)




















                