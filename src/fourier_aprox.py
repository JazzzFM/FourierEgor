#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 22:35:17 2022

@author: fernando
"""

from numpy import *
from numpy.fft import *


def fft_test(a):
    n = len(a)
    b = fft(a)
    print(linalg.norm(a))
    print(linalg.norm(b)/sqrt(float(n)))
    omega = complex(cos(-2.*pi/n), sin(-2.*pi/n))
    print(full((n,n), omega, dtype=complex))
    ind = outer(arange(0.,float(n)), arange(0.,float(n)))
    F = full((n,n), omega, dtype=complex)**ind
    c = F@a
    
    print(linalg.norm(c-b))

fft_test(random.rand(6))

#%%

def fourier_approx(f, m):
    n = m*m
    ind = arange(0, n)
    x = full(n, -float(m)/2.) + ind/float(m)
    s = power(full(n, -1.), ind)
    u = f(x)*s
    v = (s*fft(u))/float(m)
    
    return v

def test_fourier_approx_1(m):
    f = lambda x: exp(-pi*(x*x))
    g = f
    vapprox = fourier_approx(f, m)
    n = m*m
    xi = full(n, -float(m)/2.) + arange(0., float(n))/float(m)
    vexact = g(xi)
    err = max(abs(vexact - vapprox))
    
    print(err)
    
def test_fourier_approx_2(m):
    f = lambda x: exp(-pi*(x*x))
    g = lambda x: exp(-2.*pi*abs(x))
    vapprox = fourier_approx(f, m)
    n = m*m
    xi = full(n, -float(m)/2.) + arange(0., float(n))/float(m)
    vexact = g(xi)
    err = max(abs(vexact - vapprox))
    
    print(err)
    
test_fourier_approx_1(8)
test_fourier_approx_2(16)