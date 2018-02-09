#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 10:33:30 2018

Independent confirmation of the C++ and R routines.

@author: bayerc
"""

from rBergomiBen import *

xi = 0.07
H = 0.07
eta = 2.2
rho = -0.9
T = 1
K = 100
S0 = 100
N = 200
M = 50000

chol = CholeskyPricer(H, N, T, rho, xi, S0, eta)

chol.get_option_price("c", K, M)