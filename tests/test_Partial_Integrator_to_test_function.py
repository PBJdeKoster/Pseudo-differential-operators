# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:27:42 2020

@author: Pascal de Koster
"""

''' This is an algorithm to do partial integration on one specific object, here
denoted as v0. The goal is to convert the integral 
    \int f(x) d^n v/dx^n dx 
to 
    sum_{i=0}^{n-1} d^i f/dx^i * d^{n-1-i} v/dx^{n-1-i} 
        + (-1)^n \int d^n f/dx^n v dx.
We note that the function should be linear in v, i.e., in each of the additve 
terms of fun, exectly one component of V appears exactly once as 
multiplication.
'''

from sympy import *
from Modules.Partial_Integrator_to_test_function_mod \
    import Partially_integrate_to_test_function

# Variables, with u_t = K(u)
x, t, k, y = symbols('x t k y')
# q(x,t), the quantity of interest, for vNLSE,
# this is the complex wave envelope
# q0 = q, q1 = conjugate(q)
q0 = Function('q0')(x, t)
q1 = Function('q1')(x, t)
Q = Matrix([q0, q1])

# v(x,t), the eigenfunction/auxiliary function
v0 = Function('v0')(x, t)
v1 = Function('v1')(x, t)
V = Matrix([v0, v1])

q_sym = symbols('q_sym')

# fun = Matrix([[Integral(q0*diff(q1,x,2)*diff(v0,x), x)]])
fun = Integral(q0 * diff(q1, x, 2) * diff(v0, x), x)

fun_PI_check = Q[0] * V[0] * diff(Q[1], (x, 2)) \
               - Integral(Q[0] * V[0] * diff(Q[1], (x, 3)), x) \
               - Integral(V[0] * diff(Q[0], x) * diff(Q[1], (x, 2)), x)
fun_PI = Partially_integrate_to_test_function(fun, V, x)

print('\nOriginal function = ')
pprint(fun)
print('\nFunction partially integrated to v = ')
pprint(fun_PI)

if fun_PI_check == fun_PI.expand():
    print("\nPartially_integrate_to_test_function test: Passed")
else:
    print("\nPartially_integrate_to_test_function test: Failed")
