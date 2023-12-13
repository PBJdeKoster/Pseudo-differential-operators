# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 10:34:37 2020

@author: Pascal de Koster
"""

from sympy import *

from Modules.Substitute_inside_integral_mod \
    import Substitute_inside_integral



''' Some code for testing Substitute_inside_integral
'''

# Variables, with u_t = K(u)
x, t, k, y = symbols('x t k y')
# u(x,t), the quantity of interest, for NLSE,
# this is the complex wave envelope
# u1 = u, u2 = conjugate(u)
q0 = Function('q0')(x, t)
q1 = Function('q1')(x, t)

# f = q1*diff(q0, t, x)
f = Integral(q0 * diff(q0, t, x), x)

term_old = diff(q0, t)
term_new = q1

# Ordinary subs
g1 = f.subs(term_old, term_new)
# New addition to subs
g2 = Substitute_inside_integral(g1, term_old, term_new)

print('\n   f = ')
pprint(f)
print('\nThe result of g1 = f.subs(', term_old, ',', term_new, ')')
print('\n   g1 = ')
pprint(g1)
print('\nThe result of g2 = Substitute_inside_integral(g1,',
      term_old, ',', term_new, ')')
print('\n   g2 = ')
pprint(g2)

g2_check = Integral(q0 * Derivative(q1, x), x)

if g2_check == g2.expand():
    print("\nSubstitute_inside_integral test: Passed")
else:
    print("\nSubstitute_inside_integral test: Failed")
