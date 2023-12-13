# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:28:11 2020

@author: Pascal de Koster

Test the Lax pair for vNLSE (variation on Nonlinear Schroedinger equation),
which is given by
    iq_t+q_xx+2iqq^*q_x=0,
or similarly,
    q_t = iq_xx -2qq^*q_x.      (vNLSE)
Note the differenc with the standard NLSE, which is given by
    q_t = iq_xx - 2iqq^*q       (NLSE)

The Lax pair of vNLSE is apparently given by
    L = [d_x+i(q^*_x\int q dx+q^*\int q d_x dx),  
                     -i(q^*^2+q_x^*\int q^* dx -q^*\int q^*d_x dx) \\
         -i(q^2+q_x\int q dx -q\int q d_x dx),
                     -d_x-i(q_x\int q^* dx + q^*\int q^* d_x dx) ]
and 
    A = [-id_x^2-2q d_xq ^*,        -2q^* q^*_x \\
         -2q q_x           , i d_x^2 -2q^*d_x q ]
"""

from sympy import *
from Modules.PDE_initialisation import *
from Modules.LaxPairTester_mod import Lt_Psi, ALmLA, Check_Lax_equation


''' Set the PDE of interest. If the PDE is not in the list below, then
define Q_t and Q yourself.
'''


def PDE_RHS(x, t):
    #    return K_vNLSE(x,t);     #Works
    #    return K_NLSE(x,t);      #Works
    #    return Chen_Eq2(x,t);    #Works
    #    return Chen_Eq3(x,t);    #Does not work
    return K_ViscousBurgers(x, t)  # Does not work yet


#    return K_KdV(x,t);       #Works

# Custom Lax pair:
x, t = symbols('x,t')

Q_t, Q, Q_conj, Psi, foo = PDE_RHS(x, t)

'''
Test on KdV, option 1
'''
print('\nCheck the Lax pair of KdV by Peter Lax...')
# Original Lax pair by Peter Lax to validate
Q_t, Q, Q_conj, Psi, foo = K_KdV(x, t)
L_Psi = diff(Psi[0], (x, 2)) + Q[0] * Psi[0]
A_Psi = -(4 * diff(Psi[0], (x, 3)) + 6 * Q[0] * diff(Psi[0], x) + 3 * diff(Q[0], x) * Psi[0])
Lax_eq_satisfied = Check_Lax_equation(L_Psi, A_Psi, Q, Q_t, Psi, x, t)
if Lax_eq_satisfied:
    print("\nCheck_Lax_equation test 1: Passed")
else:
    print("\nCheck_Lax_equation test 1: Failed")

'''
Test on KdV, option 2
'''
print('\nCheck a Lax pair of KdV by CLL algorithm...')
# Lax pair by Chen, Lee and Liu
Q_t, Q, Q_conj, Psi, foo = K_KdV(x, t)
L_Psi = diff(Psi[0], x) + Integral(Q[0] * Psi[0], x)
A_Psi = -diff(Psi[0], (x, 3)) - 3 * Q[0] * diff(Psi[0], x)

Lax_eq_satisfied = Check_Lax_equation(L_Psi, A_Psi, Q, Q_t, Psi, x, t)

if Lax_eq_satisfied:
    print("\nCheck_Lax_equation test 2: Passed")
else:
    print("\nCheck_Lax_equation test 2: Failed")

'''
Test on KdV, option 3
'''
print('\nCheck a Lax pair of KdV by adapted CLL algorithm...')
# Lax pair 2 by Chen, Lee and Liu
Q_t, Q, Q_conj, Psi, foo = K_KdV(x, t)
L_Psi = diff(Psi[0], x, 2) + 2 * Q[0] * Psi[0] + 2 * Integral(Q[0] * diff(Psi[0], x), x)
A_Psi = -diff(Psi[0], (x, 3)) - 6 * Q[0] * diff(Psi[0], x)

Lax_eq_satisfied = Check_Lax_equation(L_Psi, A_Psi, Q, Q_t, Psi, x, t)

if Lax_eq_satisfied:
    print("\nCheck_Lax_equation test 3: Passed")
else:
    print("\nCheck_Lax_equation test 4: Failed")

'''
Test on Viscous Burgers equation
'''
print('\nCheck the Lax pair of the viscous Burgers equation...')
# Viscous Burgers equation
Q_t, Q, Q_conj, Psi, foo = K_ViscousBurgers(x, t)
alpha, beta = symbols('alpha, beta')
alpha_new = 1
Q_t = Q_t.subs(alpha, alpha_new)
alpha = alpha_new
L_Psi = diff(Psi[0], x) - alpha / (2 * beta) * Q[0] * Psi[0]
A_Psi = beta * diff(Psi[0], x, 2) - alpha * Q[0] * diff(Psi[0], x)

Lax_eq_satisfied = Check_Lax_equation(L_Psi, A_Psi, Q, Q_t, Psi, x, t)

if Lax_eq_satisfied:
    print("\nCheck_Lax_equation test 4: Passed")
else:
    print("\nCheck_Lax_equation test 4: Failed")

'''
Test on nonlinear Schroedinger equation (NLSE)
'''
print('\nCheck the Lax pair of the NLSE...')
# Viscous Burgers equation
Q_t, Q, Q_conj, Psi, foo = K_NLSE(x, t)
L_Psi = Matrix([
    [2 * Q[0] * Integral(Psi[0] * Q[1], x) - 2 * Q[0] * Integral(Psi[1] * Q[0], x) + diff(Psi[0], x)],
    [2 * Q[1] * Integral(Psi[0] * Q[1], x) - 2 * Q[1] * Integral(Psi[1] * Q[0], x) - diff(Psi[1], x)]])
A_Psi = Matrix([
    [4 * I * Psi[0] * Q[0] * Q[1] - 2 * I * Psi[1] * Q[0] ** 2 + I * diff(Psi[0], (x, 2))],
    [2 * I * Psi[0] * Q[1] ** 2 - 4 * I * Psi[1] * Q[0] * Q[1] - I * diff(Psi[1], (x, 2))]])

Lax_eq_satisfied = Check_Lax_equation(L_Psi, A_Psi, Q, Q_t, Psi, x, t)

if Lax_eq_satisfied:
    print("\nCheck_Lax_equation test 4: Passed")
else:
    print("\nCheck_Lax_equation test 4: Failed")
