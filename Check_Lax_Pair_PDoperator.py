# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 16:11:52 2020

@author: pdekoster1
"""

from Modules.Operator_to_Power_mod import *
from Modules.PDE_initialisation import *

'''Initialise the variables'''
# Set the variable
x, t = symbols('x,t')

# Set the relevant functions: u is the potential, psi and psi are test-functions
u = Function('u')(x, t)
v = Function('v')(x, t)
u1 = Function('u1')(x)
phi = Function('phi')(x)
psi = Function('psi')(x)

''' Take a pseudo-differential operator (PD-op) H of the form
    H = f_N*D_x^N+....+f_2*D_x^2 +f_1*D_x +f_0 +f_-1*D_x^{-1}+...
and represent this as vector C_H, with elements
    C = [ C[0], C[1], ..., C[N], C[-P], C[-P+1], ..., C[-2], C[-1] ], 
with 
    C = [ f_0,  f_1, ...,  f_N,  f_-P,  f_-P+1, ...,  f_-2,  f_-1 ]
such that calling C[-1] (i.e., the last elements) will indeed give us f_-1 in
python.
'''

a, b, c, d = symbols('a,b,c,d')

''' Harry Dym equation, u_t = -2uu_xxx '''
# C_L_p = [0, -2*u*Derivative(u,x), -u**2]
# C_L_m = []
# C_A_p = [                                               0,\
#         -c*u*(2*Derivative(u,x)**2+d*u*Derivative(u,x,2)),\
#                                 -b*u**2*Derivative(u,x),\
#                                                  -a*u**2 ]
# C_A_m = []
''' Sawada-Kotera equation, u_t = -u_xxxxx-5u^2u_x-5uu_xxx-5u_xu_xx '''
Q_t, Q, Q_conj, Psi, foo = K_SawadaKotera(x, t)

# C_L_p = [0, Q[0], 0, 1]
# C_L_m = []
# C_A_p = [0, 5*Q[0]**2+10*Derivative(Q[0],x,2), 15*Derivative(Q[0],x), \
#         15*Q[0], 0, 9]
# C_A_m = []

C_L_p = [3 * Q[0], Derivative(Q[0], x), 1 - Q[0]]
C_L_m = [-Derivative(Q[0], x, 1), \
         +Derivative(Q[0], x, 2), \
         -Derivative(Q[0], x, 3), \
         +Derivative(Q[0], x, 4), \
         -Derivative(Q[0], x, 5)]
C_A_p = [0, -5 * Q[0] ** 2 - 10 * Derivative(Q[0], x, 2), -10 * Derivative(Q[0], x), \
         -5 * Derivative(Q[0], x), 0, -1]
C_A_m = [0, 0, 0, 0, 0]

# Construct the C_L and C_A from the positive and negative parts,
# and extract their orders
C_L_m.reverse()
C_L = Matrix(C_L_p + C_L_m)
Order_L = len(C_L_p) - 1
Order_L_m = len(C_L_m)

C_A_m.reverse()
C_A = Matrix(C_A_p + C_A_m)
Order_A = len(C_A_p) - 1
Order_A_m = len(C_A_m)

''' Check Lax equation '''
C_AL, Order_AL = Operator_Multiplication(C_A, Order_A, C_L, Order_L, x)
C_LA, Order_LA = Operator_Multiplication(C_L, Order_L, C_A, Order_A, x)

C_ALmLA = (C_AL - C_LA).expand()

C_Lt = diff(C_L, t).expand()
C_Lt_new = C_Lt.subs(Derivative(Q[0], t), Q_t[0]).doit(manual=True).expand()

C_Lt_m_AcomL = -C_ALmLA.copy()
for i in range(0, Order_AL):
    if Order_L >= i:
        C_Lt_m_AcomL[i] = C_Lt_m_AcomL[i] + C_Lt_new[i]
for i in range(-1, -len(C_AL) + Order_AL, -1):
    if -i < len(C_L) - Order_L:
        C_Lt_m_AcomL[i] = C_Lt_m_AcomL[i] + C_Lt_new[i]

print('\nL_t = ')
pprint(C_Lt_new[Order_AL:0:-1, 0])
pprint(C_Lt_new[0:1, 0])

print('\n[A,L] = ')
pprint(C_ALmLA[Order_AL:0:-1, 0])
pprint(C_ALmLA[0:1, 0])

print('\nL_t-[A,L] = ')
pprint(C_Lt_m_AcomL[Order_AL:0:-1, 0])
pprint(C_Lt_m_AcomL[0:1, 0])
