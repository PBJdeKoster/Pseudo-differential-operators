# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:58:04 2020

@author: pdekoster1

Within this script, we take a differential operator as input, i.e., 
    L = \sum_{i=0}^{N} b_i(q,q_x,q_xx,...)(D_x)^i,
and return its inverse, which is a pseudo differential operator
    L_inv = \sum_{i=-N}^{-oo} c_i(q,q_x,q_xx,...)(D_x)^i.
Note that D_x^{-1} is an integration operator.

L_inv is determined through the relation 
    L*L_inv = I,
with I the identity operator. The coefficient c_i of L_inv are determined 
recursively, starting from c_{-N}.
"""

from sympy import symbols, Function, Matrix, diff, Derivative, integrate, \
    Integral, expand, sympify, pprint, solve

'''Initialise the variables'''
# Set the variable
x, t = symbols('x,t')

# Set the relevant functions: u is the potential, phi and psi are test-functions
u = Function('u')(x)
phi = Function('phi')(x)
psi = Function('psi')(x)

'''Set the DIFFERENTIAL operator to be inversed here'''
# Set the order
L_phi = Derivative(phi, x, 2) + u * phi

# Extract the order of the highest differential operator
Max_order = 7
for i in range(Max_order, -1, -1):
    # Check if L_phi contains a term of derivative order i
    if diff(L_phi, Derivative(phi, x, i)) != 0:
        Order_L = i
        break

'''Construct the inverse PSEUDO-DIFFERENTIAL operator'''
# Set u-dependent coefficients for the inverse operator
c0 = Function('c_0')(x)
cm1 = Function('c_-1')(x)
cm2 = Function('c_-2')(x)
cm3 = Function('c_-3')(x)
cm4 = Function('c_-4')(x)
cm5 = Function('c_-5')(x)
cm6 = Function('c_-6')(x)
cm7 = Function('c_-7')(x)
cm8 = Function('c_-8')(x)
cm9 = Function('c_-9')(x)

C = [c0, cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9]
C_name = C.copy()

# int_psi_xi stands for the i'th integral of psi, i.e. (D_x)^{-i}(psi)
int_psi_x0 = Function('I0psi')(x)
int_psi_x1 = Function('I1psi')(x)
int_psi_x2 = Function('I2psi')(x)
int_psi_x3 = Function('I3psi')(x)
int_psi_x4 = Function('I4psi')(x)
int_psi_x5 = Function('I5psi')(x)
int_psi_x6 = Function('I6psi')(x)
int_psi_x7 = Function('I7psi')(x)
int_psi_x8 = Function('I8psi')(x)
int_psi_x9 = Function('I9psi')(x)

int_psi_x_name = [int_psi_x0, int_psi_x1, int_psi_x2, int_psi_x3, int_psi_x4,
                  int_psi_x5, int_psi_x6, int_psi_x7, int_psi_x8, int_psi_x9]

# Set the order of the inverse operator, required for the expansion
Order_expansion = 6
Order_L_inv = Order_L + Order_expansion

int_psi = [None] * Order_L_inv
int_psi[0] = psi.copy()
# Create a template for the inverse operator
L_inv_psi = 0
for i in range(1, Order_L_inv):
    int_psi[i] = Integral(int_psi[i - 1], x)
    # Take the i'th integral

for i in range(Order_L, Order_L_inv):
    # Add the i'th integral with appropriate coefficient to the formula.
    L_inv_psi = L_inv_psi + C_name[i] * int_psi[i]

# Determine L*L_inv*psi, by substituting phi=L_inv*psi into L*phi
L_L_inv_psi = L_phi.subs(phi, L_inv_psi).expand().doit(manual=True)

# pprint(L_L_inv_psi)

# Replace all the integrals by placeholder names
for i in range(1, Order_L_inv):
    L_L_inv_psi = L_L_inv_psi.subs(int_psi[i], int_psi_x_name[i])

# pprint(L_L_inv_psi)

'''Solve for the C[i]'''
# C[i]=0 for i< Order_L
for i in range(0, Order_L):
    C[i] = 0

# Solve for C[Order_L] from the equation for order 0 derivative, this equation
# should result in I*psi = psi.
eq_0 = diff(L_L_inv_psi, psi) - 1
C[Order_L] = (solve(eq_0, C_name[Order_L])[0])
L_L_inv_psi = L_L_inv_psi.subs(C_name[Order_L], C[Order_L]).doit(manual=True)

# Solve for C[Order_L] from the equation for all other derivatives, these
# equation should all be 0, as the L*L_inv*psi=psi, so no other integrals are
# present.
for i in range(1, Order_expansion):
    eq_i = diff(L_L_inv_psi, int_psi_x_name[i])
    C[Order_L + i] = solve(eq_i, C_name[Order_L + i])[0]
    L_L_inv_psi = L_L_inv_psi.subs(C_name[Order_L + i], C[Order_L + i]) \
        .doit(manual=True)

pprint(C[0:Order_L_inv])
