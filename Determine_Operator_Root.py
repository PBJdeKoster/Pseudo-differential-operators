# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 11:16:18 2020

@author: pdekoster1

Within this script, we take a differential operator as input, i.e., 
    L = \sum_{i=0}^{N} b_i(q,q_x,q_xx,...)(D_x)^i,
and return its N'th root, which is a pseudo differential operator
    L_{1/M} = \sum_{i=N/M}^{-oo} c_i(q,q_x,q_xx,...)(D_x)^i.
Note that D_x^{-1} is an integration operator.

L_root is determined through the relation 
    L_root^M = L. 
The coefficients c_i of L_root are determined 
recursively, starting from c_{N/M}.

See also my document on Overleaf on some justification.
"""

from Modules.Operator_to_Power_mod import *
import sys

''' 
ONLY THE PART BELOW MAY BE TWEAKED!!!
'''
'''------------------START OF TWEAKABLE PART--------------------------------'''
'''Initialise the variables'''
# Set the variable
x, t = symbols('x,t')

# Set the relevant functions: u is the potential, psi and psi are test-functions
u = Function('u')(x)
v = Function('v')(x)
u1 = Function('u1')(x)
phi = Function('phi')(x)
psi = Function('psi')(x)

'''Set the pseudo-differential operator to be inversed here'''
# Set the differential operator, expressed as coefficients.
# C_L_p = [c_0, c_1, c_2, ..., order_L], e.g., L= u + D_x^2 => [u,0,1]
# C_L_m = [c_-1, c_-2,...]

# Set the desired power of the root, m
m_overwrite = 0
# Set the desired cutoff of L^m/M_+, where the _+ indicates that the lower part
# has been cut off. Order_cutoff indicates at which order the cutoff happens.
Order_cutoff_overwrite = 0
# Set the desired order of D_x^-{Order_expansion} up to which we want to find
# the coefficients for. Note that Order_expansion >= 0.
Order_expansion_overwrite = 0

''' Below are a couple of examples that all lead to the KdV equation. See also 
Example 2.1.26 from 'Algebaric Sturctures in Integrability'
by Sokolov, 2020.
'''
''' KdV option 1 '''
# C_L_p = [u, 0, 1]
# C_L_m = []
# m_overwrite = 1
# Order_cutoff = 0
''' KdV option 2 '''
# C_L_p = [0, 1]
# C_L_m = [u]
# m_overwrite = 3
''' KdV option 3 '''
# C_L_p = [2*u, 0, 1]
# C_L_m = [Derivative(u,x)]
''' KdV option 4 '''
# C_L_p = [0, Derivative(u,x), 2*u, 0, 1]
# C_L_m = []
''' KdV option 5 '''
##L = D_x + D_x^-1 u
# C_L_p = [0, 1]
# C_L_m = [            u     ,\
#         -Derivative(u,x,1),\
#          Derivative(u,x,2),\
#         -Derivative(u,x,3) ]
# m_overwrite = 3
# Order_expansion_overwrite = 4
''' Burgers equation '''
# C_L_p = [u, 1]
# C_L_m = []
# m_overwrite = 2
# Order_cutoff_overwrite = 1
''' Boussinesq equation '''
# C_L_p = [-(Derivative(u,x)+I*v), -2*u, 0, 1]
##C_L_p = [0, -u, 1]
# C_L_m = []
# m_overwrite = 2
# Order_cutoff_overwrite = 0
''' Harry Dym equation '''
# C_L_p = [0, -2*u*Derivative(u,x),-u**2]
##C_L_p = [0, -u, 1]
# C_L_m = []
# m_overwrite = 3
# Order_cutoff_overwrite = 1
''' Modified KdV '''
##L = 4uD_x^-1 uD_x + D_x^2
# C_L_p = [4*u**2,0,1]
# C_L_m = [-4*u*Derivative(u,x,1),\
#          4*u*Derivative(u,x,2),\
#         -4*u*Derivative(u,x,3) ]
# m_overwrite = 3
# Order_cutoff_overwrite = 0
''' Sawada-Kotera equation '''
C_L_p = [0, u, 0, 1]
C_L_m = []
m_overwrite = 5
Order_cutoff = 0
Order_expansion_overwrite = 7
''' Custom '''
# C_L_p = [1, 0, 1]
# C_L_m = []
# m_overwrite = 1
# Order_cutoff = 0
# Order_expansion_overwrite = 7

''' Set the desired order m in L^{m/M}, and Order_cutoff
'''
if m_overwrite != 0:
    m = m_overwrite
else:
    m = 1
if Order_cutoff_overwrite != 0:
    Order_cutoff = Order_cutoff_overwrite
else:
    Order_cutoff = 0

''' Set the desired root order. Choose from '1' and 'largest'
'''
Root_order = '1'
# Root_order = 'largest'

C_L_m.reverse()
C_L = Matrix(C_L_p + [0] * 9 + C_L_m)

'''Set the desired order of D_x^-{Order_expansion} up to which we want to find 
the coefficients for. Note that Order_expansion >= 0.'''
if Order_expansion_overwrite != 0:
    Order_expansion = Order_expansion_overwrite
else:
    Order_expansion = 3
'''--------------------END OF TWEAKABLE PART--------------------------------'''

# Append 9 zeros to C_L, representing the coefficients C[-1]...C[-10], wich
# are by definition all 0:
# index: [ ...,   -2,   -1,   0,   1,   2,   3, ...]
#     C: [ ..., c_-2, c_-1, c_0, c_1, c_2, c_3, ...]

# Extract the order of the differential operator L
Order_L = len(C_L_p) - 1
Order_L_m = len(C_L_m)

'''Set N in L^{1/M}. M has to be a positive integer such that it is a divisor 
of Order_L. Logical possiblities are the smallest divisor of Order_L, or 
Order_L itself.'''
if Root_order == '1':
    M = Order_L
elif Root_order == 'largest':
    if Order_L == 1:
        M = 1
    elif Order_L == 0:
        M = 0
    elif Order_L > 0:
        M = primefactors(Order_L)[0]
    else:
        sys.exit('Order_L must be an integer >= 0.')
else:
    sys.exit('Root_divisor not implemented')

# Set the order of L_root
Order_L_root = int(Order_L / M)

'''Construct the inverse PSEUDO-DIFFERENTIAL operator'''
# Set u-dependent coefficients for the inverse operator
cp9 = Function('c_9')(x)
cp8 = Function('c_8')(x)
cp7 = Function('c_7')(x)
cp6 = Function('c_6')(x)
cp5 = Function('c_5')(x)
cp4 = Function('c_4')(x)
cp3 = Function('c_3')(x)
cp2 = Function('c_2')(x)
cp1 = Function('c_1')(x)
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

# c_positive and c_negative: e.g.:
#    L^1/N = c_1 D_x^1+c_0    
#            +c_-1 D_x^-1+c_-2 D_x^-2
Cp = [cp9, cp8, cp7, cp6, cp5, cp4, cp3, cp2, cp1, c0]
Cm = [cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8, cm9]

# Concatenate Cp and Cm in such a way that C[1]=c1, C[0]=c0, C[-1]-c_-1, etc.
Cp_tmp = Cp.copy()
Cm_tmp = Cm.copy()
Cm_tmp = Cm_tmp[0:Order_expansion]
Cp_tmp.reverse()
Cm_tmp.reverse()
C_L_root = Cp_tmp[0:Order_L_root + 1] + Cm_tmp
C_name = C_L_root.copy()

''' Determine L^{M/M}, indicated by C_L_MM
'''
C_L_MM, Order_L_MM = Operator_to_Power(C_L_root, Order_L_root, M, x)

''' Solve for the coefficients C_L_MM by iteratively comparing 
L = L^{M/M} => C(L) = C(L^{M/M}) => C_L = C_L_MM
'''

# We will fill in C iteratively, starting with the highest non-zero
# coefficients, which has index Order_L. We continue downwards

Order_diff = Order_L - Order_L_root

# Note that C[Order_L_root] is solved from the equation from the equation
# for C_L_MM[Order_L], for example,
#   c_1 D_x(c_1 D_x) = c_1^2 D_x^2 + O(D_x)
# so we solve c_1 from the equation for D_x^2.
# The difference between the current unknown C and the index of the equation
# from which this C can be derived is Order_L-Order_L_root, and remains
# constant.

# Solve for positive-index coefficients of L^{1/M} (=L_root).
for i in range(Order_L_root, 0, -1):
    sol_temp = solve(C_L_MM[Order_diff + i] - C_L[Order_diff + i], C_name[i])
    sol_temp = Matrix(sol_temp).expand()
    if 1 in sol_temp:
        C_L_root[i] = 1
    else:
        C_L_root[i] = sol_temp[0]
    # Substitute the result for C_L_root[i] into C_L_MM
    for j in range(0, Order_L):
        C_L_MM[j] = C_L_MM[j].subs(C_name[i], C_L_root[i])
    for j in range(-1, -Order_expansion - 1, -1):
        C_L_MM[j] = C_L_MM[j].subs(C_name[i], C_L_root[i])

# Solve for index 0
C_L_root[0] = solve(C_L_MM[Order_diff + 0] - C_L[Order_diff + 0], C_name[0])[0] \
    .expand()
# Substitute the result for C_L_root[i] into C_L_MM
for j in range(0, Order_L):
    C_L_MM[j] = C_L_MM[j].subs(C_name[0], C_L_root[0])
for j in range(-1, -Order_expansion - 1, -1):
    C_L_MM[j] = C_L_MM[j].subs(C_name[0], C_L_root[0])

# Solve for negative indices
for i in range(-1, -Order_expansion - 1, -1):
    C_L_root[i] = solve(C_L_MM[Order_diff + i] - C_L[Order_diff + i], C_name[i])[0] \
        .expand()
    # Substitute the result for C_L_root[i] into C_L_MM
    for j in range(0, Order_L):
        C_L_MM[j] = C_L_MM[j].subs(C_name[i], C_L_root[i])
    for j in range(-1, -Order_expansion - 1, -1):
        C_L_MM[j] = C_L_MM[j].subs(C_name[i], C_L_root[i])

C_L_root = Matrix(C_L_root)
# Print the result
print('\nC_MM = \n')
pprint(C_L_root[Order_L_root:0:-1, 0])
pprint(C_L_root[0:1, 0])
pprint(C_L_root[-1:-Order_expansion - 1:-1, 0])

''' Now, determine L^{m/M}, with m mod M != 0
'''

# Compress the vector C_L_root to only include the known coefficients
C_L_root_new = zeros(Order_L_root + 1 + Order_expansion, 1)
C_L_root_new[0:Order_L_root + 1, 0] = C_L_root[0:Order_L_root + 1, 0]
if Order_expansion > 0:
    C_L_root_new[-Order_expansion:-1, 0] = C_L_root[-Order_expansion:-1, 0]
C_L_root_new[-1, 0] = C_L_root[-1, 0]

C_L_mM, Order_L_mM = Operator_to_Power(C_L_root_new, Order_L_root, m, x)

# Print the result
print('\nC_mM = \n')
pprint(C_L_mM[Order_L_mM:0:-1, 0])
pprint(C_L_mM[0:1, 0])
pprint(C_L_mM[-1:-Order_expansion - 1:-1, 0])

''' Determine the PDE corresponding to the Lax pair (L, L^{m/M})
'''

I_Psi_0 = Function('I_Psi_0')(x, t)
I_Psi_1 = Function('I_Psi_1')(x, t)
I_Psi_2 = Function('I_Psi_2')(x, t)
I_Psi_3 = Function('I_Psi_3')(x, t)
I_Psi_4 = Function('I_Psi_4')(x, t)
I_Psi_5 = Function('I_Psi_5')(x, t)
I_Psi_6 = Function('I_Psi_6')(x, t)
I_Psi_7 = Function('I_Psi_7')(x, t)
I_Psi_8 = Function('I_Psi_8')(x, t)
I_Psi_9 = Function('I_Psi_9')(x, t)

# Intialise a list of integral functions
Int_Psi = Matrix([[I_Psi_0], \
                  [I_Psi_1], \
                  [I_Psi_2], \
                  [I_Psi_3], \
                  [I_Psi_4], \
                  [I_Psi_5], \
                  [I_Psi_6], \
                  [I_Psi_7], \
                  [I_Psi_8], \
                  [I_Psi_9]])
Int_Psi_name = Int_Psi.copy()

Int_Psi[0] = psi
for i in range(0, len(Int_Psi) - 1):
    Int_Psi[i + 1, 0] = Integral(Int_Psi[i, 0], x)

# L_t_Psi, L_t_Psi_matrix = Lt_Psi(Q,Psi,Q_t)

Psi = Matrix([[psi]])

''' Old code, using A_Psi and L_Psi to determine AL and LA
L_Psi = 0
for i in range(0,Order_L+1):
    L_Psi = L_Psi + C_L[i]*Derivative(Psi[0],x,i)
for i in range(-1,-Order_L_m-1,-1):
    L_Psi = L_Psi + C_L[i]*Int_Psi[-i,0]
L_Psi = Matrix([[L_Psi]])

A_Psi = 0
for i in range(0,Order_L_mM+1):
    A_Psi = A_Psi + C_L_mM[i]*Derivative(Psi[0],x,i)
A_Psi = Matrix([[A_Psi]])

U = Matrix([[u]])

ALmLA_Psi, ALmLA_Psi_matrix =  ALmLA(L_Psi, A_Psi, U, Psi, x)

print('\nL_Psi = ')
pprint(L_Psi)
print('\nA_Psi = ')
pprint(A_Psi)
print('\n[A,L]Psi = ')
pprint(ALmLA_Psi)
'''

''' New code, using the only the PD-structure and C_L and C_A '''
Order_A = Order_L_mM
# Only take the non-negative part: A = L^{m/M}_+

C_A = C_L_mM[0:Order_A + 1, 0]
C_A[0:Order_cutoff, 0] = zeros(Order_cutoff, 1)

# Cut both operators down to the relevant size
if len(C_A) > Order_A + 1 + Order_expansion:
    C_A = list(C_A)
    C_A[Order_A + 1:-Order_expansion] = []
    C_A = Matrix(C_A)
if len(C_L) > Order_L + 1 + Order_expansion:
    C_L = list(C_L)
    C_L[Order_L + 1:-Order_expansion] = []
    C_L = Matrix(C_L)

C_AL, Order_AL = Operator_Multiplication(C_A, Order_A, C_L, Order_L, x)
C_LA, Order_LA = Operator_Multiplication(C_L, Order_L, C_A, Order_A, x)

C_ALmLA = (C_AL - C_LA).expand()

print('\n[A,L] = ')
pprint(C_ALmLA[Order_AL:0:-1, 0])
pprint(C_ALmLA[0:1, 0])
pprint(C_ALmLA[-1:-Order_expansion - 1:-1, 0])
