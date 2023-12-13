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

from sympy import symbols, Function, diff, I, simplify, \
    Matrix, expand, integrate, pprint, zeros, Integral, sin
import sys

from Partial_Integrator_to_test_function_mod \
    import Partially_integrate_to_test_function
from Substitute_inside_integral_mod \
    import Substitute_inside_integral

sys.path.append('./../')
sys.path.append('./../Modules/')


''' L_t_Psi, L_t_Psi_matrix = Lt_Psi(Q,Psi,Q_t)
Lax L-operator for NLSE (spectral operator, Lv = lambda v)
Input
L_Psi:  The L-operator applied on test function Psi. (N x 1 matrix)
Q:      The (vector) potential function which is the PDE of interest for the 
        Lax pair L,A. (N x 1 matrix)
Psi:    The (vector) test function on which AL-LA is applied. This should
        be the same size as Q. (N x 1 matrix)
x:      The spatial coordinate (scalar)
t:      The temporal coordinate (scalar)

Ouput
L_t_Psi:    Denotes dL/dt*Psi, i.e., diff(Lv,t)-L(diff(v,t)). This is the 
        definition of the time derivative of the L-operator.
L_t_Psi_matrix:     The same as L_t_Psi, but now in unsummed matrix form
'''


def Lt_Psi(L_Psi, Q, Psi, x, t):
    # Ensure all inputs are in Matrix/list form
    if not hasattr(L_Psi, '__len__'):
        L_Psi_new = Matrix([[L_Psi]])
    else:
        L_Psi_new = Matrix(L_Psi)
    if not hasattr(Q, '__len__'):
        Q = Matrix([[Q]])
    else:
        Q = Matrix(Q)
    if not hasattr(Psi, '__len__'):
        Psi = Matrix([[Psi]])
    else:
        Psi = Matrix(Psi)

    # Initialise end result
    L_t_Psi_matrix = zeros(len(Q), len(Q))

    # Loop over the elements of L_Psi
    for i in range(0, len(Q)):
        L_Psi_i = L_Psi_new[i]
        for j in range(0, len(Q)):
            L_Psi_ij = L_Psi_i.copy()
            # Set all Psi[k!=j] to 0, such that only the elements with Psi[j]
            # are left
            for k in range(0, len(Q)):
                if k != j:
                    L_Psi_ij = L_Psi_ij.subs(Psi[k], 0)
            L_Psi_ij = L_Psi_ij.expand().doit(manual=True)

            # L_t*Psi = (L*Psi)_t - L*Psi_t
            L_Psi_ij_t = diff(L_Psi_ij, t)
            L_Psi_t_ij = L_Psi_ij.subs(Psi[j], diff(Psi[j], t))
            L_t_Psi_ij = L_Psi_ij_t - L_Psi_t_ij

            #            L_t_Psi_matrix[i,j] = L_t_Psi_ij.expand()\
            #                                            .doit(manual=True)\
            #                                            .expand()
            L_t_Psi_matrix[i, j] = L_t_Psi_ij.expand()

    L_t_Psi = zeros(len(Q), 1)
    for i in range(0, len(Q)):
        for j in range(len(Q)):
            L_t_Psi[i] = L_t_Psi[i] + L_t_Psi_matrix[i, j]

    # Return L_t_Psi in the same form as the input L_Psi
    if not hasattr(L_Psi, '__len__'):
        L_t_Psi = L_t_Psi[0, 0]

    return L_t_Psi, L_t_Psi_matrix;


''' ALmLA_Psi, ALmLA_Psi_matrix =  ALmLA(L_Psi, A_Psi, Q, Psi, x)
Return the commutation of the linear operator (matrices) A and L, 
[A,L] = AL-LA, applied on test function Psi

Input
L_Psi:  A len(Q)x1 vector, the linear operator L applied on Psi
A_Psi:  A len(Q)x1 vector, the linear operator A applied on PsiQ:      The (vector) potential function which is the PDE of interest for the 
        Lax pair L,A.
Psi:    The (vector) test function on which AL-LA is applied. This should
        be the same size as Q.
x:      Scalar, the spatial coordinate

Output
ALmLA_Psi:  A len(Q)xlen(Q) matrix, respresenting (AL-LA)Psi = [A,L]Psi, i.e., 
        the communtation of A and L applied on test function Psi.
ALmLA_Psi_matrix:   The same as ALmLA_Psi, but now in unsummed matrix form
'''


def ALmLA(L_Psi, A_Psi, Q, Psi, x):
    # Ensure all inputs are in Matrix form
    if not hasattr(L_Psi, '__len__'):
        L_Psi = Matrix([[L_Psi]])
    else:
        L_Psi = Matrix(L_Psi)
    if not hasattr(A_Psi, '__len__'):
        A_Psi = Matrix([[A_Psi]])
    else:
        A_Psi = Matrix(A_Psi)
    if not hasattr(Q, '__len__'):
        Q = Matrix([[Q]])
    else:
        Q = Matrix(Q)
    if not hasattr(Psi, '__len__'):
        Psi = Matrix([[Psi]])
    else:
        Psi = Matrix(Psi)

    ALmLA_Psi_matrix = zeros(len(Q), len(Q))

    for i in range(0, len(Q)):
        for j in range(0, len(Q)):
            # The j'th column of L and A
            L_Psi_col_j = L_Psi
            A_Psi_col_j = A_Psi
            for dim_tmp in range(0, len(Q)):
                if not dim_tmp == j:
                    L_Psi_col_j = L_Psi_col_j.subs(Psi[dim_tmp], 0)
                    A_Psi_col_j = A_Psi_col_j.subs(Psi[dim_tmp], 0)
            L_Psi_col_j = L_Psi_col_j.doit(manual=True)
            A_Psi_col_j = A_Psi_col_j.doit(manual=True)

            # The i'th row of L and A, in matrix form
            L_Psi_row_i = zeros(1, len(Q))
            A_Psi_row_i = zeros(1, len(Q))
            for dim_tmp in range(0, len(Q)):
                # Remove all terms except dim_tmp
                L_Psi_row_i[dim_tmp] = (L_Psi[i] - L_Psi[i].subs(Psi[dim_tmp], 0) \
                                        .doit(manual=True)).expand()
                A_Psi_row_i[dim_tmp] = (A_Psi[i] - A_Psi[i].subs(Psi[dim_tmp], 0) \
                                        .doit(manual=True)).expand()

            #                    | A00 |
            # LA00 = [ L00 L01 ]* | A10 | = L00*A00+L01*A10 *Psi0

            # First AL
            for dim_tmp in range(0, len(Q)):
                ALmLA_Psi_matrix[i, j] = ALmLA_Psi_matrix[i, j] + \
                                         A_Psi_row_i[dim_tmp].subs(Psi[dim_tmp],
                                                                   L_Psi_col_j[dim_tmp])
            # Then LA
            for dim_tmp in range(0, len(Q)):
                ALmLA_Psi_matrix[i, j] = ALmLA_Psi_matrix[i, j] - \
                                         L_Psi_row_i[dim_tmp].subs(Psi[dim_tmp],
                                                                   A_Psi_col_j[dim_tmp])

            # Expand term by term
            for dim_tmp in range(0, len(Q)):
                ALmLA_Psi_matrix[i, j] = ALmLA_Psi_matrix[i, j] \
                    .expand().doit(manual=True)
    ALmLA_Psi_matrix[i, j] = ALmLA_Psi_matrix[i, j].expand()

    ALmLA_Psi = zeros(len(Q), 1)
    for i in range(0, len(Q)):
        for j in range(len(Q)):
            ALmLA_Psi[i] = ALmLA_Psi[i] + ALmLA_Psi_matrix[i, j]

    return ALmLA_Psi, ALmLA_Psi_matrix;


''' ALmLA(Q,Psi, A_Psi, L_Psi)
Return the commutation of the linear operator (matrices) A and L, 
[A,L] = AL-LA, applied on test function Psi

Input
L_Psi:  A len(Q)x1 vector, the linear operator L applied on Psi
A_Psi:  A len(Q)x1 vector, the linear operator A applied on Psi
Q:      The (vector) potential function which is the PDE of interest for the 
        Lax pair L,A.
Psi:    The (vector) test function on which AL-LA is applied. This should
        be the same size as Q.

Output
True/False:     Depending on whether the Lax equation is satisfied (True) or 
        not (False). (Boolean)
'''


def Check_Lax_equation(L_Psi, A_Psi, Q, Q_t, Psi, x, t):
    # Ensure all inputs are in Matrix form
    if not hasattr(L_Psi, '__len__'):
        L_Psi = Matrix([[L_Psi]])
    else:
        L_Psi = Matrix(L_Psi)
    if not hasattr(A_Psi, '__len__'):
        A_Psi = Matrix([[A_Psi]])
    else:
        A_Psi = Matrix(A_Psi)
    if not hasattr(Q, '__len__'):
        Q = Matrix([[Q]])
    else:
        Q = Matrix(Q)
    if not hasattr(Q_t, '__len__'):
        Q_t = Matrix([[Q_t]])
    else:
        Q_t = Matrix(Q_t)
    if not hasattr(Psi, '__len__'):
        Psi = Matrix([[Psi]])
    else:
        Psi = Matrix(Psi)

    foo, L_t_Psi_matrix, = Lt_Psi(L_Psi, Q, Psi, x, t)
    foo, ALv_m_LAv_matrix = ALmLA(L_Psi, A_Psi, Q, Psi, x)
    Lt_m_AcomL = zeros(len(Q), len(Q))

    # Take L_t_Psi_matrix, and substitute diff(Q,t) by Q_t for each element in
    # Q.
    for i in range(0, len(Q)):
        for j in range(0, len(Q)):
            L_t_Psi_ij = L_t_Psi_matrix[i, j]
            # Replace each element diff(Q,t) with its RHS Q_t.
            for dim_tmp in range(0, len(Q)):
                # Normal substitution
                L_t_Psi_ij = L_t_Psi_ij.subs(diff(Q[dim_tmp], t), Q_t[dim_tmp])
                # Substitution inside integrals (special function required)
                L_t_Psi_ij = Substitute_inside_integral( \
                    L_t_Psi_ij, diff(Q[dim_tmp], t), Q_t[dim_tmp])
                # Substitute the element without diff(Q,t) into L_t_Psi
                L_t_Psi_matrix[i, j] = L_t_Psi_ij

    # Check if L_t = AL-LA for each element.
    for i_el in range(0, len(ALv_m_LAv_matrix)):
        # L_t - [A,L]  (say: L t minus A commute L)
        Lt_m_AcomL_tmp = (L_t_Psi_matrix[i_el] - ALv_m_LAv_matrix[i_el]) \
            .expand()

        Lt_m_AcomL_tmp = Partially_integrate_to_test_function( \
            Lt_m_AcomL_tmp, Psi, x \
            ) \
            .expand().doit(manual=True)

        Lt_m_AcomL[i_el] = Lt_m_AcomL_tmp.expand()

    #    print('\nLt_m_AcomL[0,0] = ')
    #    print(Lt_m_AcomL[0,0])

    if Lt_m_AcomL == zeros(len(Q), len(Q)):
        print('Lax equation is satisfied')
        return True;
    else:
        print('Lax equation is NOT satisfied')
        return False;
