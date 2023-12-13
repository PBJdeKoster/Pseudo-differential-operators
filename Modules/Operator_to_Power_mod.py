# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 16:23:15 2020

@author: pdekoster1

This module contains a function which determines the N'th power of a pseudo-
differential operator
"""

from sympy import *

''' C_AB, Order_AB = Operator_Multiplication(C_A, Order_A, C_B, Order_B, x)
Take in 2 pseudo-differential (PD) operators A and B, and return the product A*B, 
which is also a pseudo-differential operator. Note that A*B != B*A.
Every operator L is represented by its coefficient vector C_L: 
    L = c2 D_x^2+ c1 D_x +c0 +c-1 D_x^-1+c_-2 D_x^-2 + ...
Only the coefficient vector C_L and Order_L is required for both operators A 
and B, where Order_L is 2 for the examplary operator L above. 
C should have the following structure:
    C = [c_0 ,c_1, c_2, ... c_Order_L, c_-N, ..., c_-2, c_-1],
such that C[0]=c_0, C[1]=c_1, C[-1]=c_-1, which is guaranteed this way. Note
that in this notation, the Order_L is required in order to know where the index
goes over from postive to negative indices. 

    Input:
C_A:            List or vector. The coefficient vector of PD-operator A.
Order_A:        Nonnegative integer. Represents the Order of L. In case the 
                order is negative, put this to 0 and fill up C_L with zeros up
                to order 0.
C_B:            List or vector. The coefficient vector of PD-operator B.
Order_B:        Nonnegative integer. Represents the Order of L. In case the 
                order is negative, put this to 0 and fill up C_L with zeros up
                to order 0.
x:              Scalar. Spatial coefficient.

    Output
C_AB:           List or vector. The coefficient vector of PD-operator A*B.
                Its order is Order_AB = Order_A + Order_B. 
Order_AB:       Nonnegative integer. The order of AB, which is just 
                Order_A + Order_B.
'''


def Operator_Multiplication(C_A, Order_A, C_B, Order_B, x):
    # The order of L^M will be Order_L*M
    Order_AB = Order_A + Order_B

    # Add extra zeros to C_L right after index Order_L, to make its size equal
    # to C_expanded:
    #   C = [c_0 ,c_1, ..., c_Order_L, c_-N, ..., c_-2, c_-1] 
    # becomes
    #   C = [c_0 ,c_1, ..., c_Order_L, 0, .........., 0, c_-N, ..., c_-2, c_-1] 
    #                               ^[0] x (M-1)*Order_L^

    # Determine the negative orders of A and B
    Order_A_n = len(C_A) - Order_A - 1
    Order_B_n = len(C_B) - Order_B - 1

    # Order_expansion is the highest order of D_x^-1 present in both A and B.
    # A further expansion is possible, but not probably not necessary.
    Order_expansion = max(Order_A_n, Order_B_n)

    ''' Apply the operator A to the operator B. 
    '''
    # The new vector will have range from D_x^Order_AB to D_x^{Order_expansion},
    # And the term for D_x^0 is also included of course, hence its size:
    C_AB = zeros(Order_AB + 1 + Order_expansion, 1)
    # Update C_AB as we loop over all operators.

    # Loop over the positive orders of D_x^i, i=0,1,2,...
    for i in range(0, Order_A + 1):
        if C_A[i] == 0:
            continue
        C_temp = zeros(len(C_AB), 1)
        # Loop over the postive right operator coefficients
        for j in range(0, Order_B + 1):
            if C_B[j] == 0:
                continue
            # Loop over the expanded effect of the D_x^i operator on the
            # coefficient
            for k in range(0, i + 1):
                C_temp[k + j] = C_temp[k + j] \
                                + binomial_coefficients_list(i)[k] * diff(C_B[j], x, i - k)
        # Loop over the negative right operator coefficients
        for j in range(-1, -Order_B_n - 1, -1):
            if C_B[j] == 0:
                continue
            # Loop over the expanded effect of the D_x^i operator on the
            # coefficient
            for k in range(0, i + 1):
                C_temp[k + j] = C_temp[k + j] \
                                + binomial_coefficients_list(i)[k] * diff(C_B[j], x, i - k)
        C_AB = C_AB + C_A[i] * C_temp

    # Loop over the strictly negative orders of D_x^i, i=-1,-2,-3, ...
    for i in range(-1, -Order_A_n - 1, -1):
        if C_A[i] == 0:
            continue
        C_temp = zeros(len(C_AB), 1)
        # Loop over the postive right operator coefficients
        for j in range(0, Order_B + 1):
            if C_B[j] == 0:
                continue
            # Loop over the expanded effect of the D_x^i operator on the
            # coefficient
            for k in range(i, -Order_expansion - j - 1, -1):
                C_temp[k + j] = C_temp[k + j] \
                                + (-1) ** (k - i) * binomial_coefficients_list(-k - 1)[-i - 1] \
                                * diff(C_B[j], x, i - k)
        # Loop over the negative right operator coefficientsc
        for j in range(-1, -Order_B_n - 1, -1):
            if C_B[j] == 0:
                continue
            # Loop over the expanded effect of the D_x^i operator on the
            # coefficient
            for k in range(i, -Order_expansion - j - 1, -1):
                C_temp[k + j] = C_temp[k + j] \
                                + (-1) ** (k - i) * binomial_coefficients_list(-k - 1)[-i - 1] \
                                * diff(C_B[j], x, i - k)
        C_AB = C_AB + C_A[i] * C_temp

    #    #Print the results for checking purposes
    #    pprint(C_AB[Order_AB:0:-1,0])
    #    pprint(C_AB[0:1,0])
    #    pprint(C_AB[-1:-Order_expansion-1:-1,0])

    return C_AB, Order_AB


''' C_LM, Order_LM = Operator_to_Power(C_L, Order_L, M, x)
Take in a Pseudo-differential operator L, and return its M'th power, which is
also a pseudo-differential operator. All these operators are represented by
their coefficient vector C: 
    L = c2 D_x^2+ c1 D_x +c0 +c-1 D_x^-1+c_-2 D_x^-2 + ...
Only the coefficient vector C_L and Order_L is required, which is 2 for the
operator above. C should have the following structure:
    C = [c_0 ,c_1, c_2, ... c_Order_L, c_-N, ..., c_-2, c_-1],
such that C[0]=c_0, C[1]=c_1, C[-1]=c_-1, which is guaranteed this way. Note
that in this notation, the Order_L is required in order to know where the index
goes over from postive to negative indices. 

    Input:
C_L:            List or vector. The coefficient vector of C_L.
Order_L:        Nonnegative integer. Represents the Order of L. In case the 
                order is negative, put this to 0 and fill up C_L with zeros up
                to order 0.
M:              Nonnegative integer. The desired power for L^M
x:              Scalar. Spatial coefficient.

    Output
C_LM:           List or vector. The coefficient vector of L^M. The order of 
                the LM is Order_L*M. 
Order_LM:       Nonnegative integer. The order of L^M, which is just Order_L*M.
'''


def Operator_to_Power(C_L, Order_L, M, x):
    """ Apply the operator L (M-1) times on itself.
    """
    # The next M-1 operations update C_new, to get
    #   L^2=>L^2=>L^3=>... ,
    # until we get to L^M. Note that the current operator L^m is always
    # represented by its C-vector, which is C_old.

    Order_new = Order_L
    C_new = C_L.copy()
    for m in range(1, M):
        C_old = C_new.copy()
        Order_old = Order_new

        C_new, Order_new = \
            Operator_Multiplication(C_L, Order_L, C_old, Order_old, x)

    C_LM = C_new.copy()
    Order_LM = Order_new

    #    #Print the results for checking purposes
    #    pprint(C_LM[Order_LM:0:-1,0])
    #    pprint(C_LM[0:1,0])
    #    pprint(C_LM[-1:-Order_expansion-1:-1,0])

    return C_LM, Order_LM
