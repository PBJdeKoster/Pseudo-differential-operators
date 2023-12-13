# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:27:42 2020

@author: Pascal de Koster
"""

''' This is an algorithm to do partial integration on one specific object, here
denotes as v0. The goal is to convert the integral 
    \int f(x) d^n v/dx^n dx 
to 
    sum_{i=0}^{n-1} d^i f/dx^i * d^{n-1-i} v/dx^{n-1-i} 
        + (-1)^n \int d^n f/dx^n v dx.
We note that the function should be linear in v, i.e., in each of the additve 
terms of fun, exectly one component of V appears exactly once as 
multiplication.
'''

''' 
Inputs
fun_original:   The (matrix) function to be partially integrated to V
V:          Vector containing all components of the test function to be
            integrated
x:          Spatial variable (scalar)

Ouputs:
fun_PI:     The partially integrated (PI) function which only contains 'V' in 
            its integrals, and not any of its derivatives
'''

from sympy import *


def Partially_integrate_to_test_function(fun_original, V, x):
    # Temporary variable
    q_sym = symbols('q_sym')

    fun = expand(fun_original)
    if hasattr(fun, '__len__'):
        N_row, N_col = fun.shape
    else:
        N_row = 1
        N_col = 1
        fun = Matrix([[fun]])

    # The final partially integrated (PI) function
    fun_PI = zeros(N_row, N_col)

    # loop over the elements of the (matrix) function to be rewritten
    for dim0 in range(0, N_row):
        for dim1 in range(0, N_col):
            fun_tmp = fun[dim0, dim1]
            if isinstance(fun_tmp, Add):
                fun_add = fun_tmp.args
            else:
                fun_add = [fun_tmp]

            # Initialise new additive element
            fun_tmp_add_new = 0
            # loop over additive elements of the current entry
            for i_add in range(0, len(fun_add)):
                fun_add_i = fun_add[i_add]
                if isinstance(fun_add_i, Mul):
                    fun_mul = fun_add_i.args
                else:
                    fun_mul = [fun_add_i]

                # Initialise new multiplicative element
                fun_tmp_mul_new = 1
                # loop over multiplicative elements of current additive element
                for i_mul in range(0, len(fun_mul)):
                    fun_mul_i = fun_mul[i_mul]
                    if not isinstance(fun_mul_i, Integral):
                        # Add to multiplicative terms
                        fun_tmp_mul_new = fun_tmp_mul_new * fun_mul_i
                    else:
                        int_args0 = fun_mul_i.args[0]
                        # Check if the current term contains any element of V
                        # Assume V in not in current argument
                        V_in_arg = 0
                        for i_dim_tmp in range(0, len(V)):
                            int_args = int_args0.subs(V[i_dim_tmp], q_sym)
                            if q_sym in int_args.free_symbols:
                                # We found the correct term, now replace it
                                dim_V = i_dim_tmp
                                # V is in this argument
                                V_in_arg = 1
                                break
                        if V_in_arg == 0:
                            # If the integral does not contain any element of V,
                            # then skip this entire iteration
                            fun_tmp_mul_new = fun_tmp_mul_new * fun_mul_i
                            continue

                        # We now continue with int_args and dim_V.
                        # We want to extract the derivative order of V
                        if isinstance(int_args0, Mul):
                            int_mul_args = int_args.args
                        else:
                            int_mul_args = [int_args]
                        for i2_mul in range(0, len(int_mul_args)):
                            int_mul_args_i = int_mul_args[i2_mul]
                            # Check if this is the element with V
                            if q_sym in int_mul_args_i.free_symbols:
                                if isinstance(int_mul_args_i, Derivative):
                                    # (dv/dx).args = (v(x,t) , (x,1) ),
                                    Order_der = int_mul_args_i.args[1][1]
                                else:
                                    Order_der = 0
                                break
                        # f_int is all terms in the integral apart from V
                        f_int = int_args / int_mul_args_i
                        ''' 
                        This is where the partial integration happens.
                        Let (.)_nx denote the n'th derivative to x
                        Change
                            \int f(x) v_nx dx 
                        to 
                            sum_{i=0}^{n-1} -(-1)^i f_ix * v_(n-1-i)x 
                                + (-1)^n \int f_nx v dx.
                        '''
                        # Sum of
                        f_PI0 = 0
                        for i_order in range(1, Order_der + 1):
                            f_PI0 = f_PI0 \
                                    - (-1) ** i_order * diff(f_int, (x, i_order - 1)) \
                                    * diff(V[dim_V], (x, Order_der - i_order))

                        f_PI1 = (-1) ** Order_der \
                                * Integral(diff(f_int, (x, Order_der)) \
                                           * V[dim_V], x)

                        f_PI = expand(f_PI0 + f_PI1)

                        fun_tmp_mul_new = fun_tmp_mul_new * f_PI

                # Add to additive terms
                fun_tmp_add_new = fun_tmp_add_new + fun_tmp_mul_new

            # After all additions are done, add it at the partially integrated
            # term at the correct spot
            fun_PI[dim0, dim1] = fun_tmp_add_new

    # Return a scalar function if the input was scalar
    if not hasattr(fun_original, '__len__'):
        fun_PI = fun_PI[0, 0]
    # Return the partially integrated expression
    return fun_PI;
