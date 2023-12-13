# Pseudo-differential-operators (PD operator)
Generating and manipulating Pseudo-differential operators relevant for Lax pairs for
analysis of Lax-intergrable systems

# Author information
Pascal de Koster, 2021
Formerly with TU Delft
Contact: pbj.dekoster@gmail.com

# Summary
This code base is designed to determine test Pseudo-differential operators, and find the inverse.
For example, the following is a Pseudo-differential operator:
	L = D_{xx} + u
This Pseudo-differential operator actually belongs to the Korteweg-de Vries equation.

This code base can do 3 things:

# Check Lax Pair of a PD operator
See also Check_Lax_Pair_PDoperator.py

Take a pseudo-differential operator (PD-op) H of the form
    H = f_N*D_x^N+....+f_2*D_x^2 +f_1*D_x +f_0 +f_-1*D_x^{-1}+...
and represent this as vector C_H, with elements
    C = [ C[0], C[1], ..., C[N], C[-P], C[-P+1], ..., C[-2], C[-1] ], 
with 
    C = [ f_0,  f_1, ...,  f_N,  f_-P,  f_-P+1, ...,  f_-2,  f_-1 ]
such that calling C[-1] (i.e., the last elements) will indeed give us f_-1 in
python.


# Determine the inverse of a PD operator
See also Determine_Operator_Inverse.py

Within this script, we take a differential operator as input, i.e., 
    L = \sum_{i=0}^{N} b_i(q,q_x,q_xx,...)(D_x)^i,
and return its inverse, which is a pseudo differential operator
    L_inv = \sum_{i=-N}^{-oo} c_i(q,q_x,q_xx,...)(D_x)^i.
Note that D_x^{-1} is an integration operator.

L_inv is determined through the relation 
    L*L_inv = I,
with I the identity operator. The coefficient c_i of L_inv are determined 
recursively, starting from c_{-N}.

# Determine the root of a PD operator
See also Determine_Operator_Root.py

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