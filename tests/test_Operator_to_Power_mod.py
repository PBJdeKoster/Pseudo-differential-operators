# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 16:23:15 2020

@author: pdekoster1

Testing the module Operator_to_Power_mod
"""

from sympy import *

from Modules.Operator_to_Power_mod import Operator_to_Power, Operator_Multiplication

x = symbols('x')
u = Function('u')(x)

''' Operator_Multiplication(C_A, Order_A, C_B, Order_B, x)
'''

# L = D_x^-3
C_L = Matrix([[0], [1], [0], [0]])
Order_L = 0
# A = u
C_A = Matrix([[u], [0], [0], [0], [0], [0], [0], [0]])
Order_A = 0

# AL = uD_x^-2
AL_check = Matrix([[0], [0], [0], [0], [0], [u], [0], [0]])
# LA is more complicated...
LA_check = Matrix([
    [0],
    [15.0 * Derivative(u, (x, 4))],
    [-10.0 * Derivative(u, (x, 3))],
    [6.0 * Derivative(u, (x, 2))],
    [-3.0 * Derivative(u, x)],
    [u],
    [0],
    [0]])

Order_AL_check = 0
Order_LA_check = 0

AL, Order_AL = Operator_Multiplication(C_A, Order_A, C_L, Order_L, x)
LA, Order_LA = Operator_Multiplication(C_L, Order_L, C_A, Order_A, x)

failed = 0
if Order_AL_check != Order_AL:
    print('Failed in Order_AL')
    failed = 1
if Order_LA_check != Order_LA:
    print('Failed in Order_LA')
    failed = 1
if AL != AL:
    print('Failed in AL')
    failed = 1
if AL != AL:
    print('Failed in LA')
    failed = 1

if failed == 0:
    print('Operator_Multiplication, Test 1: Passed')
else:
    print('Operator_Multiplication, Test 1: Failed')

''' Operator_to_Power(C_L, Order_L, M, x), Test 1
'''

# L = D_x + uD_x
C_L = Matrix([u, 1])
Order_L = 1
Power = 2

C_LM_check = Matrix([[Derivative(u, x) + u ** 2], [2 * u], [1]])

Order_LM_check = Order_L * Power

C_LM, Order_LM = Operator_to_Power(C_L, Order_L, Power, x)

failed = 0
if Order_LM_check != Order_LM:
    print('Failed in Order_LM')
    failed = 1
if (C_LM_check - C_LM).expand() != zeros(len(C_LM), 1):
    print('Failed in L^M')
    failed = 1

if failed == 0:
    print('Operator_to_Power, Test 2: Passed')
else:
    print('Operator_to_Power, Test 2: Failed')

''' Operator_to_Power(C_L, Order_L, M, x), Test 2
'''

# L = D_x + uD_x
C_L = Matrix([u, 0, 1])
Order_L = 0
Power = 2

C_LM_check = Matrix([[u ** 2], [1 - Derivative(u, x)], [2 * u]])

Order_LM_check = Order_L * Power

C_LM, Order_LM = Operator_to_Power(C_L, Order_L, Power, x)

failed = 0
if Order_LM_check != Order_LM:
    print('Failed in Order_LM')
    failed = 1
if (C_LM_check - C_LM).expand() != zeros(len(C_LM), 1):
    print('Failed in L^M')
    failed = 1

if failed == 0:
    print('Operator_to_Power, Test 3: Passed')
else:
    print('Operator_to_Power, Test 3: Failed')

''' Operator_to_Power(C_L, Order_L, M, x), Test 3
'''

# L = D_x + uD_x
C_L = Matrix([0, 1, 0, 0, 0, u])
Order_L = 1
Power = 3

C_LM_check = \
    Matrix([
        [3 * Derivative(u, x)],
        [3 * u],
        [0],
        [1],
        [0],
        [3 * u ** 2 + Derivative(u, (x, 2))]])

Order_LM_check = Order_L * Power

C_LM_temp, Order_LM = Operator_to_Power(C_L, Order_L, Power, x)

Order_LM_n = len(C_LM_check) - Order_LM - 1
C_LM = zeros(len(C_LM_check), 1)
C_LM[0:Order_LM + 1, 0] = C_LM_temp[0:Order_LM + 1, 0]
C_LM[-1, 0] = C_LM_temp[-1, 0]
C_LM[-Order_LM_n:-1, 0] = C_LM_temp[-Order_LM_n:-1, 0]

failed = 0
if Order_LM_check != Order_LM:
    print('Failed in Order_LM')
    failed = 1
if (C_LM_check - C_LM).expand() != zeros(len(C_LM), 1):
    print('Failed in L^M')
    failed = 1

if failed == 0:
    print('Operator_to_Power, Test 4: Passed')
else:
    print('Operator_to_Power, Test 4: Failed')
