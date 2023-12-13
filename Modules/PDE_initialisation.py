# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 17:42:27 2020

@author: Pascal de Koster

This file contains several Partial Differential Equation (PDEs), which are
initialised with their required attributes for the CLL algorithm.
The initialisation for each PDE takes in a spatial and temporal variable, and
sets a range of attributes appopriately
Inputs:
x:          Scalar, spatial variable
t:          Scalar, temporal variable

Outputs:    
Q_t:        Nx1 matrix, containing the RHS of the evolution equation.
Q:          Nx1 matrix, containing the undefined Functions dependent of x and t
Q_conj:     Nx1 matrix, containing the conjugates of Q expressed in terms of Q,
            e.q., if Q=[[q0],[q1]]=[[q0],[q0*]], (so q1= q0*, and q0*=q1)
            Q_conj = [[q0*],[q1*]]=[[q1],[q0]].

"""

from sympy import *


# --------------------------------
# The PDE of interest
# Right hand side of the evolution equation for NLSE
def K_vNLSE(x, t):
    q0 = Function('q0')(x, t)
    q1 = Function('q1')(x, t)
    Q = Matrix([[q0], [q1]])
    Q_conj = Matrix([[q1], [q0]])
    # in u_t = K(u)
    # vNLSE: q_t = iq_xx-2(q*)qq_x
    Q_t0 = I * diff(Q[0], x, 2) - 2 * Q[0] * Q[1] * diff(Q[0], x)
    Q_t1 = -I * diff(Q[1], x, 2) - 2 * Q[1] * Q[0] * diff(Q[1], x)
    Q_t = Matrix([[Q_t0], [Q_t1]])

    psi0 = Function('psi0')(x, t)
    psi1 = Function('psi1')(x, t)
    Psi = Matrix([[psi0], [psi1]])
    L_l_Psi = Matrix([[diff(Psi[0], x), 0], \
                      [0, -diff(Psi[1], x)]])

    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def K_NLSE(x, t):
    # kappa=+1 for focussing, kappa=-1 for defocussing
    kappa = 1
    q0 = Function('q0')(x, t)
    q1 = Function('q1')(x, t)
    Q = Matrix([[q0], [q1]])
    Q_conj = Matrix([[q1], [q0]])
    # in u_t = K(u)
    # NLSE: q_t = iq_xx+2i(qq*)q
    Q_t0 = I * diff(Q[0], x, 2) + 2 * kappa * I * Q[0] * Q[1] * Q[0]
    Q_t1 = -I * diff(Q[1], x, 2) - 2 * kappa * I * Q[1] * Q[0] * Q[1]
    Q_t = Matrix([[Q_t0], [Q_t1]])

    psi0 = Function('psi0')(x, t)
    psi1 = Function('psi1')(x, t)
    Psi = Matrix([[psi0], [psi1]])
    L_l_Psi = Matrix([[diff(Psi[0], x), 0], \
                      [0, -diff(Psi[1], x)]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi,;


def K_NLSE_defocussing(x, t):
    kappa = -1
    q0 = Function('q0')(x, t)
    q1 = Function('q1')(x, t)
    Q = Matrix([[q0], [q1]])
    Q_conj = Matrix([[q1], [q0]])
    # in u_t = K(u)
    # NLSE: q_t = iq_xx+2i(qq*)q
    Q_t0 = I * diff(Q[0], x, 2) + 2 * kappa * I * Q[0] * Q[1] * Q[0]
    Q_t1 = -I * diff(Q[1], x, 2) - 2 * kappa * I * Q[1] * Q[0] * Q[1]
    Q_t = Matrix([[Q_t0], [Q_t1]])

    psi0 = Function('psi0')(x, t)
    psi1 = Function('psi1')(x, t)
    Psi = Matrix([[psi0], [psi1]])
    L_l_Psi = Matrix([[diff(Psi[0], x), 0], \
                      [0, -diff(Psi[1], x)]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi,;


def K_Manakov(x, t):
    # Focussing: kappa =1, Defocussing: kappa = -1
    kappa = 1
    # Manakov equation, Q = [q0, q0*, q1, q1*]
    q0 = Function('q0')(x, t)
    q1 = Function('q1')(x, t)
    q2 = Function('q2')(x, t)
    q3 = Function('q3')(x, t)
    Q = Matrix([[q0], [q1], [q2], [q3]])
    Q_conj = Matrix([[q1], [q0], [q3], [q2]])
    # in u_t = K(u)
    # Manakov equation:
    #   q0_t = iq0_xx+2i(q0q0*+q1q1*)q0
    #   q1_t = iq1_xx+2i(q0q0*+q1q1*)q1
    Q_t0 = I * diff(Q[0], x, 2) + 2 * I * kappa * (Q[0] * Q[1] + Q[2] * Q[3]) * Q[0]
    Q_t1 = -I * diff(Q[1], x, 2) - 2 * I * kappa * (Q[0] * Q[1] + Q[2] * Q[3]) * Q[1]
    Q_t2 = I * diff(Q[2], x, 2) + 2 * I * kappa * (Q[0] * Q[1] + Q[2] * Q[3]) * Q[2]
    Q_t3 = -I * diff(Q[3], x, 2) - 2 * I * kappa * (Q[0] * Q[1] + Q[2] * Q[3]) * Q[3]
    Q_t = Matrix([[Q_t0], [Q_t1], [Q_t2], [Q_t3]])

    psi0 = Function('psi0')(x, t)
    psi1 = Function('psi1')(x, t)
    psi2 = Function('psi2')(x, t)
    psi3 = Function('psi3')(x, t)
    Psi = Matrix([[psi0], [psi1], [psi2], [psi3]])
    L_l_Psi = Matrix([[+diff(Psi[0], x), 0, 0, 0], \
                      [0, -diff(Psi[1], x), 0, 0], \
                      [0, 0, +diff(Psi[2], x), 0], \
                      [0, 0, 0, -diff(Psi[3], x)]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def Chen_Eq2(x, t):
    a, b = symbols('a, b')
    q0 = Function('q0')(x, t)
    q1 = Function('q1')(x, t)
    Q = Matrix([[q0], [q1]])
    Q_conj = Matrix([[q1], [q0]])
    # p. 492 of "Integrability of Nonlnear Hamiltonian Systems by
    # Inverse Scattering Method", by Chen et al.
    #    a=1
    #    b=0
    Q_t0 = I * diff(Q[0], x, 2) - a * diff(Q[0] * Q[0] * Q[1], x) + b * I * Q[0] * Q[0] * Q[1]
    Q_t1 = -I * diff(Q[1], x, 2) - a * diff(Q[1] * Q[1] * Q[0], x) - b * I * Q[1] * Q[1] * Q[0]
    Q_t = Matrix([[Q_t0], [Q_t1]])

    psi0 = Function('psi0')(x, t)
    psi1 = Function('psi1')(x, t)
    Psi = Matrix([[psi0], [psi1]])
    L_l_Psi = Matrix([[diff(Psi[0], x), 0], \
                      [0, -diff(Psi[1], x)]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def Chen_Eq3(x, t):
    a, b = symbols('a, b')
    q0 = Function('q0')(x, t)
    Q = Matrix([[q0]])
    Q_conj = Matrix([q0])
    # p. 492 of "Integrability of Nonlnear Hamiltonian Systems by
    # Inverse Scattering Method", by Chen et al.
    Q_t0 = -diff(Q[0], x, 3) - 2 * diff(Q[0], x) ** 3 - a * diff(Q[0], x) * sin(2 * Q[0]) ** 2
    Q_t = Matrix([[Q_t0]])

    psi0 = Function('psi0')(x, t)
    Psi = Matrix([[psi0]])
    L_l_Psi = Matrix([[diff(Psi[0], x)]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def K_ViscousBurgers(x, t):
    q = Function('q')(x, t)
    Q = Matrix([q])
    Q_conj = Q.copy()
    # in u_t = K(u)
    # vNLSE: q_t = iq_xx-2(q*)qq_x
    alpha, beta = symbols('alpha, beta')

    psi0 = Function('psi0')(x, t)
    Psi = Matrix([[psi0]])

    q_t = beta * diff(Q[0], x, 2) - alpha * Q[0] * diff(Q[0], x)
    Q_t = Matrix([[q_t]])
    L_l_Psi = Matrix([[diff(Psi[0], x)]])
    #    L_l_Psi = Matrix( [ [ diff(Psi[0],(x,2)) ] ])  #This never works out...
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def K_KdV(x, t):
    q = Function('q')(x, t)
    Q = Matrix([q])
    Q_conj = Q
    # in u_t = K(u)
    # vNLSE: q_t = iq_xx-2(q*)qq_x
    q_t = -diff(q, x, 3) - 6 * q * diff(q, x)
    Q_t = Matrix([q_t])
    psi0 = Function('psi0')(x, t)
    Psi = Matrix([[psi0]])

    L_l_Psi = Matrix([[diff(Psi[0], x)]])
    #    L_l_Psi = Matrix( [ [ diff(Psi[0],(x,2)) ] ])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def K_KdV_O5(x, t):
    q = Function('q')(x, t)
    Q = Matrix([q])
    Q_conj = Q
    # in u_t = K(u)
    # vNLSE: q_t = iq_xx-2(q*)qq_x
    q_t = diff(q, x, 5) + 10 * q * diff(q, x, 3) + 20 * diff(q, x) * diff(q, x, 2) \
          + 30 * q ** 2 * diff(q, x)
    Q_t = Matrix([q_t])

    psi0 = Function('psi0')(x, t)
    Psi = Matrix([[psi0]])

    #    L_l_Psi = Matrix( [ [ diff(Psi[0],x) ] ])
    L_l_Psi = Matrix([[diff(Psi[0], (x, 2))]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def K_MKdV(x, t):
    q = Function('q')(x, t)
    Q = Matrix([q])
    Q_conj = Q
    # in u_t = K(u)
    # vNLSE: q_t = iq_xx-2(q*)qq_x
    q_t = -diff(q, x, 3) - 6 * q ** 2 * diff(q, x)
    Q_t = Matrix([q_t])

    psi0 = Function('psi0')(x, t)
    Psi = Matrix([[psi0]])

    #    L_l_Psi = Matrix( [ [ diff(Psi[0],x) ] ])
    L_l_Psi = Matrix([[diff(Psi[0], (x, 2))]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def K_SawadaKotera(x, t):
    # https://arxiv.org/pdf/2012.03456.pdf,
    # Integrable systems:  From the inverse spectral transform to zero curvature
    # condition,
    # Khan et al.
    q = Function('q')(x, t)
    Q = Matrix([q])
    Q_conj = Q
    # in u_t = K(u)
    # vNLSE: q_t = -q_xxxxx-5q^2q_x-5qq_xxx-5q_xq_xx
    q_t = -5 * Q[0] ** 2 * Derivative(Q[0], x) \
          - 5 * Q[0] * Derivative(Q[0], x, 3) \
          - 5 * Derivative(Q[0], x) * Derivative(Q[0], x, 2) \
          - Derivative(Q[0], x, 5)
    Q_t = Matrix([q_t])

    psi0 = Function('psi0')(x, t)
    Psi = Matrix([[psi0]])

    #    L_l_Psi = Matrix( [ [ diff(Psi[0],x) ] ])
    L_l_Psi = Matrix([[diff(Psi[0], (x, 3))]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def K_IbragimovShabat(x, t):
    q = Function('q')(x, t)
    Q = Matrix([q])
    Q_conj = Q
    # in u_t = K(u)
    # vNLSE: q_t = iq_xx-2(q*)qq_x
    q_t = diff(q, x, 3) + 3 * q ** 2 * diff(q, x, 2) + 9 * q * diff(q, x) ** 2 + 3 * q ** 4 * diff(q, x)
    Q_t = Matrix([q_t])

    psi0 = Function('psi0')(x, t)
    Psi = Matrix([[psi0]])

    L_l_Psi = Matrix([[diff(Psi[0], x)]])
    #    L_l_Psi = Matrix( [ [ diff(Psi[0],(x,2)) ] ])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;


def Chen_Eq2_variation(x, t):
    a, b = symbols('a, b')
    q0 = Function('q0')(x, t)
    q1 = Function('q1')(x, t)
    Q = Matrix([[q0], [q1]])
    Q_conj = Matrix([[q1], [q0]])
    # p. 492 of "Integrability of Nonlnear Hamiltonian Systems by
    # Inverse Scattering Method", by Chen et al.
    #    a=1
    #    b=0
    Q_t0 = I * diff(Q[0], x, 2) - diff(Q[0] * Q[0] * Q[1], x) - 2 * Q[0] * Q[1] * diff(Q[0], x)
    Q_t1 = -I * diff(Q[1], x, 2) - diff(Q[1] * Q[1] * Q[0], x) - 2 * Q[1] * Q[0] * diff(Q[1], x)
    Q_t = Matrix([[Q_t0], [Q_t1]])

    psi0 = Function('psi0')(x, t)
    psi1 = Function('psi1')(x, t)
    Psi = Matrix([[psi0], [psi1]])
    L_l_Psi = Matrix([[diff(Psi[0], x), 0],
                      [0, -diff(Psi[1], x)]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi


def K_Boussinesq(x, t):
    q0 = Function('q0')(x, t)
    q1 = Function('q1')(x, t)
    Q = Matrix([[q0], [q1]])
    Q_conj = Q.copy()
    # in u_t = K(u)
    # Bousinesq eq.:
    #   u_t = -v_x,
    #   v_t = -u_xxx+8uu_x
    #    a = 1
    a = 3
    Q_t0 = -a * diff(Q[1], x)
    Q_t1 = -diff(Q[0], x, 3) + 8 * Q[0] * diff(Q[0], x)
    Q_t = Matrix([[Q_t0], [Q_t1]])

    psi0 = Function('psi0')(x, t)
    psi1 = Function('psi1')(x, t)
    Psi = Matrix([[psi0], [psi1]])

    L_l_Psi = Matrix([[0, diff(Psi[0], x, 2)], \
                      [a * Psi[1], 0]])
    return Q_t, Q, Q_conj, Psi, L_l_Psi;
