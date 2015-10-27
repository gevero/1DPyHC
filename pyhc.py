import numpy as np  # Numerical Python
import scipy as sp
import scipy.constants as sp_c  # Universal constants
# from numba import autojit,jit, double, int16
# import numdifftools as nd
# import time
# from random import Random                            # random generator
# import os
'''Contains fields and methods to setup and solve a 1D photonic crystal problem'''

# Some constants
c = sp_c.speed_of_light * (10.0 ** 9)  # in nm/s to be used with wavevector in 1/nm


# K calculates beta (k_parallel to the planes) as a
# function of wavelength (nm) and incident angle (rad)
def f_beta(n_inc, l, theta):
    '''calculates beta (k_parallel to the planes) as a function of wavelength (nm) and incident angle (rad)

    args:
    'n_inc' = refractive index of the incident medium
    'l,theta' = incident wavelength in nm and incident angle in rad'''

    return 2.0 * np.pi * n_inc * np.sin(theta) / l


# light angular frequency as a function of wavelength
def f_omega(l):
    '''calculates light angular frequency as a function of wavelength (nm)

    args:
    'l' = incident wavelength in nm'''

    return 2.0 * np.pi * c / l


# A,B,C,D, coefficients of the unit cell translation operator
def ABCD(a, b, n1, n2, omega, beta):
    '''A,B,C,D, coefficients of the unit cell translation operator

    args:
    'a,b'=Thicknesses of the 1D unit cell (b is the first thickness)
    'n1,n2'=Refractive index of the unit cell (n2 is the first)
    'omega'=light frequency in rad/s '''

    # Wavevectors in each medium
    k1 = sp.sqrt((n1 * omega / c) ** 2 - beta ** 2)
    k2 = sp.sqrt((n2 * omega / c) ** 2 - beta ** 2)

    # A,B,C,D
    a_out = np.exp(-1j * k1 * a) * (np.cos(k2 * b) - 0.5 * 1j *
                                    (k2 / k1 + k1 / k2) * np.sin(k2 * b))
    b_out = np.exp(1j * k1 * a) * (-0.5 * 1j *
                                   (k2 / k1 - k1 / k2) * np.sin(k2 * b))
    c_out = np.exp(-1j * k1 * a) * (0.5 * 1j *
                                    (k2 / k1 - k1 / k2) * np.sin(k2 * b))
    d_out = np.exp(1j * k1 * a) * (np.cos(k2 * b) + 0.5 * 1j *
                                   (k2 / k1 + k1 / k2) * np.sin(k2 * b))

    return a_out, b_out, c_out, d_out


# |0.5(A+D)-1|, its sign gives allowed and forbidden
# band, as well as band edges at |0.5(A+D)-1|
def AD(a, b, n1, n2, omega, beta):
    '''Evaluates |0.5(A+D)-1| to study the photonic crystal bands

    args:
    'a,b'=Thicknesses of the 1D unit cell (b is the first thickness)
    'n1,n2'=Refractive index of the unit cell (n2 is the first)
    'omega'=light frequency in rad/s '''

    A, B, C, D = ABCD(a, b, n1, n2, omega, beta)

    return np.abs(0.5 * (A + D)) - 1.0


# K bloch wavevector
def f_K(a, b, n1, n2, omega, beta):
    '''Evaluates the Bloch wavevector

    args:
    'a,b'=Thicknesses of the 1D unit cell (b is the first thickness)
    'n1,n2'=Refractive index of the unit cell (n2 is the first)
    'omega'=light frequency in rad/s '''

    A, B, C, D = ABCD(a, b, n1, n2, omega, beta)

    return np.arccos(0.5 * (A + D)) / (a + b)


# N-unit-cell translation operator
def f_A(a, b, n1, n2, N, omega, beta):
    '''Evaluates the N-unit-cell translation operator

    args:
    'a,b'=Thicknesses of the 1D unit cell (b is the first thickness)
    'n1,n2'=Refractive index of the unit cell (n2 is the first)
    'N'= number of unit cells
    'omega'=light frequency in rad/s '''

    # preliminary calculations
    A, B, C, D = ABCD(a, b, n1, n2, omega, beta)
    K = f_K(a, b, n1, n2, omega, beta)
    L = a + b

    # matrix elements
    A11 = A * np.sin(N * K * L) / np.sin(K * L) - np.sin(
        (N - 1) * K * L) / np.sin(K * L)
    A12 = B * np.sin(N * K * L) / np.sin(K * L)
    A21 = C * np.sin(N * K * L) / np.sin(K * L)
    A22 = D * np.sin(N * K * L) / np.sin(K * L) - np.sin(
        (N - 1) * K * L) / np.sin(K * L)

    return np.array([[A11, A12], [A21, A22]])


# incident medium additional operator
def f_M(n1, n_inc, omega, beta):
    '''Evaluates the additional operator for the incident medium

    args:
    'n1,n_inc'=Refractive index of the unit cell n1 and of the incident medium
    'omega'=light frequency in rad/s '''

    # preliminary calculations
    k1 = sp.sqrt((n1 * omega / c) ** 2 - beta ** 2)
    k_inc = sp.sqrt((n_inc * omega / c) ** 2 - beta ** 2)

    # print 'n1',n1
    # print 'n_inc',n_inc
    #
    # print 'k1^2',(n1*omega/c)**2-beta**2
    # print 'k_inc^2',(n_inc*omega/c)**2-beta**2
    # print 'k1',sp.sqrt((n1*omega/c)**2-beta**2)
    # print 'k_inc',sp.sqrt((n_inc*omega/c)**2-beta**2)
    # print 'k1',k1
    # print 'k_inc',k_inc

    # matrix elements
    M11 = 1.0 + k1 / k_inc
    M12 = 1.0 - k1 / k_inc
    M21 = 1.0 - k1 / k_inc
    M22 = 1.0 + k1 / k_inc

    # print 'M11,M12,M21,M22',M11,M12,M21,M22

    return 0.5 * np.array([[M11, M12], [M21, M22]])


# substrate medium additional operator
def f_N(a, b, n1, n_sub, omega, beta):
    '''Evaluates additional operator for the substrate

    args:
    'a,b'=Thicknesses of the 1D unit cell (b is the first thickness, a the second)
    'n1,n_sub'=Refractive index of the unit cell (n2 is the first, n1 the second) and of the substrate
    'omega'=light frequency in rad/s '''

    # preliminary calculations
    k1 = sp.sqrt((n1 * omega / c) ** 2 - beta ** 2)
    k_sub = sp.sqrt((n_sub * omega / c) ** 2 - beta ** 2)
    L = a + b

    # matrix
    N11 = sp.exp(1j * k_sub * L) * (1.0 + k1 / k_sub)
    N12 = sp.exp(1j * k_sub * L) * (1.0 - k1 / k_sub)
    N21 = sp.exp(-1j * k_sub * L) * (1.0 - k1 / k_sub)
    N22 = sp.exp(-1j * k_sub * L) * (1.0 + k1 / k_sub)

    return 0.5 * np.linalg.inv(np.array([[N11, N12], [N21, N22]]))


# total stack operator
def f_T(a, b, n1, n2, n_inc, n_sub, N, omega, beta):
    '''Evaluates additional operator for the substrate

    args:
    'a,b'=Thicknesses of the 1D unit cell (b is the first thickness, a the second)
    'n1,n2,n_inc,n_sub'=Refractive index of the unit cell (n2 is the first, n1 the second) and of incident medium and substrate
    'N'=total number of unit cells
    'omega'=light frequency in rad/s '''

    # preliminary calculations
    m_M = f_M(n1, n_inc, omega, beta)
    m_A = f_A(a, b, n1, n2, N, omega, beta)
    m_N = f_N(a, b, n1, n_sub, omega, beta)

    # print 'm_M',m_M,m_M.shape
    # print 'm_A',m_A,m_A.shape
    # print 'm_N',m_N,m_N.shape
    # print 'np.dot(m_M,m_A)',np.dot(m_M,m_A),np.dot(m_M,m_A).shape

    # matrix elements
    m_T = np.dot(np.dot(m_M, m_A), m_N)

    # print 'omega',omega
    # print 'beta', beta
    # print 'm_T',m_T

    return m_T


# reflectance for symmetric structures
def rN_sym(a, b, n1, n2, N, omega, beta):
    '''Evaluates the reflectance for a symmetric structure

    args:
    'a,b'=Thicknesses of the 1D unit cell (b is the first thickness, a the second)
    'n1,n2'=Refractive index of the unit cell (n2 is the first, n1 the second)
    'N'=total number of unit cells
    'omega'=light frequency in rad/s '''

    # preliminary calculations
    A, B, C, D = ABCD(a, b, n1, n2, omega, beta)
    K = f_K(a, b, n1, n2, omega, beta)
    L = a + b

    return (np.abs(C) ** 2) / (np.abs(C) ** 2 + np.abs(
        (np.sin(K * L) / np.sin(N * K * L))) ** 2)


# reflectance in the omega beta space
def rN(a, b, n1, n2, n_inc, n_sub, N, omega, beta):
    '''Evaluates the reflectance for a general structure in the omega beta space

    args:
    'a,b'=Thicknesses of the 1D unit cell (b is the first thickness, a the second)
    'n1,n2,n_inc,n_sub'=Refractive index of the unit cell (n2 is the first, n1 the second) and of incident medium and substrate
    'N'=total number of unit cells
    'omega'=light frequency in rad/s '''

    # preliminary calculations
    m_T = f_T(a, b, n1, n2, n_inc, n_sub, N, omega, beta)
    T21 = m_T[1, 0]
    T11 = m_T[0, 0]

    return np.abs(T21 / T11) ** 2


# reflectance in the lambda theta space
def rN_lambda_theta(a, b, n1, n2, n_inc, n_sub, N, wl, theta):
    '''Evaluates the reflectance for a general structure in the lambda theta space

    args:
    'a,b'=Thicknesses of the 1D unit cell (a is the first thickness, b the second)
    'n1,n2,n_inc,n_sub'=Refractive index of the unit cell (n1 is the first, n2 the second) and of incident medium and substrate
    'N'=total number of unit cells
    'wl'=light wavelength in nm
    'theta'=incident angle in degrees '''

    # Calculations
    m_T = f_T(b, a, n2, n1, n_inc, n_sub, N, f_omega(wl), f_beta(
        n_inc, wl, np.pi * theta / 180.0))
    T21 = m_T[1, 0]
    T11 = m_T[0, 0]

    return np.abs(T21 / T11) ** 2
