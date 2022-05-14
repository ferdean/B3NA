"""
Bending of Bernoulli beams project (numerical analysis) library of functions.
=============================================================================
 Last reviewed:  05.2022
=============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import fixed_quad
from scipy import sparse
import scipy as sp
import sympy as sm
def get_phi(grid, i, derivative = 0):
    """
    Computes the functions that form the ansatz space Vh, where Vh is a space 
    of pecewise cubic polynomials.
    
    Input
    ----------
    grid: {array}
        Vector with gridpoints.
    i: {int}
        Index of the considered node in the grid.

    Returns
    -------
    phi: {function}
        Basis functions
        
    """
        
    h = grid[1:] - grid[:-1]
    
    def SF(xi, derivative = derivative):
        """
        Generates the cubic shape functions.

        Parameters
        ----------
        xi: {vector, int}
            (Normalized) absolute coordinate in he grid

        Returns
        -------
        phi_bar: {tuple}
            Collection of shape functions

        """
        
        if derivative == 0:
            phi_1 = 1 - 3 * xi**2 + 2*xi**3
            phi_2 = xi * (xi - 1)**2
            phi_3 = 3 * xi**2 - 2 * xi**3
            phi_4 = xi**2 * (xi - 1)
        
        elif derivative == 1:
            phi_1 = -6 * xi + 6 * xi**2
            phi_2 = 3 * xi**2 - 4 * xi + 1
            phi_3 = -6 * (xi - 1) * xi
            phi_4 = xi * (3 * xi - 2)
            
        elif derivative == 2:
            phi_1 = -6 + 12 * xi
            phi_2 = 6 * xi - 4
            phi_3 = 6 - 12 * xi
            phi_4 = 6 * xi - 2
            
        phi_bar = (phi_1, phi_2, phi_3, phi_4)
        
        return phi_bar

    
    def phi(y):
        """
        Computes the basis functions defined inside the ansatz space

        Parameters
        ----------
        y: {vector, int}
            Absolute position in the grid.

        Returns
        -------
        funct_odd: {vector, int}
            Basis function related with position.
        funct_even: {vector, int}
            Basis function related with derivative.

        """
        
        if type(y) == int or float: 
            y = np.array(y)
            
        funct_odd  = np.zeros(y.shape)
        funct_even = np.zeros(y.shape)   
    
        if i == 0:
            funct_odd[ (y >= grid[0]) & (y <= grid[1])] = SF(y[ (y >= grid[0]) & (y <= grid[1])]/h[0])[0]
            funct_even[(y >= grid[0]) & (y <= grid[1])] = h[0] * SF(y[ (y >= grid[0]) & (y <= grid[1])]/h[0])[1]
           
        elif i == grid.shape[0]-1:
            funct_odd[ (y >= grid[-2]) & (y <= grid[-1])] = SF((y[ (y >= grid[-2]) & (y <= grid[-1])] - grid[-2])/h[0])[2]
            funct_even[(y >= grid[-2]) & (y <= grid[-1])] = h[0] * SF((y[ (y >= grid[-2]) & (y <= grid[-1])] - grid[-2])/h[0])[3]
           
        else:
            funct_odd[ (y>= grid[i - 1]) & (y <= grid[i])] = SF((y[ (y>= grid[i - 1]) & (y <= grid[i])] - grid[i - 1])/h[i])[2]
            funct_odd[ (y>= grid[i]) & (y <= grid[i + 1])] = SF((y[ (y>= grid[i]) & (y <= grid[i + 1])] - grid[i])/h[i])[0]
            
            funct_even[(y>= grid[i - 1]) & (y <= grid[i])] = h[i] * SF((y[(y>= grid[i - 1]) & (y <= grid[i])] - grid[i - 1])/h[i])[3]
            funct_even[(y>= grid[i]) & (y <= grid[i + 1])] = h[i] * SF((y[(y>= grid[i]) & (y <= grid[i + 1])] - grid[i])/h[i])[1]
    
        
        return funct_odd, funct_even
    
    return phi


def get_sol(grid, coeffs):
    """
    Returns the solution function

    Parameters
    ----------
    grid: {array}
        Vector with gridpoints.
    coeffs: {array}
        Solution vector (the odd positions represent deformation
                         and the even ones, derivative).
    Returns
    -------
    w: {function}
        Solution function.
    """
    
    # This for loop will be removed using the matrices ex and dx from computeMatrices()
    def w(y):
        beam = np.zeros(y.shape)
        
        for idx in range(int(coeffs.shape[0]/2)):
            phi   = get_phi(grid, idx)
            beam += coeffs[2 * idx] * phi(y)[0] + coeffs[2 * idx + 1] * phi(y)[1]

        return beam
    
    return w

def plotBeam(grid, coeffs, nData, *argv):
    """
    Plots the deformed beam

    Parameters
    ----------
    grid: {array}
        Vector with gridpoints.
    coeffs: {array}
        Solution vector (the odd positions represent deformation
                         and the even ones, derivative).
    nData: {int, optional}
        Number of plotting datapoints. The default is 200.

    Returns
    -------
    None.

    """
    nN     = len(grid)
    x_plot = np.linspace(grid.min(), grid.max(), nData) 
    
    beam = get_sol(grid, coeffs)
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'font.size' : 9})
    
    fig, ax = plt.subplots(figsize=(5, 3), dpi = 200)
    
    ax.plot(x_plot, beam(x_plot) * 1e3, color= '#808080')
    ax.plot([grid.min(), grid.max()], [0, 0], color= '#959595', linestyle= '--')
    
    for arg in argv: 
        ax.plot(x_plot, arg(x_plot) * 1e3, color = 'r', linestyle = '-.')
    
    ax.axvline(x=0, color="black", linestyle="-", linewidth = 5)
    
    ax.set_xlabel('x-direction (-)')
    ax.set_ylabel('deformation (mm)')
    
    ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
        top=True, right= False, left=True, width = 1)
    
    low, high = plt.ylim()
    bound = max(abs(low), abs(high))
    plt.ylim(-bound, bound)
    
    plt.text(0.875, 0.425,'undeformed', ha='center', va='center', transform=ax.transAxes, color= '#959595')
    plt.show()

def computeMatrices(grid, q, E, I, n_quad = 40):
    """
    Stiffness, inhomogeneity and physical/natural BC matrices computation
    
    Parameters
    ----------
    grid: {array}
        Vector with gridpoints.
    q:  {function}
        Right-hand-side of the differential equation.
    E: {function} or {scalar}
        Young modulus [N/mm2]
    I: {function} or {scalar}
        Area moment of inertia.
    n_quad: {scalar, optional}
        Order of the quadrature rule (numerical integration). Default is 40.

    Returns
    -------
    S: {array}
        Stiffness matrix.
    RHS: {array}
        Right-hand-side matrix
    ex: {array}
        Vector of basis functions
    dx: {array}
        Vector of derivatives of basis functions

    """
        
    if callable(E): 
        constantprops = False 
    else: 
        constantprops = True;
    
    N   = grid.shape[0]
    L   = grid[-1]            # in 1D grid
    S   = np.zeros((2*N, 2*N))
    RHS = np.zeros((2*N,))
    
    e0  = np.zeros((2*N,))
    eL  = np.zeros((2*N,))
    d0  = np.zeros((2*N,))
    dL  = np.zeros((2*N,))
       
    # Implementation with for loops (should be avoided) and unoptimized (dense matrices). 
    # Will be updated in further versions:
    for j in range(0, 2*N, 2):
        
        phi_j         = get_phi(grid, j//2, derivative = 0)
        phi_j_prime   = get_phi(grid, j//2, derivative = 1)
                
        RHS[j], _     = fixed_quad(lambda x: q(x) * phi_j(x)[0], 0, L, n = n_quad)
        RHS[j + 1], _ = fixed_quad(lambda x: q(x) * phi_j(x)[1], 0, L, n = n_quad)
        
        e0[j]         = phi_j(0)[0];   eL[j]      = phi_j(L)[0]
        e0[j + 1]     = phi_j(0)[1];   eL[j + 1]  = phi_j(L)[1]
        
        d0[j]         = phi_j_prime(0)[0];   dL[j]      = phi_j_prime(L)[0]
        d0[j + 1]     = phi_j_prime(0)[1];   dL[j + 1]  = phi_j_prime(L)[1]
        
        for k in range(0, 2*N, 2):
            
            phi_j = get_phi(grid, j//2, derivative = 2)
            phi_k = get_phi(grid, k//2, derivative = 2)
                        
            if constantprops: 
                S[j, k], _         = fixed_quad(lambda x: E * I * phi_j(x)[0] * phi_k(x)[0], 0, L, n = n_quad)
                S[j + 1, k], _     = fixed_quad(lambda x: E * I * phi_j(x)[1] * phi_k(x)[0], 0, L, n = n_quad)
                S[j, k + 1], _     = fixed_quad(lambda x: E * I * phi_j(x)[0] * phi_k(x)[1], 0, L, n = n_quad)
                S[j + 1, k + 1], _ = fixed_quad(lambda x: E * I * phi_j(x)[1] * phi_k(x)[1], 0, L, n = n_quad)
                
            else: 
                S[j, k], _         = fixed_quad(lambda x: E(x) * I(x) * phi_j(x)[0] * phi_k(x)[0], 0, L, n = n_quad)
                S[j + 1, k], _     = fixed_quad(lambda x: E(x) * I(x) * phi_j(x)[1] * phi_k(x)[0], 0, L, n = n_quad)
                S[j, k + 1], _     = fixed_quad(lambda x: E(x) * I(x) * phi_j(x)[0] * phi_k(x)[1], 0, L, n = n_quad)
                S[j + 1, k + 1], _ = fixed_quad(lambda x: E(x) * I(x) * phi_j(x)[1] * phi_k(x)[1], 0, L, n = n_quad)
                
    return sparse.csr_matrix(S), RHS, (e0, eL), (d0, dL)

def get_local_matrix():
    # <><><><><><>
    x = sm.Symbol('x')
    print("done")
    S = np.zeros((4,4))
    phi1 = 1 - 3* x**2 + 2 * x**3
    phi2 = x* (x-1)**2
    phi3 = 3*x**2 - 2*x**3
    phi4 = x**2 * (x-1)
    phi  = [phi1,phi2,phi3,phi4]
    for i in range(4):
        for j in range(4):
            S[i,j] = float(sm.integrate(sm.diff(phi[i],x,2)*sm.diff(phi[j],x,2),(x,0,1)))
    #a = 
    #a = sm.diff(phi1,x,4)
    #print(S)
    return S


def get_global_matrices(grid, E, I, loc_S):

    N = grid.shape[0]*2
    S = np.zeros((N,N))
    #S = loc_S/(h**3)
    first_mode = np.array([ [1,0,1,0],
                            [0,0,0,0],
                            [1,0,1,0],
                            [0,0,0,0] ])

    second_mode = np.array([ [0,1,0,1],
                             [1,0,1,0],
                             [0,1,0,1],
                             [1,0,1,0] ])

    third_mode = np.array([ [0,0,0,0],
                            [0,1,0,1],
                            [0,0,0,0],
                            [0,1,0,1] ])

    
    h_array = grid[1:] - grid[0:-1]
    h_odd  = np.take(h_array, np.arange(0,len(h_array),2))
    h_even = np.take(h_array, np.arange(1,len(h_array),2))
    A = np.kron(np.diag(h_odd**-3),first_mode) + np.kron(np.diag(h_odd**-2),second_mode) + np.kron(np.diag(h_odd**-1),third_mode)
    A = A * np.kron(np.eye(len(h_odd)),loc_S)

    B = np.kron(np.diag(h_even**-3),first_mode) + np.kron(np.diag(h_even**-2),second_mode) + np.kron(np.diag(h_even**-1),third_mode)
    B = B * np.kron(np.eye(len(h_even)),loc_S)
    S[0:A.shape[0],0:A.shape[0]] += A
    S[2:B.shape[0]+2,2:B.shape[0]+2] += B
    #print(h_odd)
    #print(h_even)
    #print(S)
    return S*E*I

def get_RHS(grid,q):  
    # q(x) is a function that has to work with numpy arrays
    # Everything here is done in the matrix vector multiplication form(not super obvious, but less for loops)
    # Since the test functions are 4th order, I decided to do 4-point gaus quadrature, not to waste precision
     
    N = len(grid)
    h = grid[1:] - grid[0:-1]
    # weights and nodes of Gaus quadrature for 4 points (scaled to integrate from 0 to 1)(original values from wikipedia)
    # https://en.wikipedia.org/wiki/Gaussian_quadrature
    nodes = (np.sort(np.array([(3/7 - (2/7)*(6/5)**0.5)**0.5, -(3/7 - (2/7)*(6/5)**0.5)**0.5 , 
                              (3/7 + (2/7)*(6/5)**0.5)**0.5, -(3/7 + (2/7)*(6/5)**0.5)**0.5])) +1) /2
    weights = np.array([(18-30**0.5)/36, (18+30**0.5)/36, (18+30**0.5)/36, (18-30**0.5)/36]) /2

    # basis functions on the unit interval
    phis = [lambda x :1-3*x*x+2*x*x*x, lambda x :x*(x-1)**2, lambda x :3*x*x - 2*x*x*x, lambda x :x*x*(x-1)]
    # initialize and fill local matrix
    s = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            s[i,j] = weights[j] * phis[i](nodes[j]) 
    # assemble global matrix
    G = np.zeros((N*2,(N-1)*4))
    for i in range(N-1):
        G[2*i:2*i+4, 4*i:4*i+4] = s

    # locations on the beam where q has to be evaluated ( 4 points per element)
    q_nodes = np.tile(nodes,(N-1,1)) * np.tile(h,(4,1)).T + np.tile(grid[:-1],(4,1)).T
    # compute q at those locations and scale each q value by the element length( the element at which this node happened to be)
    q_vec = (q(q_nodes) * np.tile(h,(4,1)).T).flatten() 

    LHS  = G @ q_vec
    #np.set_printoptions(precision=3)
    #print(LHS)
    return LHS


def fixBeam(S, RHS, e, d, BC):
    
    
    e0, eL = e
    d0, dL = d
    a, b, QL, ML = BC
    
    basisVector = sparse.csr_matrix(np.vstack((e0, d0)))
        
    RHSe = RHS + QL * eL + ML * dL
    RHSe = np.hstack((RHSe, np.array([a, b])))
    
    Se   = sparse.vstack((S, basisVector))
    Se   = sparse.hstack((Se, (sparse.hstack((basisVector, sparse.csr_matrix((2, 2))))).T))
    
    return Se, RHSe

