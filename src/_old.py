"""
Old versions of functions 
=============================================================================
 Last reviewed:  05.2022
=============================================================================
"""
import numpy as np
from lib import *

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