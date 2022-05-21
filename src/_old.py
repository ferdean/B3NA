"""
Old versions of functions 
=============================================================================
 Last reviewed:  05.2022
=============================================================================
"""
import numpy as np
from lib import *

# +++++++++++++++++++++++++++++++++++++++++++++++
# +    Constant material properties (Vova's)    +
# +++++++++++++++++++++++++++++++++++++++++++++++

def get_local_matrix(verbose = False):
    """
    Computes the local stiffness matrix on an isoparameterized element of length
    equal to 1. Works for constant material properties.     
    """
    # <><><><><><>
      
    xi = sm.Symbol('xi')
    
    if verbose: print("done")
    
    S = np.zeros((4,4))
    
    phi1 = 1 - 3 * xi**2 + 2 * xi**3
    phi2 = xi * (xi - 1)**2
    phi3 = 3 * xi**2 - 2 * xi**3
    phi4 = xi**2 * (xi - 1)
    phi  = [phi1, phi2, phi3, phi4]
    
    for i in range(4):
        for j in range(4):
            S[i,j] = float(sm.integrate(sm.diff(phi[i], xi, 2)*sm.diff(phi[j], xi ,2), (xi, 0, 1)))
 
    M = np.zeros((4,4))

    for i in range(4):
        for j in range(4):
            M[i,j] = float(sm.integrate(phi[i]*phi[j], (xi, 0, 1)))
            
    return S, M

def get_global_matrices(grid, E, I, loc_S, loc_M, mu):  
    """
    Computes the local stiffness matrix on an isoparameterized element of length
    equal to 1. Works for constant material properties.    
    """
    N = grid.shape[0]*2
    S = np.zeros((N,N))
    M = np.zeros((N,N))
   
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
    
    A = np.kron(np.diag(h_odd**-3),first_mode) + np.kron(np.diag(h_odd**-2),second_mode) + np.kron(np.diag(h_odd**-1), third_mode)
    A = A * np.kron(np.eye(len(h_odd)), loc_S)

    B = np.kron(np.diag(h_even**-3),first_mode) + np.kron(np.diag(h_even**-2),second_mode) + np.kron(np.diag(h_even**-1),third_mode)
    B = B * np.kron(np.eye(len(h_even)), loc_S)
    
    S[0:A.shape[0],0:A.shape[0]] += A
    S[2:B.shape[0]+2,2:B.shape[0]+2] += B
    
    # defining M matrix for timedependent case
    A = np.kron(np.diag(h_odd),first_mode) + np.kron(np.diag(h_odd**2),second_mode) + np.kron(np.diag(h_odd**3),third_mode)
    A = A * np.kron(np.eye(len(h_odd)),loc_M)

    B = np.kron(np.diag(h_even),first_mode) + np.kron(np.diag(h_even**2),second_mode) + np.kron(np.diag(h_even**3),third_mode)
    B = B * np.kron(np.eye(len(h_even)),loc_M)
    M[0:A.shape[0],0:A.shape[0]] += A
    M[2:B.shape[0]+2,2:B.shape[0]+2] += B    

    return S*E*I, mu*M