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