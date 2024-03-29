"""
Bending of Bernoulli beams project (numerical analysis) library of functions.
=============================================================================
 Last reviewed:  05.2022
=============================================================================
"""

# from types import NoneType
import numpy as np
from numpy import radians as rad
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, RegularPolygon
from scipy.integrate import fixed_quad
from scipy import sparse
from scipy.sparse.linalg import eigsh
from numpy.linalg import inv
from numpy.linalg import eig
import scipy as sp
import sympy as sm


# ++++++++++++++++++++++++++++++++++++++++++++++
# +             BASIS FUNCTIONS                +
# ++++++++++++++++++++++++++++++++++++++++++++++
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

def getGaussParams():
    nodes = (np.sort(np.array([(3/7 - (2/7)*(6/5)**0.5)**0.5, -(3/7 - (2/7)*(6/5)**0.5)**0.5 , 
                              (3/7 + (2/7)*(6/5)**0.5)**0.5, -(3/7 + (2/7)*(6/5)**0.5)**0.5])) +1) /2
    weights = np.array([(18-30**0.5)/36, (18+30**0.5)/36, (18+30**0.5)/36, (18-30**0.5)/36]) /2
    
    return nodes, weights


# ++++++++++++++++++++++++++++++++++++++++++++++
# +             MATRICES COMPUT.               +
# ++++++++++++++++++++++++++++++++++++++++++++++
    
    
def getMatrices(grid, E, I, mu, quadrature = True):
    """
    Computes mass and stiffness matrices for a 1D problem. 

    Parameters
    ----------
    grid: {array}
        Vector with gridpoints.
    E: {function} or {scalar}
        Young modulus [N/mm2]
    I: {function} or {scalar}
        Area moment of inertia.
    mu: {function} or {scalar}
        Density.
    quadrature: {boolean, optional}
        States if the integration is performed numercially (True) or 
        anallytically (False). The default is True.

    Returns
    -------
    S: {array}
        Stiffness matrix.
    M: {array}
        Mass matrix.

    """      
    def getLocalMatrices(idx):
        
        if callable(E): 
            constantprops = False 
        else: 
            constantprops = True;
        
        x0 = grid[idx]
        x1 = grid[idx + 1]
        
        S = np.zeros((4,4))        
        M = np.zeros((4,4))
        
        h  = (x1 - x0)
                    
        first_mode = np.array([ [1,0,1,0],
                                [0,0,0,0],
                                [1,0,1,0],
                                [0,0,0,0] ])  * h**-3

        second_mode = np.array([ [0,1,0,1],
                                 [1,0,1,0],
                                 [0,1,0,1],
                                 [1,0,1,0] ]) * h**-2
    
        third_mode = np.array([ [0,0,0,0],
                                [0,1,0,1],
                                [0,0,0,0],
                                [0,1,0,1] ]) * h**-1
        
        h_array = first_mode + second_mode + third_mode
            
        if not quadrature:
            xi = sm.Symbol('xi')
            x  = (x1 - x0) * xi + x0

            phi1 = 1 - 3 * xi**2 + 2 * xi**3
            phi2 = xi * (xi - 1)**2
            phi3 = 3 * xi**2 - 2 * xi**3
            phi4 = xi**2 * (xi - 1)
            phi  = [phi1, phi2, phi3, phi4]
    
            if constantprops:
                for i in range(4):
                    for j in range(4):
                    
                        S[i,j] = float(sm.integrate(sm.diff(phi[i], xi, 2)*sm.diff(phi[j], xi ,2), (xi, 0, 1)))
                        M[i,j] = float(sm.integrate(phi[i]*phi[j],                                 (xi, 0, 1))) 
                                        
                S = S*E*I
                M = mu*M
                        
            else:
                for i in range(4):
                    for j in range(4):
                        
                        # TBD: Check if we need to multiply with Jacobian (x1 - x0)
                        
                        S[i,j] = float(sm.integrate(E(x) * I(x) * sm.diff(phi[i], xi, 2)*sm.diff(phi[j], xi ,2), (xi, 0, 1)))
                        M[i,j] = float(sm.integrate(mu(x) * phi[i]*phi[j],                                       (xi, 0, 1)))      

        else:
                             
            phi   = [lambda x :1-3*x*x+2*x*x*x, lambda x :x*(x-1)**2, lambda x :3*x*x - 2*x*x*x, lambda x :x*x*(x-1)]
            ddphi = [lambda x :-6 + 12 * x, lambda x :6 * x - 4, lambda x :6 - 12 * x, lambda x :6 * x - 2]
            
            nodes, weights = getGaussParams()
            
            if constantprops: 
                 for i in range(4):
                    for j in range(4):
                        for k in range(4):
                            S[i,j] += weights[k] * E * I * ddphi[i](nodes[k])  * ddphi[j](nodes[k]) 
                            M[i,j] += weights[k] * mu * phi[j](nodes[k]) 
            else:             
                for i in range(4):
                    for j in range(4):
                        for k in range(4):
                            S[i,j] += weights[k] * E((x1 - x0) * nodes[k] + x0) * I((x1 - x0) * nodes[k] + x0) * ddphi[i](nodes[k])  * ddphi[j](nodes[k])
                            M[i,j] += weights[k] * mu(nodes[k]) * phi[i](nodes[k])  * phi[j](nodes[k])
                        
        return np.multiply(S, h_array), np.multiply(M, h_array)*h**4
        
    # ++++++++++++++++++++++++++++++++++++++++++++++
    # +            Function starts here            +
    # ++++++++++++++++++++++++++++++++++++++++++++++
       
    nN    = grid.shape[0]  # Number of nodes
    nE    = nN - 1         # Number of elements
    nNe   = 2              # Number of nodes per element
    nDOFn = 2              # Number of DOF per node
    nDOFg = nDOFn * nN     # Global number of DOF
    nDOFe = nDOFn * nNe    # Number of DOF per element
    
    S = np.zeros((nDOFg, nDOFg))
    M = np.zeros((nDOFg, nDOFg))
    
    Ig = np.zeros((nDOFe**2 + (nE - 1)*nDOFe**2,))
    Jg = np.zeros((nDOFe**2 + (nE - 1)*nDOFe**2,))
    Mg = np.zeros((nDOFe**2 + (nE - 1)*nDOFe**2,))
    Sg = np.zeros((nDOFe**2 + (nE - 1)*nDOFe**2,))
    
    top = np.repeat(range(nN), 2)[1:-1]  # Topology matrix
    top = top.reshape((-1, nNe))
        
    for idx in range(nE):
        e_k = np.array([top[idx, 0] + idx, top[idx, 0] + idx + 1, 
                        top[idx, 0] + idx + 2, top[idx, 0] + idx + 3])  # Local DOF
        
        S_loc, M_loc = getLocalMatrices(idx)
        
        # Matrix fast assembly [ref: F. Cuvelier, et al, arXiv preprint arXiv:1305.3122]     
        t     = np.array(range(nDOFe))
        Tt, T = np.meshgrid(t, t)
        
        ii  = T.flatten()
        jj  = Tt.flatten()
        kk  = np.array(range(nDOFe**2))
        kkk = kk + nDOFe**2 * (idx)
        
        Ig[kkk] = e_k[ii]
        Jg[kkk] = e_k[jj]
        Mg[kkk] = M_loc.flatten() 
        Sg[kkk] = S_loc.flatten()
    
    M = sp.sparse.csr_matrix((Mg, (Ig, Jg)), shape = (nDOFg, nDOFg))
    S = sp.sparse.csr_matrix((Sg, (Ig, Jg)), shape = (nDOFg, nDOFg))
       
    return S, M 


def getRHS(grid, q, t = None): 
    """
    Computes inhomogeneity. Everything is done in the matrix-vector multiplication 
    form (not super obvious, but less for loops and increased efficiency). Since 
    the test functions are 4th order, it has been decided to do 4-point gauss 
    quadrature, not to waste precision

    Parameters
    ----------
    grid: {array}
        Vector with gridpoints.
    q: {function}
        Applied force. It has to work with numpy arrays

    Returns
    -------
    RHS: {vector}
        Right-hand-side.

    """
    
    N = len(grid)
    h = grid[1:] - grid[0:-1]
    
    nodes, weights = getGaussParams() # 4-point Gauss quadrature (https://en.wikipedia.org/wiki/Gaussian_quadrature)

    # Basis functions on the unit interval
    phis = [lambda x :1-3*x*x+2*x*x*x, lambda x :x*(x-1)**2, lambda x :3*x*x - 2*x*x*x, lambda x :x*x*(x-1)]
    
    # Initialize and fill local matrix
    s = np.zeros((4,4))
    
    for i in range(4):
        for j in range(4):
            s[i,j] = weights[j] * phis[i](nodes[j]) 
            
    # Global matrix assembly
    G = np.zeros((N*2,(N-1)*4))
    
    for i in range(N-1):
        G[2*i:2*i+4, 4*i:4*i+4] = s

    # Locations on the beam where q has to be evaluated (4 points per element)
    q_nodes = np.tile(nodes,(N-1,1)) * np.tile(h,(4,1)).T + np.tile(grid[:-1],(4,1)).T
    
    # Compute q at those locations and scale each q value by the element length(the element at which this node happened to be)
    print(t)
    try:
        q_vec = (q(q_nodes) * np.tile(h,(4,1)).T).flatten() 
    except TypeError:
        def RHS(t):
            q_vec = (q(q_nodes,t) * np.tile(h,(4,1)).T).flatten()
            return G@q_vec
        return RHS

    RHS  = G @ q_vec
    
    return RHS


def getPointForce(grid, nodeID, forces, loadType = 'force'):
    
    nN    = grid.shape[0]  # Number of nodes
    nDOFn = 2              # Number of DOF per node
    nDOFg = nDOFn * nN     # Global number of DOF
    nNe   = 2
    
    top = np.repeat(range(nN), 2)[1:-1]  # Topology matrix
    top = top.reshape((-1, nNe))
    
    RHS = np.zeros((nDOFg,))
    
    if loadType == 'force':
        RHS[np.multiply(2, nodeID)] = forces
    
    elif loadType == 'moment':
        RHS[np.multiply(2, nodeID) + 1] = forces
        
    else: 
        raise ValueError("Not implemented type of load")
             
    return RHS


def fixBeam(M, S, RHS, e, d, BC, BCtype):
    """
    Applies boundary conditions to the problem

    Parameters
    ----------
    S: {array}
        Free stiffness matrix.
    RHS: {vector} 
        Free right-hand-side of the equation (inhomogeneities).
    e: {array}
        Basis functions evaluated at x = 0 and x = L.
    d: {array} 
        Derivatives of basis functions evaluated at x = 0 and x = L.
    BC: {Tuple}
        Boundary conditions, such that
            *BC[0] = Deformation at x = 0 (essential bc)
            *BC[1] = Derivative of the deformation at x = 0 (essential bc)
            *BC[2] = Shear force at x = L (physical bc)
            *BC[3] = Bending moment at x = L (physical bc)
    BCtype: {string}
        Type of beam:
            * 'cantilever'
            * 'fixed'

    Returns
    -------
    Me: {array}
        Constrained mass matrix.
    Se: {array}
        Constrained stiffness matrix.
    RHSe: {vector}
        Constrained right-hand-side.
    """   
    
    nDOF = M.shape[0]
    
    e0, eL = e
    d0, dL = d
    
    if BCtype == 'cantilever':
    
        a, b, QL, ML = BC
        
        basisVector = sparse.csr_matrix(np.vstack((e0, d0)))
                
        try:
            RHSe = RHS + QL * eL + ML * dL
            RHSe = np.hstack((RHSe, np.array([a, b])))
        except TypeError:
            def RHSe(t):
                RHSe = RHS(t) + QL * eL + ML * dL
                return np.hstack((RHSe, np.array([a, b]))) 
                
    elif BCtype == 'fixed':
        
        a0, aL, M0, ML = BC
        
        basisVector = sparse.csr_matrix(np.vstack((e0, -eL)))
        
        try:
            RHSe = RHS - M0 * d0 + ML * dL
            RHSe = np.hstack((RHSe, np.array([a0, -aL])))
        except TypeError:
            def RHSe(t):
                RHSe = RHS(t) - M0 * d0 + ML * dL
                return np.hstack((RHSe, np.array([a0, -aL])))

    else:
        raise ValueError("Not implemented type of beam")
           
    Se   = sparse.vstack((S, basisVector))
    Se   = sparse.hstack((Se, (sparse.hstack((basisVector, sparse.csr_matrix((2, 2))))).T))
    
    Me   = sparse.vstack((sparse.hstack((M, sparse.csr_matrix((nDOF, 2)))), sparse.csr_matrix((2, nDOF + 2))))   
    
    return Me, Se, RHSe


# ++++++++++++++++++++++++++++++++++++++++++++++
# +               POSTPROCESS                  +
# ++++++++++++++++++++++++++++++++++++++++++++++


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
    
    def w(y):
        beam = np.zeros(y.shape)
        
        for idx in range(int(coeffs.shape[0]/2)):
            phi   = get_phi(grid, idx)
            beam += coeffs[2 * idx] * phi(y)[0] + coeffs[2 * idx + 1] * phi(y)[1]

        return beam
    
    return w


def plotBeam(grid, coeffs, ylim, nData = 200, BCtype = 'cantilever', exact = None):
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
    ylim: {tuple}
    BCtype: {string}
        Type of beam:
            * 'cantilever'
            * 'fixed'
    Returns
    -------
    None.

    """
    x_plot = np.linspace(grid.min(), grid.max(), nData) 
    
    beam = get_sol(grid, coeffs)
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'font.size' : 9})
    
    fig, ax = plt.subplots(figsize=(5, 3), dpi = 150)
    
    ax.plot(x_plot, beam(x_plot) * 1e3, color= '#808080', label = 'numerical')
    ax.plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '--')
    
    if exact != None:
        ax.plot(x_plot, exact(x_plot) * 1e3, color = 'r', linestyle = '-.', label = 'analytical')
        plt.legend(loc = 'lower left')
    
    ax.axvline(x=0, color="black", linestyle="-", linewidth = 5)
    
    if BCtype == 'fixed':
        ax.axvline(x = grid.max(), color="black", linestyle="-", linewidth = 5) 
        
    ax.set_xlabel('x-direction (-)')
    ax.set_ylabel('deformation (mm)')
    
    ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
        top=True, right= False, left=True, width = 1)
    
    if ylim == -1:
        low, high = plt.ylim()
        bound = np.argmax((abs(low), abs(high)))
        if bound == 1:
            plt.ylim(beam(x_plot)[0]*2e3 - high * 1.25, high * 1.25)
            bound = high * 1.25
        else:
            plt.ylim(low * 1.25, beam(x_plot)[0]*2e3 - low * 1.25)
            bound = beam(x_plot)[0]*2e3 - low * 1.25
    else:
        plt.ylim(ylim[0], ylim[1])
        bound = ylim[1]
    
    
    plt.show()
    
    return fig, bound


def plotMesh(grid, nData = 100, BCtype = 'cantilever'):
    
    if BCtype != 'cantilever' and BCtype != 'fixed':
        raise ValueError("Not implemented type of beam")
         
    yData  = np.zeros(grid.shape)
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'font.size' : 9})
    
    color = ["black"]
    color = color * grid.size
    color[0] = 'white'
            
    fig, ax = plt.subplots(figsize=(5, 3), dpi = 150)
    
    if BCtype == 'fixed':
        ax.axvline(x = grid.max(), color="black", linestyle="-", linewidth = 5) 
        color[-1] = 'white'
   
    ax.plot([grid.min(), grid.max()], [0, 0], color= '#959595', linestyle= '-', zorder = 1)  
    ax.scatter(grid, yData, c = color, marker = 'x', s = 10, alpha = 1, zorder = 10) 
    ax.axvline(x=0, color="black", linestyle="-", linewidth = 5) 
        
    ax.set_xlabel('x-direction (-)')
    ax.set_ylabel('deformation (mm)')
    
    ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
        top=True, right= False, left=True, width = 1)

    plt.show()
    
    return fig


def drawCirc(ax, radius, centX, centY, angle_, theta2_, color_='black'):
    arc = Arc([centX,centY],radius,radius,angle=angle_,
          theta1=0,theta2=theta2_,capstyle='round',linestyle='-',lw=1.5,color=color_)
    ax.add_patch(arc)

    endX=centX+(radius/2)*np.cos(rad(theta2_+angle_)) #Do trig to determine end position
    endY=centY+(radius/2)*np.sin(rad(theta2_+angle_))

    ax.add_patch(                    #Create triangle as arrow head
        RegularPolygon(
            (endX, endY),            # (x,y)
            3,                       # number of vertices
            radius/9,                # radius
            rad(angle_+theta2_),     # orientation
            color=color_
        )
    )


# ++++++++++++++++++++++++++++++++++++++++++++++
# +              NEWMARK METHOD                +
# ++++++++++++++++++++++++++++++++++++++++++++++

def Newmarkmethod_step(u, u_1, u_2, h, M, S, p, beta = 1/4, gamma = 1/2, t = None):
    """
    Calculates one iterate of the Newmark method

    Parameters
    ----------
    u: {array}
        Vector with the initial value.
    u_1: {array}
        Vector with the initial value for the first derivative w.r.t. time.
    u_2: {array}
        Vector with the initial value for the second derivative w.r.t. time.
    h: {float}
        Time stepsize.
    M: {array}
        global Mass matrix.
    S: {array}
        global Stiffness matrix.
    p: {vector}
        Right-hand-side.

    Returns
    -------
    u: {array}
        Vector with the value at time = initial time + h.
    u_1: {array}
        Vector with the value at time = initial time + h for the first derivative w.r.t. time.
    u_2: {array}
        Vector with the value at time = initial time + h for the second derivative w.r.t. time.

    """

    u_star = u + u_1*h+(0.5-beta)*u_2*h**2
    u_1_star = u_1 + (1-gamma)*u_2*h
    
    A   = M + beta*h**2*S
    try: 
        b   = p - S@u_star
    except TypeError:
        b   = p(t) - S@u_star

    u_2 = sparse.linalg.spsolve(A, b)

    u   = u_star + beta*h**2*u_2
    u_1 = u_1_star + gamma*h*u_2

    return u, u_1, u_2

def newmarkMethod(M, S, RHSe, initialConds, h, t0, T, verbose = False):
    
    nS    = int((T - t0)//h)
    time  = np.linspace(t0, T, nS)
    
    u   = initialConds[0]
    u_1 = initialConds[1]
    u_2 = initialConds[2]
    
    sol = np.zeros((u.shape[0], nS))
    
    for idx in range(nS):
        u, u_1, u_2 = Newmarkmethod_step(u, u_1, u_2, h, M, S, RHSe, beta = 1/4, gamma = 1/2,t = time[idx]+h)
        sol[:, idx] = u
        
        if verbose:
            print("Epoch: " + str(idx + 1) +"/" + str(nS))
        
    return sol, time

# +++++++++++++++++++++++++++
# +    EIGENVALUE METHOD    +
# +++++++++++++++++++++++++++

def eigenvalue_method(Me,Num,Se):
    """
    OLD FUNCTION DO NOT USE...
    Calculates and sorts the first N eigenvalues and eigenmodes of the 
    generalized eigenvalue problem Me x = lambda Se x

    Parameters
    ----------
    Me: {array}
        Matrix in the left hand side of the generalized eigenvalue problem
    Se: {array}
        Matrix in the right hand side of the generalized eigenvalue problem
    Num: {integer}
        number of eigenfrequencies and eigenmodes. 

    Returns
    -------
    eigfreq: {array}
        Vector of size N with the eigenfrequencies.
    eigvec: {array}
        matrix of size N by two times the size of the grid + 4 (or +2 need to find this out :( ),  
        the column i gives the coefficient/weights for the shape functions of eigenmode i.
    """

    eigval, eigvec = eigsh(Me,Num,M = Se)
    
    idx = eigval.argsort()[::-1]   
    eigval = eigval[idx]
    eigvec = eigvec[:,idx]
    eigfreq = 1/np.sqrt(eigval)
    return eigfreq,eigvec

def eigenvalue_method_2(Me,N,Se):
    """
    Calculates and sorts (from low eigenfreq to high eigenfreq) the first N eigenvalues and eigenmodes of the 
    generalized eigenvalue problem Me x = lambda Se x

    Parameters
    ----------
    Me: {array}
        Matrix in the left hand side of the generalized eigenvalue problem
    Se: {array}
        Matrix in the right hand side of the generalized eigenvalue problem
    N: {integer}
        number of eigenfrequencies and eigenmodes to be calculated. 

    Returns
    -------
    eigfreq: {array}
        Vector of size N with the eigenfrequencies.
    eigvec: {array}
        matrix of size N by two times the size of the grid + 2,  
        the column i gives the coefficient/weights for the shape functions of eigenmode i.
    """

    A = inv(np.array(Se))@np.array(Me)
    eigval, eigvec = eig(A)
    eigval = eigval.real
    eigvec = eigvec.real

    #deleting the four eigenvalues with corresponding eigenvector which are zero
    idx = eigval > 10e-16
    eigval = eigval[idx]
    eigvec = eigvec[:,idx]

    #sorting the eigenvalues with corresponding eigenvector from big to small     
    idx = eigval.argsort()[::-1]   
    eigval = eigval[idx]
    eigvec = eigvec[:,idx]
    eigfreq = 1/np.sqrt(eigval)
    return eigfreq[:N],eigvec[:,:N],eigval

def eigenvalue_method_exact(grid, E, I, mu, L, N, BCtype = "Cantilever"):
    """
    Calculates the first N eigenvalues and eigenmodes of the cantilever and supported beam problem exactly.

    Parameters
    ----------
    grid: {array}
        Vector with gridpoints.
    E: {function} or {scalar}
        Young modulus [N/mm2]
    I: {function} or {scalar}
        Area moment of inertia.
    mu: {function} or {scalar}
        Density.
    L: {scalar}
        Length of the beam.
    N: {integer}
        number of frequencies calculated and number of first N eigenmodes superpositioned. (maybe introduce two seperate integers for this)
    BCtype: {string}
        Specify which beam problem to be solved

    Returns
    -------
    omega_j: {array}
        Vector with N natural frequencies.
    w_x_t: {array}
        the first N eigenmodes evaluated on the grid.    
    
    """
    if BCtype == "Cantilever":
        j = np.linspace(1,N,N)
        x_j = (j - 0.5)*np.pi
        if N >= 0:
            x_j[0] = 1.8751
        if N >= 1:
            x_j[1] = 4.6941
        if N >= 2:
            x_j[2] = 7.8548
        k_j = x_j/L
        eigfreq = np.sqrt(E*I/mu)*k_j**2

        def w_j(k_j,x_j,x):
            return 1/np.sqrt(L)*(np.cosh(k_j*x)-np.cos(k_j*x) - (np.cosh(x_j)+np.cos(x_j))/(np.sinh(x_j)+np.sin(x_j))*(np.sinh(k_j*x) - np.sin(k_j*x)))
        
        eigfuncs = np.zeros((grid.shape[0],N))
        for i in range(N):
            eigfuncs[:,i] = w_j(k_j[i],x_j[i],grid)

        return eigfreq,eigfuncs
        
    elif BCtype == "fixed":
        j = np.linspace(1,N,N)
        k_j = j*np.pi/L
        eigfreq = np.sqrt(E*I/mu)*k_j**2

        def w_j(k_j,x):
            return 1/np.sqrt(L)*np.sin(k_j*x)
        
        eigfuncs = np.zeros((grid.shape[0],N))
        for i in range(N):
            eigfuncs[:,i] = w_j(k_j[i],grid)

        return eigfreq,eigfuncs
    else:
        return "not a known BCtype"

def eigenvalue_method_dynamic(t_0,t_f,Nt,w_0,w_diff_0,M,S,modes,Fourier = True):

    """
    If Fourier = True then the eigenvaluemethod (i.e. all modes are taken into accordance) is computed. 
    else only the modes specified are taken into account and the weights are manually set to 1 for all modes.

    Parameters
    ----------
    t_0:{scalar} 
        initial time for the simulation
    t_f:{scalar}
        final time for the simulation
    Nt: {integer}
        Number of timesteps.
    Me: {array} 
        Matrix at the left hand side of the generalized eigenvalue problem (i.e. extended mass matrix Me).
    Se: {array} 
        Matrix at right hand side of the generalized eigenvalue problem (i.e.extended stiffnes matrix Se).
    Modes: {array}
        Array that contains the modenumber taken into accordance. Note that if Fourier is true then this argument is absolete.

    Returns
    -------
    superposition_t: {array}
        Matrix of size Nt by two times the size of the grid + 2,     
        the column i gives the coefficient/weights for the shape functions at time t = i*(t_f-t_0)/Nt.
    
    """

    Num = np.max(modes)
    a_k = np.zeros(np.max(modes))

    for i in modes:
        a_k[i-1] = 1

    b_k = np.copy(a_k)

    w_k,eigvec,_ = eigenvalue_method_2(M,Num,S)

    if Fourier:
        a_k_star = np.diag((((eigvec.T)@M)@(np.array([w_0,]*Num).T))/(((eigvec.T)@M)@eigvec))
        b_k_star = np.diag((((eigvec.T)@M)@(np.array([w_diff_0,]*Num).T))/(((eigvec.T)@M)@eigvec))
        a_k *= a_k_star
        b_k *= b_k_star

    dt = (t_f - t_0)/Nt

    superposition_t = np.zeros((M.shape[0],Nt))
    
    def superposition(t):
        return ((a_k*np.cos(w_k*t)+b_k/w_k*np.sin(w_k*t))*eigvec).sum(axis = 1)

    for i in range(Nt):
        superposition_t[:,i] = superposition(t_0+i*dt)
    
    return superposition_t


# +++++++++++++++++++++++++++++++++
# +    2D FRAME AND STRUCTURES    +
# +++++++++++++++++++++++++++++++++


def getLongitudinalMatrices(grid, E, A):
    """
    Computes mass and stiffness matrices for a 1D problem in the longitudinal
    axis (version of getMatrices()).
    
    In further versions this function will be merged with getMatrices().
    """
    
    S_loc = E * A * np.array([[1, -1], [-1, 1]]) 
    
    nN    = grid.shape[0]  # Number of nodes
    nE    = nN - 1         # Number of elements
    nNe   = 2              # Number of nodes per element
    nDOFn = 1              # Number of DOF per node       
    nDOFg = nDOFn * nN     # Global number of DOF
    nDOFe = nDOFn * nNe    # Number of DOF per element
    
    S = np.zeros((nDOFg, nDOFg))
    # M = np.zeros((nDOFg, nDOFg))
    
    Ig = np.zeros((nDOFe**2 + (nE - 1)*nDOFe**2,))
    Jg = np.zeros((nDOFe**2 + (nE - 1)*nDOFe**2,))
    # Mg = np.zeros((nDOFe**2 + (nE - 1)*nDOFe**2,))
    Sg = np.zeros((nDOFe**2 + (nE - 1)*nDOFe**2,))
    
    top = np.repeat(range(nN), 2)[1:-1]  # Topology matrix
    top = top.reshape((-1, nNe))
        
    for idx in range(nE):
        
        e_k   = top[idx]
        
        # Matrix fast assembly [ref: F. Cuvelier, et al, arXiv preprint arXiv:1305.3122]     
        t     = np.array(range(nDOFe))
        Tt, T = np.meshgrid(t, t)
        
        ii  = T.flatten()
        jj  = Tt.flatten()
        kk  = np.array(range(nDOFe**2))
        kkk = kk + nDOFe**2 * (idx)
        
        Ig[kkk] = e_k[ii]
        Jg[kkk] = e_k[jj]
        # Mg[kkk] = M_loc.flatten() 
        Sg[kkk] = S_loc.flatten()
    
    # M = sp.sparse.csr_matrix((Mg, (Ig, Jg)), shape = (nDOFg, nDOFg))
    S = sp.sparse.csr_matrix((Sg, (Ig, Jg)), shape = (nDOFg, nDOFg))
       
    return S


def getLongitudinalDef(grid, E, A, f):
    """
    Returns the solution function in a longitudinal problem

    Parameters
    ----------
    grid: {array}
        Vector with gridpoints.
    E: {function} or {scalar}
        Young modulus [N/mm2]
    A: {scalar}
        Section area of the beam.
    f: {array}
        Applied force.

    Returns
    -------
    u: {function}
        Solution function.
    """
    
    from scipy.interpolate import interp1d
    
    S = getLongitudinalMatrices(grid, E, A)
    S = np.delete(np.delete(S.toarray(), 0, 1), 0, 0)    
    S = sparse.csr_matrix(S)

    nN = grid.size 
    u  = sparse.linalg.spsolve(S, f)
    
    for idx in range(1, nN):
        u[idx - 1] = u[idx - 1] + grid[idx]
        
    u = np.insert(u, 0, grid[0], axis=0)
    
    u = interp1d(grid, u)
    
    return u

def rotateBeam(locDef, L, x0, theta):    
    """
    Returns the rotated solution function

    Parameters
    ----------
    locDef: {tuple}
        Contains the deformation functions in both directions.
            * v(x): longitudinal deformation
            * w(x): transversal deformation
    L: {scalar}
        Beam length.
    x0: {array}
        (x,y) position of the initial point.
    theta: {scalar}
        Rotation angle.

    Returns
    -------
    rotbeam: {function}
        Parametric representation of the solution.
    """
    v, w = locDef
    R    = np.array([[np.cos(theta), - np.sin(theta)],
                    [np.sin(theta),   np.cos(theta)]])

    def rotbeam(x):
        beta    = np.array([v(x) + x, w(x)])      
        globDef = x0.reshape((-1,1)) + R @ beta
        
        return globDef
    
    return rotbeam



