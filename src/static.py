import numpy as np
import matplotlib.pyplot as plt

from lib import *

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +     Static solution with eneral material properties     +
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# %% Problem characteristics

### Material properties (constant)
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 1e3      # [N]
L  = 1        # [m]
nN = 3        # [-]
mu = 0.1      # [kg/m]

### Material properties (variable)
# def E(x): return 210 * (x + 1)
# def I(x): return 3.3e7
# def mu(x): return 0.1 * (1 + x/10)

### Boundaries
# BCtype = 'cantilever'
BCtype = 'fixed'

BC     = (0, 0, 0, 0)

# Mesh
grid = np.linspace(0, L, nN)
plotMesh(grid, BCtype = BCtype)

# %% CASE 1: Linear load distribution
# def q(x):
#     return k * x      

# def exact(x):
#     return (20*k*L**3*x**2 - 10*k*L**2*x**3 + k*x**5)/(120*E*I)

# RHS   = getRHS(grid, q)

# %% CASE 2: Nodal force

node   = np.array([1])   # ID of nodes where force is applied
force  = np.array([- k])   # Applied nodal forces

RHS = getPointForce(grid, node, force)

def exact(x):    
    return np.piecewise(x, [x < L/2, x>= L/2], [lambda x: (-k*x / (48 * E * I)) * (3*L**2 - 4*x**2), lambda x: (-k*(L-x) / (48 * E * I)) * (3*L**2 - 4*(L-x)**2)])

# %% Solver

S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-2] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-1] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC, BCtype)

# Solve
sol     = sparse.linalg.spsolve(Se, RHSe)
fig, _  = plotBeam(grid, sol[:-2], -1, nData = 100, BCtype = BCtype, exact = exact)
