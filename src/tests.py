import numpy as np
import matplotlib.pyplot as plt

from lib import *
from _old import *

# %% Matrix computation and solver tests

"""
# +++++++++++++++++++++++++++++++++++++++++++++++
# +    Constant material properties (Vova's)    +
# +++++++++++++++++++++++++++++++++++++++++++++++

# Material properties
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 1000     # [N]
L  = 1        # [m]
nN = 30       # [-]
mu = 0.1      # [kg/m]

# Problem characteristics (mesh, BCs and applied force)
grid = np.linspace(0, L, nN)

BC   = (0, 0, 0, 0)

def q(x):
    return k * x

def exact(x):
    return (20*k*L**3*x**2 - 10*k*L**2*x**3 + k*x**5)/(120*E*I)


# Compute matrices
loc_S, loc_M = get_local_matrix()
S, M         = get_global_matrices(grid, E, I, loc_S, loc_M, mu)

RHS = get_RHS(grid,q)

e0 = np.zeros(nN*2)
e0[0]= 1.0

eL = np.zeros(nN*2)
eL[-1]= 1.0

d0 = np.zeros(nN*2)
d0[1] = 1.0

dL = np.zeros(nN*2)
dL[-2] = 1.0

# Apply BCs
Se, RHSe = fixBeam(S, RHS, (e0, eL), (d0, dL), BC)

# Solver
sol      = sparse.linalg.spsolve(Se, RHSe)

plotBeam(grid, sol[:-2], 100, exact)
"""

# +++++++++++++++++++++++++++++++++++++++++++++++
# +    General material properties (Ferran's)   +
# +++++++++++++++++++++++++++++++++++++++++++++++

# Material properties (constant)
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 1000     # [N]
L  = 1        # [m]
nN = 30       # [-]
mu = 0.1      # [kg/m]


# Material properties (variable)
# def E(x): return 210 * (x + 1)
# def I(x): return 3.3e7
# def mu(x): return 0.1 * (1 + x/10)

# Problem characteristics (mesh, BCs and applied force)
grid = np.linspace(0, L, nN)

BC   = (0, 0, 0, 0)

def q(x):
    return k * x

def exact(x):
    return (20*k*L**3*x**2 - 10*k*L**2*x**3 + k*x**5)/(120*E*I)

S, M  = getMatrices(grid, E, I, mu, quadrature = True)
RHS   = getRHS(grid, q)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-1] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-2] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC)

# Solve
sol      = sparse.linalg.spsolve(Se, RHSe)
plotBeam(grid, sol[:-2], 100, exact)

# %% Check computational cost 

import time
from lib import *

start = time.time()
sym_S, _  = getMatrices(grid, E, I, mu, quadrature = False)
print("--- Symbolic: %2.1f seconds" % (time.time() - start))

start    = time.time()
quad_S, _ = getMatrices(grid, E, I, mu, quadrature = True)
print("--- Quadrature: %.2E seconds" % (time.time() - start))

error = np.linalg.norm(sym_S.toarray() - quad_S.toarray(), ord = 2)
print("--- Error: %.2E" % (error))


#%% Timedependent problem uses the variables defined in the fixed case

u = sol
u_1 = np.zeros(np.shape(sol))
u_2 = np.copy(u_1)
h = 0.1
for i in range(100):
    u,u_1,u_2 = Newmarkmethod_step(u,u_1,u_2,h,M,Se,RHSe)
    

