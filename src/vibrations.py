import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh

from lib import *

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +     Static solution with eneral material properties     +
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# %% Problem characteristics

# Material properties (constant)
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 1e3      # [N]
L  = 1        # [m]
nN = 30       # [-]
mu = 0.1      # [kg/m]

# Material properties (variable)
# def E(x): return 210 * (x + 1)
# def I(x): return 3.3e7
# def mu(x): return 0.1 * (1 + x/10)

# Mesh
grid = np.linspace(0, L, nN)

# Boundary
BC   = (0, 0, 0, 0)
# %% Solver

S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-1] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-2] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, np.zeros(S.shape[0]), (e0, eL), (d0, dL), BC)

# Solve
N = 2*grid.shape[0]
K = Me.shape[0] - N

a_k = np.ones(N-K)
b_k = np.ones(N-K)
#a_k[0] = 1
#b_k[0] = 1
t = 1

nat_freq, eigenmode = eigenvalue_method_exact(grid, t, E, I, mu, L, a_k, b_k, 1)
plt.plot(grid,eigenmode)
plt.show()
#eigenval, eigenmode = eigenvalue_method(N-K,Me, Se, t, a_k, b_k)
#exact_eigenval = exact_eigenvalue(E,I,mu,L,100)

