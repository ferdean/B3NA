import numpy as np
import matplotlib.pyplot as plt
from torch import eig

from lib import *

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +     Eigenmodes/superposition simulation with general material properties     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# %% Problem characteristics

# Material properties (constant)
E  = 10     # [N/mm2]
I  = 1       # [mm4]
k  = 1       # [N]
L  = 1       # [m]
nN = 30      # [-]
mu = 1       # [kg/m]

# Material properties (variable)
# def E(x): return 210 * (x + 1)
# def I(x): return 3.3e7
# def mu(x): return 0.1 * (1 + x/10)

# Mesh
grid = np.linspace(0, L, nN)

#++++++++++++++++++++++++++++++++++++++++++++
#+   Cantilever beam vibrational analysis   +
#++++++++++++++++++++++++++++++++++++++++++++

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

# Exact eigenfreq and eigenfunctions
eigfreq_exact, eigfunc_exact = eigenvalue_method_exact(grid, E, I, mu, L, 1)
plt.figure()
plt.plot(grid,eigfunc_exact)
plt.show()

#Numerical eigenfreq and eigenfunctions
eigfreq_numerical, eigvec_numerical = eigenvalue_method(N-K,Me, Se)

plotBeam(grid, eigvec_numerical[:-2,0], 100, -1)

plt.figure()
plt.plot(eigfreq_exact[:5],"*",label = "exact")
plt.plot(eigfreq_numerical[:5],"o",label = "numerical")
plt.xlabel("ith eigenfrequency")
plt.legend(loc = "upper right")
plt.show()