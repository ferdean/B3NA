import numpy as np
import matplotlib.pyplot as plt
from lib import *

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +     Eigenmodes/superposition simulation with general material properties     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# %% Problem characteristics

# Material properties (constant)
E  = 1     # [N/mm2]
I  = 1       # [mm4]
k  = 1       # [N]
L  = 1       # [m]
nN = 30      # [-]
mu = 1       # [kg/m]

# Mesh
grid = np.linspace(0, L, nN)

# Boundary
BC   = (0, 0, 0, 0)

#get matrices
S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-1] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-2] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, np.zeros(S.shape[0]), (e0, eL), (d0, dL), BC)

# Solving generalized eigenvalue problem exactly and numerically
from scipy.sparse.linalg import eigsh

eigfreq_num, eigvec = eigenvalue_method(Me,Se)
eigfreq_exact, eigfunc = eigenvalue_method_exact(grid, E, I, mu, L, 10)

#comparing eigenfunction of exact problem with the eigenvector of discretized problem
plotBeam(grid, eigvec[:-2,2], 100, -1)
plt.figure()
plt.plot(grid,eigfunc[:,2])
plt.show()


plt.figure()
plt.plot(eigfreq_exact[:5],"*",label = "exact")
plt.plot(eigfreq_num[:5],"o",label = "numerical")
plt.xlabel("ith eigenfrequency")
plt.legend(loc = "upper right")
plt.show()

#Simulating superpositions of eigenvectors
modes = np.zeros(6)
modes[4] = 1e-3 #activated modes
modes[5] = 1e-3 #activated modes
t_0 = 0
t_f = 100
Nt = 100

superposition_dynamic = eigenvalue_method_dynamic(t_0,t_f,Nt,Me,Se,modes)

print(superposition_dynamic.shape)

#still need to do the animation for now I have subplots at fixed timestamps :))
for i in range(4):
    plotBeam(grid,superposition_dynamic[:-2,i*25],100,-1)
