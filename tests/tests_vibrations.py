# %%
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('C:/Users/Ferhat/Desktop/TUBerlin/Project numerical analysis/beam-num-analysis/src') #This only works on my laptop 
#, to fix it we need to make a __init__.py file...
from lib import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as Tk
import matplotlib.animation as animation

# +++++++++++++++++++++++++++++++++
# +     1D Cantilever Problem     +
# +++++++++++++++++++++++++++++++++

# %% Problem characteristics

# Material properties (constant)
E  = 1       # [N/mm2]
I  = 1       # [mm4]
k  = 1       # [N]
L  = 1       # [m]
nN = 100     # [-]
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
Me, Se, RHSe = fixBeam(M, S, np.zeros(S.shape[0]), (e0, eL), (d0, dL), BC, BCtype = "cantilever")

# Solving generalized eigenvalue problem exactly and numerically
n = 10 #number of eigenfreq/vectors
eigfreq_num, eigvec = eigenvalue_method(Me,n,Se)
eigfreq_exact, eigfunc = eigenvalue_method_exact(grid, E, I, mu, L, n)

#Comparing the numerical and exact eigenfrequencies
plt.figure()
plt.plot(eigfreq_exact,"*",label = "Exact")
plt.plot(eigfreq_num,"o",label = "Numerical")
plt.xlabel("j-th Eigenfrequency")
plt.legend(loc = "upper left")
plt.title("Comparison of the exact and numerical eigenfrequencies")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(np.arange(0, (eigfreq_num.shape)[0], step=1))
plt.grid(linestyle='-.', linewidth=0.7)
plt.tight_layout()
#plt.savefig("Exact_Numerical_Eigenfrequency.png",dpi = 300)
plt.show()

#comparing eigenfunctions and eigenvectors 

plotBeam(grid,eigvec[:-2,0],-1,flag = False,scaling=True)
plt.plot(eigfunc[:,0]/np.max(np.abs(eigfunc[:,0])),"--")
plt.show()

# %%
