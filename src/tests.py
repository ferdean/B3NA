import numpy as np
import matplotlib.pyplot as plt

from lib import *


# %% Shape functions validation
"""
colors = ('#fcfdbf', '#fed395', '#fea973', '#fa7d5e', '#e95462', '#c83e73',
          '#a3307e', '#7e2482', '#59157e', '#331067', '#120d31', '#000004')

x_plot =  np.linspace(0, 1, 1000)
grid = np.linspace(0, 1, 3)

fig, ax = plt.subplots(2, 1, dpi = 330)

ax[0].scatter(grid, np.zeros(grid.shape), marker= 'x', color = 'k')
ax[1].scatter(grid, np.zeros(grid.shape), marker= 'x', color = 'k')

for i in range(3):
    
    phi = get_phi(grid, i, derivative = 2)
    
    ax[0].plot(x_plot, phi(x_plot)[0], color= colors[i])
    ax[1].plot(x_plot, phi(x_plot)[1], color= colors[i])
    

ax[0].set_title('Odd shape functions (position)')
ax[1].set_title('Even shape functions (derivative)')

plt.tight_layout()
plt.show()

# %% Beam plotting test

grid = np.linspace(0, 1, 5)
#solution = np.zeros((2*len(grid),))

#plotBeam(grid, solution)

# %% Matrix computation and solver test - not working yet
"""
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 10       # [N]
L  = 1        # [m]
nN = 30       # [-]
mu = 0.1      # [kg/m]

grid = np.linspace(0, L, nN)

BC = (0, 0, 0, 0)

def q(x):
    return k * x

def exact(x):
    return (20*k*L**3*x**2 - 10*k*L**2*x**3 + k*x**5)/(120*E*I)

#S, RHS, (e0, eL), (d0, dL) = computeMatrices(grid, q, E, I, n_quad = 50)


loc_S,loc_M = get_local_matrix()
S,M = get_global_matrices(grid, E, I, loc_S, loc_M, mu)

RHS = get_RHS(grid,q)
e0 = np.zeros(nN*2)
e0[0]= 1.0
eL = np.zeros(nN*2)
eL[-1]= 1.0
d0 = np.zeros(nN*2)
d0[1] = 1.0
dL = np.zeros(nN*2)
dL[-2] = 1.0

Se, RHSe = fixBeam(S, RHS, (e0, eL), (d0, dL), BC)

sol      = sparse.linalg.spsolve(Se, RHSe)

plotBeam(grid, sol[:-2], 100, exact)
#print(sol)

#%% Timedependent problem uses the variables defined in the fixed case

u = sol
u_1 = np.zeros(np.shape(sol))
u_2 = np.copy(u_1)
h = 0.1
for i in range(100):
    u,u_1,u_2 = Newmarkmethod_step(u,u_1,u_2,h,M,Se,RHSe)
    
# %%
