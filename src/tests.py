import numpy as np
import matplotlib.pyplot as plt

from lib import *

# %% Test #1:
#
# +++++++++++++++++++++++++++++++++++++++++++++++
# +    Cantilever beam with distributed load    +
# +++++++++++++++++++++++++++++++++++++++++++++++

# Material and mesh properties

E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 1e3      # [N]
L  = 1        # [m]
nN = 30       # [-]
mu = 20       # [kg/m]

### Boundaries
BCtype = 'cantilever'
BC     = (0, 0, 0, 0)   # (w(0), w'(0), Q(L), M(L))

### Mesh
grid = np.linspace(0, L, nN)

### Load
def q(x):
    return k * x      

def analytical(x):
    return - (20*k*L**3*x**2 - 10*k*L**2*x**3 + k*x**5)/(120*E*I)

RHS   = getRHS(grid, q)

### Main solver
S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-2] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-1] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC, BCtype)

# Solve
sol     = sparse.linalg.spsolve(Se, - RHSe)

### Plot solution
# General definitions
nData  = 100
x_plot = np.linspace(grid.min(), grid.max(), nData) 
    
beam = get_sol(grid, sol[:-2])

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

fig, ax = plt.subplots(2, 1, figsize=(5, 3), dpi = 150, sharex = True)

# Load plot
ax[0].fill_between(x_plot, q(x_plot), color = 'r', alpha = 0.5, label = 'distributed load')
ax[0].plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '-')
ax[0].axvline(x=0, color="black", linestyle="-", linewidth = 5)

ax[0].legend(loc = 'lower right')

ax[0].set_ylabel('load (N)')

ax[0].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
       top=True, right= False, left=True, width = 1)
ax[0].ticklabel_format(style = 'sci', scilimits = (-1, 1))

low, high = ax[0].set_ylim()
ax[0].set_ylim(beam(x_plot)[0]*2e3 - high * 1.25, high * 1.25)

# Beam plot
ax[1].plot(x_plot, beam(x_plot) * 1e3, color= '#808080', label = 'numerical')
ax[1].plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '--')
ax[1].plot(x_plot, analytical(x_plot) * 1e3, color = 'r', linestyle = '-.', label = 'analytical')
ax[1].axvline(x=0, color="black", linestyle="-", linewidth = 5)

plt.legend(loc = 'upper right')
    
ax[1].set_xlabel('x-direction (-)')
ax[1].set_ylabel('deformation (mm)')

ax[1].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
       top=True, right= False, left=True, width = 1)

low, high = ax[1].set_ylim()
ax[1].set_ylim(low * 1.25, beam(x_plot)[0]*2e3 - low * 1.25)

plt.show()

# %% Test #2:
#
# +++++++++++++++++++++++++++++++++++++++++++++++
# +  Simply supported beam with point force     +
# +++++++++++++++++++++++++++++++++++++++++++++++

# Material and mesh properties
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 1e3      # [N]
L  = 1        # [m]
nN = 3        # [-]
mu = 0.1      # [kg/m]

### Boundaries
BCtype = 'fixed'
BC     = (0, 0, 0, 0)   # (w(0), w(L), M(0), M(L))

### Mesh
grid = np.linspace(0, L, nN)

### Load
node   = np.array([1])   # ID of nodes where force is applied
force  = np.array([-k])  # Applied nodal forces

RHS = getPointForce(grid, node, force)

def analytical(x):    
    return np.piecewise(x, [x < L/2, x>= L/2], [lambda x: (-k*x / (48 * E * I)) * (3*L**2 - 4*x**2), lambda x: (-k*(L-x) / (48 * E * I)) * (3*L**2 - 4*(L-x)**2)])

### Main solver
S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-2] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-1] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC, BCtype)

# Solve
sol     = sparse.linalg.spsolve(Se, RHSe)

### Plot solution
# General definitions
nData  = 100
x_plot = np.linspace(grid.min(), grid.max(), nData) 
    
beam = get_sol(grid, sol[:-2])

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

fig, ax = plt.subplots(2, 1, figsize=(5, 3), dpi = 150, sharex = True)

# Load plot
ax[0].plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '-')
ax[0].axvline(x=0, color="black", linestyle="-", linewidth = 5)
ax[0].axvline(x = grid.max(), color="black", linestyle="-", linewidth = 5) 

ax[0].arrow(L/2, k, 0, -k, length_includes_head = True, head_width = 0.02, head_length = 200, fc='r', ec='r', label = 'point force', zorder=10)

ax[0].legend(loc = 'lower right')

ax[0].set_ylabel('load (N)')

ax[0].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
       top=True, right= False, left=True, width = 1)
ax[0].ticklabel_format(style = 'sci', scilimits = (-1, 1))

low, high = ax[0].set_ylim()
ax[0].set_ylim(beam(x_plot)[0]*2e3 - high * 1.25, high * 1.25)

# Beam plot
ax[1].plot(x_plot, beam(x_plot) * 1e3, color= '#808080', label = 'numerical')
ax[1].plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '--')
ax[1].plot(x_plot, analytical(x_plot) * 1e3, color = 'r', linestyle = '-.', label = 'analytical')
ax[1].axvline(x=0, color="black", linestyle="-", linewidth = 5)
ax[1].axvline(x = grid.max(), color="black", linestyle="-", linewidth = 5) 

plt.legend(loc = 'upper right')
    
ax[1].set_xlabel('x-direction (-)')
ax[1].set_ylabel('deformation (mm)')

ax[1].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
       top=True, right= False, left=True, width = 1)

low, high = ax[1].set_ylim()
ax[1].set_ylim(low * 1.25, beam(x_plot)[0]*2e3 - low * 1.25)

plt.show()

# %% Test #3:
#
# +++++++++++++++++++++++++++++++++++++++++++++++
# +       Cantilever with different loads       +
# +++++++++++++++++++++++++++++++++++++++++++++++

# Material and mesh properties
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 1e3      # [N]
L  = 1        # [m]
nN = 21       # [-]
mu = 0.1      # [kg/m]

### Boundaries
BCtype = 'cantilever'
BC     = (0, 0, 0, 0)   # (w(0), w(L), M(0), M(L))

### Mesh
grid = np.linspace(0, L, nN)
plotMesh(grid)

### Loads
node   = np.array([7, -1])   # ID of nodes where force is applied
force  = np.array([5*k, -k])   # Applied nodal forces

def q(x):
    return - k/3  

def analytical(x):
    def_1 = - k * x**2 /(6 * E * I) * (3 * L - x)
    def_2 = np.piecewise(x, [x < grid[7], x>= grid[7]], [lambda x: (5*k*x**2 / (6 * E * I)) * (3*grid[7] - x), lambda x: (5*k*grid[7]**2 / (6 * E * I)) * (3*x - grid[7])])
    def_3 = - k/3 * x**2 /(24*E*I) * (6*L**2 - 4*L*x + x**2)
    
    return def_1 + def_2 + def_3

RHS_1 = getPointForce(grid, node, force)
RHS_2 = getRHS(grid, q)

RHS   = RHS_1 + RHS_2

### Main solver
S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-2] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-1] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC, BCtype)

# Solve
sol     = sparse.linalg.spsolve(Se, RHSe)

### Plot solution
# General definitions
nData  = 100
x_plot = np.linspace(grid.min(), grid.max(), nData) 
    
beam = get_sol(grid, sol[:-2])

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

fig, ax = plt.subplots(2, 1, figsize=(5, 3), dpi = 150, sharex = True)

# Load plot
ax[0].plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '-')
ax[0].axvline(x=0, color="black", linestyle="-", linewidth = 5)

ax[0].arrow(grid[7], - 2 * k, 0, 2 * k, length_includes_head = True, head_width = 0.02, head_length = 600, fc='r', ec='r', label = 'point forces', zorder=10)
ax[0].arrow(grid[-1], k, 0, - k, length_includes_head = True, head_width = 0.01, head_length = 300, fc='r', ec='r', zorder=10)
ax[0].fill_between(x_plot, - 1.2*q(x_plot), color = 'b', alpha = 0.5, label = 'distributed load')

ax[0].legend(loc = 'lower right')

ax[0].set_ylabel('load (N)')

ax[0].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
       top=True, right= False, left=True, width = 1)
ax[0].ticklabel_format(style = 'sci', scilimits = (-1, 1))

low, high = ax[0].set_ylim()
ax[0].set_ylim(low * 1.25, beam(x_plot)[0]*2e3 - low * 1.25)

# Beam plot
ax[1].plot(x_plot, beam(x_plot) * 1e3, color= '#808080', label = 'numerical')
ax[1].plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '--')
ax[1].plot(x_plot, analytical(x_plot) * 1e3, color = 'r', linestyle = '-.', label = 'analytical')
ax[1].axvline(x=0, color="black", linestyle="-", linewidth = 5)

plt.legend(loc = 'upper right')
    
ax[1].set_xlabel('x-direction (-)')
ax[1].set_ylabel('deformation (mm)')

ax[1].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
       top=True, right= False, left=True, width = 1)

low, high = ax[1].set_ylim()
ax[1].set_ylim(low * 1.25, beam(x_plot)[0]*2e3 - low * 1.25)

plt.show()

# %% Test #4:
#
# +++++++++++++++++++++++++++++++++++++++++++++++
# +        Cantilever with applied moment       +
# +++++++++++++++++++++++++++++++++++++++++++++++

# Material and mesh properties
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = -1e5     # [Nm]
L  = 1        # [m]
nN = 5        # [-]
mu = 0.1      # [kg/m]

### Boundaries
BCtype = 'fixed'
BC     = (0, 0, 0, 0)   # (w(0), w(L), M(0), M(L))

### Mesh
grid = np.linspace(0, L, nN)

### Load
node   = np.array([2])   # ID of nodes where force is applied
force  = np.array([-k])  # Applied nodal forces

RHS = getPointForce(grid, node, force, loadType = 'moment')

def analytical(x):    
    return np.piecewise(x, [x < L/2, x>= L/2], [lambda x: (k*x / (24 * L * E * I)) * (L**2 - 4 * x**2), lambda x: (-k*(L - x) / (24 * L * E * I)) * (L**2 - 4 * (L - x)**2)])

### Main solver
S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-2] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-1] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC, BCtype)

# Solve
sol     = sparse.linalg.spsolve(Se, RHSe)

### Plot solution
# General definitions
nData  = 100
x_plot = np.linspace(grid.min(), grid.max(), nData) 
    
beam = get_sol(grid, sol[:-2])

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

fig, ax = plt.subplots(2, 1, figsize=(5, 3), dpi = 150, sharex = True)

# Load plot
ax[0].plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '-', zorder = 0)
ax[0].axvline(x=0, color="black", linestyle="-", linewidth = 5)
ax[0].axvline(x = grid.max(), color="black", linestyle="-", linewidth = 5) 

drawCirc(ax[0], 0.1, L/2, 0, 60, 270, color_= 'red')
ax[0].scatter(L/2, 0, marker = 'o', color = 'r', s = 5, label = 'applied moment')

ax[0].legend(loc = 'lower right')

ax[0].set_ylabel('load $\cdot 10^{6}$ (Nm) ')

ax[0].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
       top=True, right= False, left=True, width = 1)

# low, high = ax[0].set_ylim()
# ax[0].set_ylim(beam(x_plot)[0]*2e3 - high * 1.25, high * 1.25)

ax[0].set_ylim([-0.18, 0.18])

# Beam plot
ax[1].plot(x_plot, beam(x_plot) * 1e3, color= '#808080', label = 'numerical')
ax[1].plot([grid.min(), grid.max()], [beam(x_plot)[0]*1e3, beam(x_plot)[0]*1e3], color= '#959595', linestyle= '--')
ax[1].plot(x_plot, analytical(x_plot) * 1e3, color = 'r', linestyle = '-.', label = 'analytical')
ax[1].axvline(x=0, color="black", linestyle="-", linewidth = 5)
ax[1].axvline(x = grid.max(), color="black", linestyle="-", linewidth = 5) 

plt.legend(loc = 'lower right')
    
ax[1].set_xlabel('x-direction (-)')
ax[1].set_ylabel('deformation (mm)')

ax[1].tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
       top=True, right= False, left=True, width = 1)
ax[1].ticklabel_format(style = 'sci', scilimits = (-1, 1))

low, high = ax[1].set_ylim()
ax[1].set_ylim(low * 1.25, beam(x_plot)[0]*2e3 - low * 1.25)

plt.show()
