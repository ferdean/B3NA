from multiprocessing.sharedctypes import Value
import numpy as np
import matplotlib.pyplot as plt


from lib import *
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

# %% CASE 1: Linear load distribution
#def q(x):
#    return k * x      

#RHS   = getRHS(grid, q)

# %% CASE 2: Nodal force

# node   = np.array([10, -1])   # ID of nodes where force is applied
# force  = np.array([k, -k/5])  # Applied nodal forces

# RHS = getPointForce(grid, node, force)

# %% CASE 3: Timedependent load distribution

def q(x,t):
    return k * x * np.exp(-10*t)    

RHS   = getRHS(grid, q)

# %% Get steady solution (works as initial conditions) 

S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-1] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-2] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC, "cantilever")

# Solve
#steadySol      = sparse.linalg.spsolve(Se, RHSe) #non timedependent forcing
steadySol      = sparse.linalg.spsolve(Se, RHSe(0)) #timedependent forcing


# %% Time simulation

u_1_0 = np.zeros(steadySol.shape)
u_2_0 = np.zeros(steadySol.shape)

initialConds = (steadySol, u_1_0, u_2_0)

# Simulation characteristics
#RHSe  = np.zeros(steadySol.shape)   # Free vibration case

h     = 1e-3
t0    = 0.0
T     = 10.0

sol, time     = newmarkMethod(Me, Se, RHSe, initialConds, h, t0, T, verbose = False)

# %% Generate animation

import matplotlib.animation as animation

nData  = 100

nN     = len(grid)
x_plot = np.linspace(grid.min(), grid.max(), nData) 
ylim   = (-2e-5, 2e-5)


plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

fig = plt.figure(figsize=(5, 3), dpi = 150)
ax  = fig.add_subplot(111)

# fig, ax = plt.subplots()

def animation_frame(i): 
    ax.clear()
    beam = get_sol(grid, sol[0:-2, i])
    
    ax.plot(x_plot, beam(x_plot) * 1e3, color= '#808080', label = 'numerical')
    ax.plot([grid.min(), grid.max()], [0, 0], color= '#959595', linestyle= '--')
    
    
    ax.axvline(x=0, color="black", linestyle="-", linewidth = 5)
    
    ax.set_xlabel('x-direction (-)')
    ax.set_ylabel('deformation (mm)')
    
    ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
        top=True, right= False, left=True, width = 1)
    
    # low, high = plt.ylim()
    # bound = max(abs(low), abs(high))
    # plt.ylim(-bound, bound)
    
    plt.ylim(ylim[0], ylim[1])
    
    plt.title('t = %.2f s'%(time[i]))
    
    plt.text(0.875, 0.425,'undeformed', ha='center', va='center', transform=ax.transAxes, color= '#959595')
    
    plt.legend(loc = 'lower left')
    
    return fig

ani = animation.FuncAnimation(fig, animation_frame,np.arange(0, 300), interval=10)

ani.save('temporal_2.gif', writer='imagemagick', fps= 50)