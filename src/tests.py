import numpy as np
import matplotlib.pyplot as plt

from lib import *

# %% Matrix computation and solver tests

# +++++++++++++++++++++++++++++++++++++++++++++++
# +    Constant material properties (Vova's)    +
# +++++++++++++++++++++++++++++++++++++++++++++++

# Material properties
E  = 210      # [N/mm2]
I  = 3.3e7    # [mm4]
k  = 1000     # [N]
L  = 1        # [m]
nN = 30       # [-]
mu = 20       # [kg/m]

# Problem characteristics (mesh, BCs and applied force)
grid = np.linspace(0, L, nN)

BC   = (0, 0, 0, 0)

def q(x):
    return k * x

def exact(x):
    return (20*k*L**3*x**2 - 10*k*L**2*x**3 + k*x**5)/(120*E*I)

S, M  = getMatrices(grid, E, I, mu, quadrature = True)

RHS   = getRHS(grid, q)
# RHS  = getPointForce(grid, [2, 10, 15], [-k, 5*k, -2.5*k])



e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-1] = 1.0
d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-2] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC)

# Solve
sol      = sparse.linalg.spsolve(Se, RHSe)
figura = plotBeam(grid, sol[:-2], 100, -1)

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

        
# %% Newmark method

# Initial conditions
u_0   = sol
u_1_0 = np.zeros(sol.shape)
u_2_0 = np.zeros(sol.shape)

initialConds = (u_0, u_1_0, u_2_0)

# Simulation characteristics
RHSe  = np.zeros(sol.shape) 
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

ani = animation.FuncAnimation(fig, animation_frame, interval=10)

# ani.save('temporal_2.gif', writer='imagemagick', fps= 50)

# %% Save clicks

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

try:
    import Tkinter as tk
except:
    import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

### Figure object definition
fig = plt.Figure(figsize = (5, 3), dpi = 150)
ax     = fig.add_subplot(111)

### Tkinter window
root   = Tk.Tk()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(column=0,row=1)

### Function definition
coords = []

def on_click(event):
    if event.inaxes is not None:
        # print(event.xdata, event.ydata)
        global coords
        # coords.append((event.xdata, event.ydata))
        
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return array[idx]
        
        x = find_nearest(grid, event.xdata)
        print(x)
        
        plt.clf()
        ax.scatter(x, np.sin(x), color = 'red', s= 5)
        
    else:
        print('Clicked ouside axes bounds but inside plot window')
    

### Plot definition

x_plot = np.linspace(0, 2*np.pi, 500)
y_plot = np.sin(x_plot) # Toy frame

line, = ax.plot(x_plot, y_plot, color= '#808080')

ax.plot(x_plot, np.zeros(x_plot.shape), color= '#959595', linestyle = '--')
ax.axvline(x = 0, color="black", linestyle="-", linewidth = 5)

ax.set_ylabel('y-dimension (-)')
ax.set_xlabel('x-dimension (-)')

ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
    top=True, right= False, left=True, width = 1)

clicks = fig.canvas.mpl_connect('button_press_event', on_click)
plt.show()


### Some postprocess
grid = np.linspace(0, 6, 5)

### Main loop
root.mainloop()

