# %%
import numpy as np
import matplotlib.pyplot as plt
from lib import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as Tk
import matplotlib.animation as animation

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +     Eigenmodes/superposition simulation with general material properties     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# %% Problem characteristics

# Material properties (constant)
E  = 1       # [N/mm2]
I  = 1       # [mm4]
k  = 1       # [N]
L  = 1       # [m]
nN = 10     # [-]
mu = 1       # [kg/m]

# Mesh
grid = np.linspace(0, L, nN)

# Boundary
BC   = (0, 0, 0, 0)

#get matrices
S, M  = getMatrices(grid, E, I, mu, quadrature = False)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-1] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-2] = 1.0

# Apply BCs
node   = np.array([5, -1])   # ID of nodes where force is applied
force  = np.array([k, -k/5])  # Applied nodal forces
RHS = getPointForce(grid, node, force)
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC, BCtype = "cantilever")
steadySol  = sparse.linalg.spsolve(Se, RHSe)

# Solving generalized eigenvalue problem exactly and numerically
#eigfreq_num, eigvec = eigenvalue_method(Me,6,Se)
eigfreq_num, eigvec, eigenval = eigenvalue_method_2(Me,6,Se)

print(eigenval)
plt.figure()
x_plot = np.linspace(grid.min(), grid.max(), 200)
beam = get_sol(grid, eigvec[:-2,0])
y_1 = beam(x_plot)/(np.max(np.abs(beam(x_plot))))
plt.plot(x_plot,y_1)
plt.show()

eigfreq_exact, eigfunc = eigenvalue_method_exact(grid, E, I, mu, L, 10)

#comparing exact and numerical eigenvalues 

#plt.figure()
#plt.plot(eigfunc[:,0])
#plt.show()

#plt.figure()
#plotBeam(grid,eigvec[:-2,0],-1)
#plt.show()

#Comparing the numerical and exact eigenfrequencies
plt.figure()
plt.plot(eigfreq_exact[:5],"*",label = "Exact")
plt.plot(eigfreq_num[:5],"o",label = "Numerical")
plt.xlabel("I-th eigenfrequency")
plt.legend(loc = "upper left")
plt.title("Comparison of the exact and numerical eigenfrequencies")
plt.show()

#%%
#Simulating superpositions of eigenvectors
modes = np.array([2,3]) #The mode numbers that will be in the superpositions

t_0 = 0
t_f = 10000
Nt = 1000

w_0 = steadySol
w_diff_0 = np.zeros(w_0.shape)
#w_0 = np.ones(w_0.shape)
superposition_dynamic = eigenvalue_method_dynamic(t_0,t_f,Nt,w_0,w_diff_0,Me,Se,modes,Fourier = False)
    
sol = superposition_dynamic
ylim = (-100, 100)

### Figure object definition
fig = plt.Figure(figsize = (5, 3), dpi = 150)

### Tkinter window
root   = Tk.Tk()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(column=0,row=1)

### Animation definition
ax     = fig.add_subplot(111)

nData  = 100

nN     = len(grid)
x_plot = np.linspace(grid.min(), grid.max(), nData) 

beam = get_sol(grid, sol[0:-2, 0])
line, = ax.plot(x_plot, beam(x_plot) * 1e3, color= '#808080', label = 'numerical')

ax.plot(x_plot, np.zeros(x_plot.shape), color= '#959595', linestyle = '--')
ax.axvline(x = 0, color="black", linestyle="-", linewidth = 5)

ax.set_ylabel('y-dimension (-)')
ax.set_xlabel('x-dimension (-)')

ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
    top=True, right= False, left=True, width = 1)

ax.set_ylim(ylim[0], ylim[1])

def animate(i):
    beam = get_sol(grid, sol[0:-2, i])
    line.set_ydata( beam(x_plot) * 1e3)  # update the data
    return line, 


ani = animation.FuncAnimation(fig, animate, np.arange(0, sol.shape[1]), interval = 5, blit=False)

### Main loop
root.mainloop()
# %%
