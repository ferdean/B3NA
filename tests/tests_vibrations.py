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


# Material and mesh properties
E  = 1    # [N/mm2]
I  = 1     # [mm4]
k  = 1      # [N]
L  = 1        # [m]
nN = 100      # [-]
mu = 1        # [kg/m]

### Boundaries
BCtype = 'cantilever'
BC     = (0, 0, 0, 0)   # (w(0), w(L), M(0), M(L))

### Mesh
grid = np.linspace(0, L, nN)

### Main solver
S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-2] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-1] = 1.0

# Apply BCs
Me, Se, _ = fixBeam(M, S,  np.zeros(S.shape[0]), (e0, eL), (d0, dL), BC, BCtype)

# Solving generalized eigenvalue problem exactly and numerically
n = 10 #number of eigenfreq/vectors
eigfreq_num, eigvec,_ = eigenvalue_method_2(Me,n,Se)
x_plot = np.linspace(grid.min(), grid.max(), 200)
eigfreq_exact, eigfunc = eigenvalue_method_exact(x_plot, E, I, mu, L, n)

#Comparing the numerical and exact eigenfrequencies
plt.figure()
plt.plot(eigfreq_exact,"*",color= 'red',label = "Exact")
plt.plot(eigfreq_num,"o",color= '#808080',label = "Numerical")
plt.xlabel("j-th Eigenfrequency")
plt.legend(loc = "upper left")
plt.title("Comparison of the exact and numerical eigenfrequencies")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(np.arange(0, (eigfreq_num.shape)[0], step=1))
plt.grid(linestyle='-.', linewidth=0.7)
plt.tight_layout()
plt.savefig("Exact_Numerical_Eigenfrequency_1D_cantilever_beam.png",dpi = 300)
plt.show()

#comparing eigenfunctions and eigenvectors 
n = 4
fig, ax = plt.subplots(int(n/2), 2)
k = 0
for i in range(int(n/2)):
    for j in range(2):
        beam = get_sol(grid, eigvec[:-2,k])
        y_1 = beam(x_plot)/(np.max(np.abs(beam(x_plot))))
        y_2 = eigfunc[:,k]/np.max(np.abs(eigfunc[:,]))
        if np.linalg.norm(y_1-y_2) < np.linalg.norm(y_1+y_2):
            sign = 1
        else:
            sign = -1

        ax[i, j].plot(x_plot, y_1, color= '#808080', label = 'eigenvector')
        ax[i, j].plot(x_plot, sign*y_2,"--", color= 'red', label = 'eigenfunction')
        ax[i, j].set_title('i = '+str(k))
        if k == 1:
            ax[0][1].legend(loc = (1.05,0.75))
        k+=1

ax[0][0].set_xticklabels([])
ax[0][1].set_xticklabels([])
ax[0][0].set_yticks(np.arange(-1,1,step = 0.5))
ax[1][0].set_yticks(np.arange(-1,1,step = 0.5))
ax[0][1].set_yticklabels([])
ax[1][1].set_yticklabels([])

fig.suptitle("Comparison of the eigenfunctions and eigenvectors")
fig.supxlabel("x-direction(-)")
fig.supylabel("Deformation (mm)")
fig.tight_layout()
plt.savefig("Eigenfunction_Eigenvector_1D_cantilever_beam.png",dpi = 300)
plt.show()

# %%
# ++++++++++++++++++++++++++++++++++++++++++++
# +     1D Simply supported beam Problem     +
# ++++++++++++++++++++++++++++++++++++++++++++

# %% Problem characteristics

# Material and mesh properties
E  = 1     # [N/mm2]
I  = 1     # [mm4]
k  = 1      # [N]
L  = 1        # [m]
nN = 100      # [-]
mu = 1        # [kg/m]

### Boundaries
BCtype = 'fixed'
BC     = (0, 0, 0, 0)   # (w(0), w(L), M(0), M(L))

### Mesh
grid = np.linspace(0, L, nN)

### Main solver
S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-2] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-1] = 1.0

# Apply BCs
Me, Se, _ = fixBeam(M, S,  np.zeros(S.shape[0]), (e0, eL), (d0, dL), BC, BCtype)

# Solving generalized eigenvalue problem exactly and numerically
n = 10 #number of eigenfreq/vectors
eigfreq_num, eigvec,_ = eigenvalue_method_2(Me,n,Se)
x_plot = np.linspace(grid.min(), grid.max(), 200)
eigfreq_exact, eigfunc = eigenvalue_method_exact(x_plot, E, I, mu, L, n,BCtype = "fixed")

#Comparing the numerical and exact eigenfrequencies
plt.figure()
plt.plot(eigfreq_exact,"*",color= 'red',label = "Exact")
plt.plot(eigfreq_num,"o",color= '#808080',label = "Numerical")
plt.xlabel("j-th Eigenfrequency")
plt.legend(loc = "upper left")
plt.title("Comparison of the exact and numerical eigenfrequencies")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xticks(np.arange(0, (eigfreq_num.shape)[0], step=1))
plt.grid(linestyle='-.', linewidth=0.7)
plt.tight_layout()
plt.savefig("Exact_Numerical_Eigenfrequency_1D_supported_beam.png",dpi = 300)
plt.show()

#comparing eigenfunctions and eigenvectors 
n = 4
fig, ax = plt.subplots(int(n/2), 2)
k = 0
for i in range(int(n/2)):
    for j in range(2):
        beam = get_sol(grid, eigvec[:-2,k])
        y_1 = beam(x_plot)/(np.max(np.abs(beam(x_plot))))
        y_2 = eigfunc[:,k]/np.max(np.abs(eigfunc[:,]))
        if np.linalg.norm(y_1-y_2) < np.linalg.norm(y_1+y_2):
            sign = 1
        else:
            sign = -1

        ax[i, j].plot(x_plot, y_1, color= '#808080', label = 'eigenvector')
        ax[i, j].plot(x_plot,sign*y_2,"--", color= 'red', label = 'eigenfunction')
        ax[i, j].set_title('i = '+str(k))
        if k == 1:
            ax[0][1].legend(loc = (1.05,0.75))
        k+=1

ax[0][0].set_xticklabels([])
ax[0][1].set_xticklabels([])
ax[0][0].set_yticks(np.arange(-1,1,step = 0.5))
ax[1][0].set_yticks(np.arange(-1,1,step = 0.5))
ax[0][1].set_yticklabels([])
ax[1][1].set_yticklabels([])

fig.suptitle("Comparison of the eigenfunctions and eigenvectors")
fig.supxlabel("x-direction(-)")
fig.supylabel("Deformation (mm)")
fig.tight_layout()
plt.savefig("Eigenfunction_Eigenvector_1D_supported_beam.png",dpi = 300)
plt.show()

#%%
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +     1D Cantilever Problem movies of eigenmodes and superpositions     +
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# %% Problem characteristics

# Material properties (constant)
E  = 1       # [N/mm2]
I  = 1       # [mm4]
k  = 1       # [N]
L  = 1       # [m]
nN = 10      # [-]
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

#Simulating superpositions of eigenvectors
modes = np.array([4]) #The mode numbers that will be in the superpositions

t_0 = 0
t_f = 1
Nt = 30*10

w_0 = steadySol
w_diff_0 = np.zeros(w_0.shape)
#w_0 = np.ones(w_0.shape)
superposition_dynamic = eigenvalue_method_dynamic(t_0,t_f,Nt,w_0,w_diff_0,Me,Se,modes,Fourier = False)
    
sol = superposition_dynamic
ylim = (-10, 10)

### Figure object definition
fig = plt.Figure(figsize = (5, 3), dpi = 150)

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

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

ax.set_ylabel('deformation (mm)')
ax.set_xlabel('x-dimension (-)')

ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
    top=True, right= False, left=True, width = 1)

ax.set_ylim(ylim[0], ylim[1])

def animate(i):
    beam = get_sol(grid, sol[0:-2, i])
    line.set_ydata( beam(x_plot) * 1e3)  # update the data
    return line, 


ani = animation.FuncAnimation(fig, animate, np.arange(0, sol.shape[1]), interval = 5, blit=False)

f = r"Cantilever_beam.gif" 
writergif = animation.PillowWriter(fps=30) 
ani.save(f, writer=writergif)

### Main loop
root.mainloop()
#%%
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +     1D Supported Beam movies of eigenmodes and superpositions     +
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# %% Problem characteristics

# Material properties (constant)
E  = 1       # [N/mm2]
I  = 1       # [mm4]
k  = 1       # [N]
L  = 1       # [m]
nN = 10      # [-]
mu = 1       # [kg/m]

### Boundaries
BCtype = 'fixed'
BC     = (0, 0, 0, 0)   # (w(0), w(L), M(0), M(L))

### Mesh
grid = np.linspace(0, L, nN)

### Load
node   = np.array([1])   # ID of nodes where force is applied
force  = np.array([-k])  # Applied nodal forces
RHS = getPointForce(grid, node, force)

### Main solver
S, M  = getMatrices(grid, E, I, mu, quadrature = True)

e0 = np.zeros(nN*2);    e0[0]  = 1.0
eL = np.zeros(nN*2);    eL[-2] = 1.0

d0 = np.zeros(nN*2);    d0[1]  = 1.0
dL = np.zeros(nN*2);    dL[-1] = 1.0

# Apply BCs
Me, Se, RHSe = fixBeam(M, S, RHS, (e0, eL), (d0, dL), BC, BCtype)

# Solve
steadySol     = sparse.linalg.spsolve(Se, RHSe)

#Simulating superpositions of eigenvectors
modes = np.array([2,3]) #The mode numbers that will be in the superpositions

t_0 = 0
t_f = 1
Nt = 30*10

w_0 = steadySol
w_diff_0 = np.zeros(w_0.shape)
#w_0 = np.ones(w_0.shape)
superposition_dynamic = eigenvalue_method_dynamic(t_0,t_f,Nt,w_0,w_diff_0,Me,Se,modes,Fourier = False)
    
sol = superposition_dynamic
ylim = (-25, 25)

### Figure object definition
fig = plt.Figure(figsize = (5, 3), dpi = 150)

### Tkinter window
root   = Tk.Tk()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(column=0,row=1)

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size' : 9})

### Animation definition
ax     = fig.add_subplot(111)

nData  = 100

nN     = len(grid)
x_plot = np.linspace(grid.min(), grid.max(), nData) 

beam = get_sol(grid, sol[0:-2, 0])
line, = ax.plot(x_plot, beam(x_plot) * 1e3, color= '#808080', label = 'numerical')

ax.plot(x_plot, np.zeros(x_plot.shape), color= '#959595', linestyle = '--')
ax.axvline(x = 0, color="black", linestyle="-", linewidth = 5)

ax.set_ylabel('deformation (mm)')
ax.set_xlabel('x-dimension (-)')

ax.tick_params(direction= 'in', which= 'major', length= 4, bottom= True,
    top=True, right= False, left=True, width = 1)

ax.set_ylim(ylim[0], ylim[1])

def animate(i):
    beam = get_sol(grid, sol[0:-2, i])
    line.set_ydata( beam(x_plot) * 1e3)  # update the data
    return line, 


ani = animation.FuncAnimation(fig, animate, np.arange(0, sol.shape[1]), interval = 5, blit=False)

f = r"Supported_beam.gif" 
writergif = animation.PillowWriter(fps=30) 
ani.save(f, writer=writergif)

### Main loop
root.mainloop()